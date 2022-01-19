#include <vector>
#include <string>
#include <array>
#include <cmath>
#include "Solid.h"
#include "Point.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "FMM.h"
#include "RTable.h"
#include "input.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> // set_difference
#include <set> // erase duplicates



extern const double M_PI = atan(1) * 4;
extern const double Const = 1 / (4 * M_PI) ;

void buildMotherMesh(std::vector<Vertex*>& vertices, 
                    std::vector<Edge*>& edges, 
                    std::vector<Face*>& elements,
                    std::vector<HalfEdge*>& heds, 
                    const std::vector<Point>& coords, 
                    std::vector<std::vector<int>>& inc)
{
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = new Vertex(coords[i], i);
  }

  for (int i = 0; i < elements.size(); i++)
  {
    std::vector<HalfEdge*> hedsFace(3);
    for (int j = 0; j < 3; j++)
    {
      int inc_i[2] = { inc[i][j], inc[i][(j + 1) % 3] };
      Vertex* p0 = vertices[inc_i[0]];
      Vertex* p1 = vertices[inc_i[1]];
      HalfEdge* hedi = new HalfEdge(inc_i, heds.size(), p0, p1);
      heds.push_back(hedi);
      p0->m_star.push_back(p1);
      p0->m_hed_star.push_back(hedi);
      hedsFace[j] = hedi;
    }
    double q[3] = { NAN, NAN, NAN };
    Face* facei = new Face(hedsFace[0]->m_id, q, 1);
    facei->m_heds = hedsFace;
    for (int j = 0; j < 3; j++)
    {
      hedsFace[j]->m_heNext = hedsFace[(j + 1) % 3];
      hedsFace[j]->m_el = facei;
      facei->m_points[j] = hedsFace[j]->getP0();
    }
    elements[i] = facei;
  }
  for (Vertex* pt : vertices)
  {
    for (int i = 0; i < pt->m_star.size(); i++)
    {
      Vertex* pt_star_i = pt->m_star[i];
      HalfEdge* hed_star_i = pt->m_hed_star[i];
      if (hed_star_i->m_edge == nullptr)
      {
        for (int j = 0; j < pt_star_i->m_star.size(); j++)
        {
          Vertex* pt_star_j = pt_star_i->m_star[j];
          HalfEdge* hed_star_j = pt_star_i->m_hed_star[j];
          if (pt_star_j == pt)
          {
            Edge* edgei = new Edge(hed_star_i, hed_star_j, heds.size());
            edges.push_back(edgei);
            hed_star_i->m_edge = edgei;
            hed_star_j->m_edge = edgei;
          }
        }
      }
    }
  }
}

void BuildConcaveAdjacency(std::vector<Face*>& elements, double fac)
{
    for (Face* element : elements)
    {
        std::vector<Face*> concaveAdjacency;
        //std::vector<Face*> convexAdjacency = element->getAdjacentElements();
        //std::sort(convexAdjacency.begin(), convexAdjacency.end());
        //std::set_difference(elements.begin(), elements.end(), convexAdjacency.begin(), convexAdjacency.end(), std::inserter(posible_concaveAdjacency, posible_concaveAdjacency.begin()));
        
        //element->buildConcaveAdjacency(elements,fac);

        double radius2 = element->getRadius(fac);
        for (int i = 0; i < 3; i++)
        {
            Point center = element->m_points[i]->m_coord;
            for (Face* element_j : elements)
            {
                Point d = element_j->getYc() - center;
                if ( (d * d) < radius2 )
                {
                    concaveAdjacency.push_back(element_j);
                }
            }
        }
        // erase duplicates
        std::set<Face*> s(concaveAdjacency.begin(), concaveAdjacency.end());
        concaveAdjacency.assign(s.begin(), s.end());
        element->m_concave_adjacency = concaveAdjacency;

    }
}

int main() {
    std::string fileName = "entradaToro.txt";
    std::vector<Point> coords;
    std::vector<std::vector<int>> inc;
    ReadInput input(fileName, coords, inc);
    const int nPts = coords.size();
    const int nEl = inc.size();
    const int nEdges = nPts + nEl - 2;
    const int nHeds = 2 * nEdges;
    std::ofstream fw("output.txt", std::ofstream::out);
    std::vector<Vertex*> vertexs(nPts);
    std::vector<Face*> elements(nEl);
    std::vector<Edge*> edges;
    std::vector<HalfEdge*> heds;
    std::vector<Vertex*> vertexsBd(1);
    double concave_fac_radius = 1.0;
    buildMotherMesh(vertexs, edges, elements, heds, coords, inc);
    BuildConcaveAdjacency(elements, concave_fac_radius);
    // start FMM
    int nL = 7;
    bool concave_allow = true;
    
    Solid solid(vertexs, edges, heds, elements, nL, vertexsBd, concave_allow, concave_fac_radius);
    // Level

    int N = 10; // truncation term
    int NG = 11; // gauss quadrature
    int Nc = 2; // level of adjacency

    FMM fmm( & solid, N, NG, Nc, 1,&fw);

    return 0;


}

//int main() {
//  std::string fileName = "entradaToro.txt";
//  std::vector<Point> coords;
//  std::vector<std::vector<int>> inc;
//  ReadInput input(fileName, coords, inc);
//  const int nPts = coords.size();
//  const int nEl = inc.size();
//  const int nEdges = nPts + nEl - 2;
//  const int nHeds = 2* nEdges;
//  //open file for writing
//  std::ofstream fw("output.txt", std::ofstream::out);
//int NG = 12;
//
//std::vector<std::string> fields = { std::string("linear"), std::string("quadratico"), std::string("cubic") };
//
//for (int field = 0; field < 3; field++)
//{
//    std::cout << "field" << fields[field] << "\n";
//    fw << "field" << fields[field] << "\n";
//    fw.precision(15);
//    for (int Nc = 4; Nc < 5; Nc++)
//    {
//        std::cout << "Nc  " << Nc << "\n";
//        fw << "Nc  " << Nc << "\n";
//
//        int nLmin = 2 + Nc - 1;
//        int nL = nLmin;
//        for (int j = nL; j < 8; j++)
//        {
//            nL++;
//            nL = 8;
//            std::cout << "nivel" << nL << "\n";
//            fw << "nivel" << nL << "\n";
//
//            for (int i = 2; i < 11; i = i + 2)
//            {
//                std::vector<Vertex*> vertexs(nPts);
//                std::vector<Face*> elements(nEl);
//                std::vector<Edge*> edges;
//                std::vector<HalfEdge*> heds;
//                std::vector<Vertex*> vertexsBd(1);
//                double concave_fac_radius = 1.0;
//                buildMotherMesh(vertexs, edges, elements, heds, coords, inc);
//                BuildConcaveAdjacency(elements, concave_fac_radius);
//                // start FMM
//                bool concave_allow = false;
//
//                Solid solid(vertexs, edges, heds, elements, nL, vertexsBd, concave_allow, concave_fac_radius);
//
//                int N = i;
//                FMM fmm(&solid, N, NG, Nc, field, &fw);
//
//
//            }
//        }
//    }
//}
//    //fmm.showB();
//     //fmm.showX();
//
//    return 0;
//}
//
