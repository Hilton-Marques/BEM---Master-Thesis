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

int main() {

#pragma omp parallel
    {
        printf("Hello, world.\n");
    }
  std::string fileName = "entradaTetra.txt";
  std::vector<Point> coords;
  std::vector<std::vector<int>> inc;
  ReadInput input(fileName, coords, inc);
  const int nPts = coords.size();
  const int nEl = inc.size();
  const int nEdges = nPts + nEl - 2;
  const int nHeds = 2* nEdges;


    //// Domain
    //const int nHeds = 12;
    //const int nPts = 4;
    //const int nEdges = 6;
    //const int nEl = 4;
    //const int nL = 7;

    //// Initialize

    //std::vector<Vertex*> vertexs(4);
    //std::vector<Edge*> edges(6);
    //std::vector<HalfEdge*> heds(12);
    ////std::vector<Vertex*> vertexsBd(1);

    //// Boundary condiition per face
    //double q1[3] = { 0.0, 0.0, 0.0 };
    //double q2[3] = { 0.0, 0.0, 0.0 };
    //double q3[3] = { 1.366972252, 1.366972252,1.366972252 };
    //double q4[3] = { -0.936329, -0.936329,-0.936329 };


    //// Type of element
    //std::string typeEl = "T3";

    //// Coordinates of piramid
    //Point p1 = Point(-3.0,-1.5,0);
    //Point p2 = Point(3.0 , -1.5, 0);
    //Point p3 = Point(0 , 1.5, 0);
    //Point p4 = Point(0, 0, 4);

    //// Create Edges
    //
    //edges[0] = new Edge(0);
    //edges[1] = new Edge(1);
    //edges[2] = new Edge(2);
    //edges[3] = new Edge(3);
    //edges[4] = new Edge(4);
    //edges[5] = new Edge(5);

    //// Create Hed
    //

    //vertexs[0] = new Vertex (p1, 0);
    //vertexs[1] = new Vertex (p2, 1);
    //vertexs[2] = new Vertex (p3, 2);
    //vertexs[3] = new Vertex (p4, 3);   
    //vertexs[0]->m_bd = true;
    //vertexs[0]->m_u = -4.5;
    //// vertexs with boundary
    //vertexsBd[0] = vertexs[0];
    //// Heds

    //int inc1[2] = { 0,1 };
    //int inc2[2] = { 1,2 };
    //int inc3[2] = { 2,0 };
    //int inc4[2] = { 3,0 };
    //int inc5[2] = { 0,2 };
    //int inc6[2] = { 2,3 };
    //int inc7[2] = { 2,1 };
    //int inc8[2] = { 1,3 };
    //int inc9[2] = { 3,2 };
    //int inc10[2] = { 1,0 };
    //int inc11[2] = { 0,3 };
    //int inc12[2] = { 3,1 };

    //heds[0] = new HalfEdge(inc1, 0, vertexs[inc1[0]], vertexs[inc1[1]], edges[0]);
    //heds[1] = new HalfEdge(inc2, 1, vertexs[inc2[0]], vertexs[inc2[1]], edges[1]);
    //heds[2] = new HalfEdge(inc3, 2, vertexs[inc3[0]], vertexs[inc3[1]], edges[2]);
    //heds[3] = new HalfEdge(inc4, 3, vertexs[inc4[0]], vertexs[inc4[1]], edges[3]);
    //heds[4] = new HalfEdge(inc5, 4, vertexs[inc5[0]], vertexs[inc5[1]], edges[2]);
    //heds[5] = new HalfEdge(inc6, 5, vertexs[inc6[0]], vertexs[inc6[1]], edges[5]);
    //heds[6] = new HalfEdge(inc7, 6, vertexs[inc7[0]], vertexs[inc7[1]], edges[1]);

    //heds[7] = new HalfEdge( inc8,  7, vertexs[inc8[0]], vertexs[inc8[1]], edges[4]);
    //heds[8] = new HalfEdge(inc9, 8, vertexs[inc9[0]], vertexs[inc9[1]], edges[5]);
    //heds[9] = new HalfEdge(inc10, 9, vertexs[inc10[0]], vertexs[inc10[1]], edges[0]);
    //heds[10] = new HalfEdge(inc11, 10, vertexs[inc11[0]], vertexs[inc11[1]], edges[3]);
    //heds[11] = new HalfEdge(inc12, 11, vertexs[inc12[0]], vertexs[inc12[1]], edges[4]);

    ////Update Edges;
    //edges[0]->m_hed1 = heds[0];
    //edges[0]->m_hed2 = heds[9];

    //edges[1]->m_hed1 = heds[1];
    //edges[1]->m_hed2 = heds[8];

    //edges[2]->m_hed1 = heds[2];
    //edges[2]->m_hed2 = heds[4];

    //edges[3]->m_hed1 = heds[3];
    //edges[3]->m_hed2 = heds[10];

    //edges[4]->m_hed1 = heds[7];
    //edges[4]->m_hed2 = heds[11];

    //edges[5]->m_hed1 = heds[5];
    //edges[5]->m_hed2 = heds[8];

    //
    //// Create Elements
    //std::vector<Face*> elements(4);

    //elements[0] = new Face(0,q1, 1);
    //elements[0]->m_bd = true;
    //elements[1] = new Face(3,q2, 1);
    //elements[1]->m_bd = true;
    //elements[2] = new Face(6,q3, 1);
    //elements[2]->m_bd = true;
    //elements[3] = new Face(9,q4, 1);
    //elements[3]->m_bd = true;

    //elements[0]->m_adjacentFaces = std::vector<Face*>{ elements[1], elements[2], elements[3] };
    //elements[1]->m_adjacentFaces = std::vector<Face*>{ elements[0], elements[2], elements[3] };
    //elements[2]->m_adjacentFaces = std::vector<Face*>{ elements[0], elements[1], elements[3] };
    //elements[3]->m_adjacentFaces = std::vector<Face*>{ elements[0], elements[1], elements[2]};
    //
    //int count = 0;
    //for (int i = 0; i < elements.size(); i++)
    //{
    //    elements[i]->m_heds[0] = heds[count];
    //    count++;
    //    elements[i]->m_heds[1] = heds[count];
    //    count++;
    //    elements[i]->m_heds[2] = heds[count];
    //    count++;
    //}

    


int NG = 10;

std::vector<std::string> fields = { std::string("linear"), std::string("quadratico"), std::string("cubic") };

//for (int field = 0; field < 3; field++)
//{
//    std::cout << "field" << fields[field] << "\n";
//    for (int nL = 3; nL < 8; nL++)
//    {
//        std::vector<Vertex*> vertexs(nPts);
//        std::vector<Face*> elements(nEl);
//        std::vector<Edge*> edges;
//        std::vector<HalfEdge*> heds;
//        std::vector<Vertex*> vertexsBd(1);
//        buildMotherMesh(vertexs, edges, elements, heds, coords, inc);
//        // start FMM
//
//        Solid solid(vertexs, edges, heds, elements, nL, vertexsBd);
//        // Level
//
//        FMM fmm(&solid, 2, NG, 1, field);
//
//    }
//}
for (int field = 1; field < 3; field++)
{
    std::cout << "field" << fields[field] << "\n";

    for (int Nc = 1; Nc < 4; Nc++)
    {
        std::cout << "Nc  " << Nc << "\n";
        int nLmin = 2 + Nc - 1;
        int nL = nLmin;
        for (int j = nL; j < 10 - nLmin; j++)
        {
            nL++;
            std::cout << "nivel" << nL << "\n";

            for (int i = 0; i < 9; i = i + 2)
            {
                std::vector<Vertex*> vertexs(nPts);
                std::vector<Face*> elements(nEl);
                std::vector<Edge*> edges;
                std::vector<HalfEdge*> heds;
                std::vector<Vertex*> vertexsBd(1);
                buildMotherMesh(vertexs, edges, elements, heds, coords, inc);
                // start FMM

                Solid solid(vertexs, edges, heds, elements, nL, vertexsBd);
                // Level

                int N = 2 + i;
                FMM fmm(&solid, N, NG, Nc, field);


            }
        }
    }
}
    //fmm.showB();
     //fmm.showX();

    return 0;
}

