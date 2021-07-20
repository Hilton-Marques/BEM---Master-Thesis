//
// Created by hilto on 27/01/2021.
//

#include "Solid.h"
#include "HalfEdge.h"
#include <cmath> // pow


Solid::Solid(std::vector<Vertex*> verts, std::vector<Edge*> edges, std::vector<HalfEdge*> heds, std::vector<Face*> elements,
             int nL, std::vector<Vertex*> vertsBd):
m_verts(verts), m_edges (edges), m_heds(heds), m_elementsMother(elements) , m_nL(nL), m_vertsBd(vertsBd)

{
    m_nHeds = m_heds.size();
    m_nPts = m_verts.size();
    m_nEdges = m_edges.size();
    m_nEl   = m_elementsMother.size();
    m_elementsLevel.reserve(m_nL);
    const int nTotalEl = m_nEl * (std::pow(4, m_nL) - 1) / 3; // sum of P.G.
    const int maxNEdges = calculateNMaxEdges(m_nL, m_nEdges, m_nEl);
    const int maxNHeds = 2 * maxNEdges;
    const int maxNVerts = calculateNMaxVerts(m_nL, m_nEdges, m_nEl, m_nPts);
    
    //Initialize
    m_elements.reserve(nTotalEl);
    m_verts.reserve(maxNVerts);
    m_edges.reserve(maxNHeds);
    m_heds.reserve(maxNHeds);
    Start();
}

int Solid::calculateNMaxVerts(const int nL, int nEdges, const int nEl, const int nVerts)
{
    int m_nL = 2;
    int maxNVerts = nVerts + nEdges;
    while (m_nL < nL)
    {
        nEdges = nEdges * 2 + nEl * std::pow(4, m_nL - 2) * 3;
        maxNVerts = maxNVerts + nEdges;
        m_nL++;
    }
    return maxNVerts;
}

int Solid::calculateNMaxEdges(const int nL, int nEdges, const int nFaces)
{
    int m_nL = 1;
    int maxEdges = nEdges;
    while (m_nL < nL)
    {
        nEdges = nEdges * 2 + nFaces * std::pow(4, m_nL - 1) * 3;
        maxEdges = maxEdges + nEdges;
        m_nL++;
    }
    return maxEdges;
}

Solid::~Solid()
{
    for (Edge* edge : m_edges) {
        delete edge;
    }
    for (HalfEdge* hed : m_heds)
    {
        delete hed;
    }
    for (Vertex* vert : m_verts)
    {
        delete vert;
    }
    for (int i = 0; i < m_elements.size(); i++)
    {
      //Face* face = m_elements[i];
      //delete face;
    }
}

void Solid::Start() {
    m_elementsLevel.push_back(m_elementsMother);
    m_elements.insert(m_elements.end(), std::begin(m_elementsMother), std::end(m_elementsMother));
    
    for (int i = 1; i < m_nL; i++)
    {
        std::vector<Face*> elements = m_elementsLevel[i - 1];
        int nElNv = elements.size();
        std::vector<Face*> newElements;
        newElements.reserve(4 * nElNv);
        
        for (Face* element : elements)
        {
            
            std::vector<HalfEdge*> localHeds;
            std::vector<Vertex*> localPoints;
            localHeds.reserve(6);
            localPoints.reserve(3);

            for (HalfEdge* hed_i : element->m_heds)
            {
              
             
                Edge* edge_i = hed_i->m_edge;

                if (!edge_i->m_isSplited) {

                    // create midPoint on edge
                    Vertex* p1 = m_verts[hed_i->m_inc[0]];
                    Vertex* p2 = m_verts[hed_i->m_inc[1]];
                   
                    Point newCoord = (p1->m_coord + p2->m_coord) / 2.0; // Mean

                    // create new vertex

                    m_nPts++;
                    //Vertex newVertex = Vertex(newCoord, m_nPts-1);
                    m_verts.push_back(new Vertex(newCoord, m_nPts - 1));
                    Vertex* newVertex = m_verts.back();
                    localPoints.push_back(newVertex);

                    // Create new heds(4 new heds
                    for (int k = 0; k < 2; k++)
                    {
                        m_nEdges++;
                        m_edges.push_back(new Edge(m_nEdges - 1));
                        Edge* newEdge = m_edges.back();
                        edge_i->m_childrenId[k] = newEdge;

                        int inc1[2] = { hed_i->m_inc[k],m_nPts -1 };
                        int inc2[2] = { m_nPts-1 , hed_i->m_inc[k] };
                        int incs[2][2] = { {inc1[0],inc1[1]} , {inc2[0],inc2[1]} };

                        //first hed
                        m_nHeds++;
                        m_heds.push_back(new HalfEdge(incs[k], m_nHeds - 1, m_verts[incs[k][0]], m_verts[incs[k][1]],  newEdge));
                        HalfEdge* newHed = m_heds.back();
                        //HalfEdge newHed = HalfEdge(incs[k], m_nHeds -1 , m_nEdges -1 );

                        // twin
                        m_nHeds++;
                        m_heds.push_back(new HalfEdge(incs[(k + 1) % 2], m_nHeds - 1, m_verts[incs[(k + 1) % 2][0]], m_verts[incs[(k + 1) % 2][1]] , newEdge));
                        HalfEdge* newHedTwin = m_heds.back();
                        
                        
                       
                        // update

                        localHeds.push_back(newHed);
                        newEdge->m_hed1 = newHed;
                        newEdge->m_hed2 = newHedTwin;
                        
                    }

                    edge_i->m_isSplited = true;
                    edge_i->m_mid = newVertex;

                }
                else {
                    localPoints.push_back(edge_i->m_mid);
                    for (int k = 1; k > -1; k--)
                    {
                        localHeds.push_back(edge_i->m_childrenId[k]->m_hed2);
                    }
                 }
            }

            std::vector<HalfEdge*> internalHeds(3);
            for (int k = 0; k < 3; k++)
            {
                m_nEdges++;
                m_edges.push_back(new Edge(m_nEdges - 1));
                Edge* newEdge = m_edges.back();

                int p2 = localPoints[(k + 1) % 3]->m_id;
                int p1 = localPoints[k % 3]->m_id;

                //newHed
                m_nHeds++;
                int inc1[2] = { p1, p2 };
                m_heds.push_back(new HalfEdge(inc1, m_nHeds - 1, m_verts[inc1[0]], m_verts[inc1[1]], newEdge));
                HalfEdge* newHed = m_heds.back();
                internalHeds[k] = newHed;

                //HalfEdge newHed = HalfEdge(std::vector<int> {p1, p2}, m_nHeds -1 , m_nEdges -1);
                //newHedTwin
                m_nHeds++;
                int inc2[2] = { p2, p1 };
                m_heds.push_back(new HalfEdge(inc2, m_nHeds - 1, m_verts[inc2[0]], m_verts[inc2[1]], newEdge));
                HalfEdge* newHedTwin = m_heds.back();

                //HalfEdge newHedTwin = HalfEdge(std::vector<int> {p2, p1}, m_nHeds -1 , m_nEdges -1);             
                //newEdge    

                //Edge newEdge = Edge(newHed, newHedTwin, m_nEdges - 1);

                newEdge->m_hed1 = newHed;
                newEdge->m_hed2 = newHedTwin;
                
            }
            // reorder
            internalHeds = { internalHeds[2], internalHeds[0], internalHeds[1] };

            //update hNext for border heds and create new Elements
            int count1 = 0;
            for (int k = 0; k < 5; k = k + 2)
            {
                m_nEl++;                
                int firstId = localHeds[k]->m_id;
                int idMidle = internalHeds[count1]->m_edge->m_hed2->m_id;
                localHeds[k]->m_heNext = m_heds[idMidle];                
                int idBefore = (6 + k - 1) % 6;
                m_heds[idMidle]->m_heNext = localHeds[idBefore];
                
                localHeds[idBefore]->m_heNext = localHeds[k];
                
                // create new Face
                newElements.push_back(new Face(localHeds[k]->m_id, element->m_q, element->m_nL + 1));
                Face* newFace = newElements.back();
                newFace->m_heds[0] = localHeds[k];
                newFace->m_heds[1] = m_heds[idMidle];
                newFace->m_heds[2] = localHeds[idBefore];
                //newFace->m_points = { m_verts[localHeds[k]->m_inc[0]], m_verts[m_heds[idMidle]->m_inc[0]], m_verts[localHeds[idBefore]->m_inc[0]] };
                newFace->m_points[0] = m_verts[localHeds[k]->m_inc[0]];
                newFace->m_points[1] = m_verts[m_heds[idMidle]->m_inc[0]];
                newFace->m_points[2] = m_verts[localHeds[idBefore]->m_inc[0]];
                newFace->m_fatherElement = element;
                //Update new Face
                localHeds[k]->m_el = newFace;
                m_heds[idMidle]->m_el = newFace;
                localHeds[idBefore]->m_el = newFace;
                //Update Element
                element->m_childrenElement[count1] = newFace;

                
                
                count1++;
                

            }
            // update hNext for internal heds and create central Element
            m_nEl++;
            newElements.push_back(new Face(internalHeds[0]->m_id, element->m_q, element->m_nL + 1));
            Face* newFace = newElements.back();
            newFace->m_fatherElement = element;
            element->m_childrenElement[3] = newFace;
            for (int k = 0; k < 3; k++)
            {
                internalHeds[k]->m_heNext = internalHeds[(k + 1) % 3];
                internalHeds[k]->m_el = newFace;
                newFace->m_heds[k] = internalHeds[k];
                newFace->m_points[k] = m_verts[internalHeds[k]->m_inc[0]];
            }  

            
          
        }

        m_elementsLevel.push_back(newElements);
        m_elements.insert(m_elements.end(), std::begin(newElements), std::end(newElements));
    }


}

