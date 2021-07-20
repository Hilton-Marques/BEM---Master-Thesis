//
// Created by hilto on 27/01/2021.
//

#ifndef MAIN_HALFEDGE_H
#define MAIN_HALFEDGE_H

#include <array>
class Face; // forward declaration
class Edge; // forward declaration
class Vertex; // forward declaration
class HalfEdge {
protected:
public:

    int m_inc[2];
    int m_id = -1;
    Edge* m_edge = nullptr;
    HalfEdge* m_heNext;
    Face* m_el ;
    Vertex* m_p2;
    Vertex* m_p1;

    HalfEdge();
   // HalfEdge(std::vector<int> cInc, int cId, int cEdgeId, int cElId, HalfEdge* heNext);
    HalfEdge(int cInc[2], int cId, Vertex* p1, Vertex* p2, Edge* cEdge = nullptr);
    Vertex* getP0() { return m_p1; };
    Vertex* getP1() { return m_p2; };
    HalfEdge* getTwin();
};


#endif //MAIN_HALFEDGE_H
