//
// Created by hilto on 27/01/2021.
//

#include "HalfEdge.h"
#include "Edge.h"

HalfEdge::HalfEdge() {

}
HalfEdge::HalfEdge(int cInc[2], int cId, Vertex* p1, Vertex* p2, Edge* cEdge)
{
  m_inc[0] = cInc[0], m_inc[1] = cInc[1];
  m_id = cId;
  m_edge = cEdge;
  m_p1 = p1, m_p2 = p2;
}

HalfEdge* HalfEdge::getTwin()
{
  return m_edge->getTwin(this->m_id);
}



//HalfEdge::HalfEdge(std::vector<int> cInc, int cId)
//{
//    m_inc = cInc;
//    m_id = cId;
//    m_edgeId = cEdgeId;
//    m_heNext = heNext;
//    m_elId = cElId;
//}

