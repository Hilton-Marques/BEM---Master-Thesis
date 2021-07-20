//
// Created by hilto on 27/01/2021.
//

#include "Vertex.h"
#include "Face.h"
#include <iostream>
extern const double M_PI;
Vertex::Vertex() {

}

Vertex::Vertex(Point c_point, int c_id)
{
    m_coord = c_point;
    m_id = c_id;
}

double Vertex::calculateSolidAngle()
{
  double fi = getDihedralSum();

  // just for normal inside the surface (the value should negativate for another orientatnion)
  // Or u could subtract by 4* Pi and get the integer part
  
  return -(2 * M_PI - fi );


}

double Vertex::getDihedralSum()
{
  // getFirstHed
  HalfEdge* hed0 = nullptr;
  for (HalfEdge* hed : m_elements[0]->m_heds)
  {
    if (hed->getP0() == this)
    {
      hed0 = hed;
      break;
    }
  }
  if (hed0 == nullptr)
    return 0.0;

  Vertex* v0_star = hed0->getP1();  
  HalfEdge* hed_i = hed0;
  Vertex* v1_star = nullptr;
  Vertex* v0_star_bef = v0_star;
  double fi = 0.0; // dihedral angle
  while (v1_star != v0_star)
  {
    Face* element_i = hed_i->m_el;
    Point n_i = element_i->getNormalUnitary();
    HalfEdge* hed_i_twin = hed_i->getTwin();
    Face* element_j = hed_i_twin->m_el;
    Point n_j = element_j->getNormalUnitary();
    Point edge = v0_star_bef->m_coord - m_coord;
    edge.normalize();
    fi += std::atan2(Point::crossProd(n_i, n_j) * edge, n_i * n_j); // compute dihedral
    hed_i = hed_i_twin->m_heNext;
    v1_star = hed_i->getP1();
    v0_star_bef = v1_star;
  }


    return fi;
}





