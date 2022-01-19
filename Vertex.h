//
// Created by hilto on 27/01/2021.
//

#ifndef MAIN_VERTEX_H
#define MAIN_VERTEX_H

#include "Point.h"
#include <vector>

class Face;
class HalfEdge;

class Vertex {

public:
    Vertex();
    Vertex(Point c_point, int c_id);
    double calculateSolidAngle();
    double getDihedralSum();
    Point m_coord;
    int m_id;
    double m_u = 1.0;
    bool markTemp = false;
    bool m_bd = false;
    std::vector <Face*> m_elements;
    std::vector<Face*> m_closeElements;
    std::vector<Vertex*> m_star;
    std::vector<HalfEdge*> m_hed_star;
    std::vector<double> Gt;
};


#endif //MAIN_VERTEX_H
