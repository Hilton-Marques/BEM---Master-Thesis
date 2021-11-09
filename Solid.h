//
// Created by hilto on 27/01/2021.
//

#ifndef MAIN_SOLID_H
#define MAIN_SOLID_H

#include <vector>
#include "Point.h"

#include "Vertex.h"
#include "Edge.h"
#include "HalfEdge.h"
#include "Face.h"

class Solid {
public:

    Solid(std::vector<Vertex*> verts, std::vector<Edge*> edges, std::vector<HalfEdge*> heds,
          std::vector<Face*> elements, int nL, std::vector<Vertex*> vertsBd, bool concave_allow, double concave_radiu_fac);

    static int calculateNMaxVerts(const int nL, int nEdges, const int nEl, const int nVerts);
    static int calculateNMaxEdges(const int nL, int nEdges, const  int nFaces);
    ~Solid();


    void Start(bool concave_allow, double concave_radiu_fac);


    std::vector<Vertex*> m_verts;
    std::vector<Vertex*> m_vertsBd;
    std::vector<HalfEdge*> m_heds;
    std::vector<Edge*> m_edges;
    std::vector<Face*> m_elementsMother;
    std::vector< std::vector<Face*> > m_elementsLevel;
    std::vector<Face*> m_elements;

    int m_nL;
    int m_nHeds;
    int m_nPts;
    int m_nEdges;
    int m_nEl;



};


#endif //MAIN_SOLID_H
