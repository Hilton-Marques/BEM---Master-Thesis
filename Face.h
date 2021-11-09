//
// Created by hilto on 27/01/2021.
//

#ifndef MAIN_FACE_H
#define MAIN_FACE_H



#include <vector>
#include <string>
#include <complex>

#include "Vertex.h"
#include "HalfEdge.h"
#include "Edge.h"
#include "matrix.h"

//#include "Solid.h"

class Face {

public:
    int m_hedInc;
    int m_id;
    double m_q[3] = { NAN, NAN, NAN };
    std::string m_typeEl;
    std::vector <HalfEdge*> m_heds;
    Vertex* m_points[3];
    bool markTemp = false;
    std::complex<double>* m_MEH;
    
    std::complex<double>* m_MEG;

    std::vector<std::complex<double>* > m_MEG_p;

    std::vector<std::complex<double>* > m_MEH_p;

    bool mark_concave = false;


    Point m_normal;
    double m_jacobian; 
    Face* m_adjacentElementsEdges[3];
    Face* m_fatherElement;
    Face* m_childrenElement[4] = { nullptr, nullptr, nullptr, nullptr };
    int m_N;
    std::vector<Vertex*> m_collectVertexes;
    bool m_isLeaf = false;
    std::vector<Face*> m_adjacentFaces;
    std::vector<Face*> m_concave_adjacency;

    bool m_bd = false;
    short int m_nL;
    int m_nCloseNodes;
    Matrix m_H{ m_nCloseNodes,3 };
    Matrix m_G{ m_nCloseNodes,3 };

    Face();
    Face(int cHedInc ,double q[], short int level);
    std::vector<Face*> getAdjacentElements(bool flag = false);
    Point getNormalUnitary();
    void setNormal();
    void setFMMSizes(const int& N);
    Point getYc();
    void getLeafElements(std::vector<Face*>* leafElements);
    void getLevelElements(std::vector<Face*>* leafElements, short int level);
    void accumulate();
    void reset();
    std::vector<Vertex*> getBoundaryPoints();
    void showMEG();
    void showMEH();
    void set_FMM_ME_sizes();
    Matrix exactIntG(Vertex* source);
    Matrix exactIntH(Vertex* source);
    std::vector<Face*> getChildrenElement();
    double getRadius(double fac);
    void buildConcaveAdjacency(std::vector<Face*> & range_elements, double fac);

    ~Face();

    



};


#endif //MAIN_FACE_H
