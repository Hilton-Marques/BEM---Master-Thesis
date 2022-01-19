#ifndef FMM_H
#define FMM_H

#include <vector>
#include "Face.h"
#include "Solid.h"
#include "Solver.h"
#include "Solid.h"
#include "gmres.h"
#include "matrix.h"

#include <iostream>
#include <stack>
#include <ctime>
#include <fstream>

class FMM
{
public:

	FMM();
	~FMM();
	FMM(Solid* solid, int n, int ng, int nc, int field, std::ofstream* fw);
	std::vector<Vertex*> getLeafNodesFromElements(std::vector<Face*> adjacentElements);
	std::vector<Face*> getFarElements(Face* element);
	std::vector< std::vector<Face*> > *m_elementsLevel;
	void checkWithConventionalBem(Solid* solid, std::vector<Face*> elements, Matrix &fmmAx);
	void checkWithConventionalBemH(Solid* solid, std::vector<Face*> elements);
	Matrix computeBvector(double & erro);
	Matrix matrixVectorMulti(Matrix &x);
	Matrix getHarmonicTemperature();
	std::stack<clock_t> m_tictoc_stack;
	int m_nLMax;
	int m_nTotalVerts;
	int m_N;
	int m_nDb;
	int m_tot;
	int m_Nc;
	int m_field;
	double m_resid;
	Solver* m_solver;
	Solid* m_solid;
	Matrix m_b{ m_nTotalVerts , 1 };
	Matrix m_x{ m_nTotalVerts , 1 };
	Matrix m_Hdb{ m_nTotalVerts , m_nDb }; // prescribed displacements
	Matrix m_db{ m_nDb , 0 }; // prescribed displacements
	void prepareForGMRES();
	void showB();
	void showX();
	void setFMM_ME_sizes();
	double trivialField(double x, double y, double z);
	double trivialGrad(double x, double y, double z, Point& n);
	void changeBoundaryCondition(std::vector<Vertex*> points, std::vector<Face*> elements);
	void tic() { m_tictoc_stack.push(clock()); };
	void toc();
	Matrix Gt(Solid* solid, std::vector<Face*> elements);
	Matrix Hd(Solid* solid, std::vector<Face*> elements);
	std::ofstream* m_fw;
	};



#endif // !FMM_H
