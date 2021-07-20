#ifndef GMRES_H
#define GMRES_H

#include "matrix.h"
#include "FMM.h"
class FMM;
class GMRES
{
public:
	GMRES(FMM* fmm);
	FMM* m_fmm;
	double m_resid;
	void Update(Matrix& x, int k, Matrix& h, Matrix& s, Matrix& Q);
	Matrix Solver(Matrix& b,  int m = 10,  int max_iter = 10,  double tol = 0.000001);
	void GeneratePlaneRotation(double &dx, double &dy, double &cs, double & sn);
	void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn);
	Matrix SolveM(Matrix &A, Matrix& b, int m = 10, int max_iter = 10, double tol = 0.0001);
	Matrix BiCG(Matrix& b, int m = 10, int max_iter = 10, double tol = 0.0001);
	Matrix BiCGM(Matrix &A , Matrix& b, int m = 30, int max_iter = 30, double tol = 0.0001);

};


#endif // !GMRES_H

