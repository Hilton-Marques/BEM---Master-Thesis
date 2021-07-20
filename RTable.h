#ifndef RTABLE_H
#define RTABLE_H

#include <complex>
#include <vector>
#include "Point.h"

class RTable
{
public:

	RTable();
	~RTable();
	int m_N;
	int m_NR;
	int m_totNR;
	int m_NS;
	int m_totNS;
	std::vector<std::complex<double>> m_Ry;
	std::complex<double>* m_RecRy;
	std::complex<double>* m_RecSx;
	std::complex<double>* m_RecRyFull;
	std::complex<double>* m_RecSxFull;


	//std::vector<std::complex<double>> m_Ry;
	

	//std::complex<double>* m_Ry;


	void setSizeRTable(const int& N);
	void setSizeRecTableR(const int& N);
	void setSizeRecTableS(const int &N);
	std::vector<std::complex<double>> evaluateTableR(const Point& y);
	std::vector<std::complex<double>> evaluateRecursiveTableR(const Point& y );
	std::vector<std::complex<double>>  evaluateRecursiveTableS(const Point& x);
	std::complex<double> getPos(const int & n, const int &m );
	std::complex<double> getPosCdR(const int n, const int m, std::vector<std::complex<double>>& R);
	std::complex<double> getProjdR(const int n, const int m, std::vector<std::complex<double>> &R);
	int getPosZdR(const int n, const int m);
	static const int getIdNM(const int& N, const int& M);
	

private:

};

#endif // !RTABLE_H
