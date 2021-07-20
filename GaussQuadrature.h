#ifndef GAUSS_QUADRATURE
#define GAUSS_QUADRATURE	

#include <vector>
class GaussQuad
{
public:
	GaussQuad(int N);
	GaussQuad();
	~GaussQuad();
	void init(const int &N);

	int m_N;
	int m_totN;
	double* m_gaussPointsXi;
	double* m_gaussPointsEta;
	double* m_gaussWeights;
	



private:

};


#endif // !GAUSS_QUADRATURE
