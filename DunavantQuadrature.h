#ifndef DUNAVANT_QUADRATURE
#define DUNAVANT_QUADRATURE

class DunavantQuad
{
public:
	DunavantQuad();
	~DunavantQuad();



private:
	// the following matrixes correspond to 20x79 where 20 is the hightest degree 
	double m_weight[79][20];
	double m_xi[79][20];
	double m_ksi[79][20];
	int m_n[20];
};



#endif
