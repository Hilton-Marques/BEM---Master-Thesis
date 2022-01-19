#include "DunavantQuadrature.h"

DunavantQuad::DunavantQuad()
{
}

DunavantQuad::DunavantQuad(int p)
{
	init(p);
}

DunavantQuad::~DunavantQuad()
{
}

void DunavantQuad::init(int p)
{
	
	m_totN = m_n[p];
}
