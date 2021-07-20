//
// Created by hilto on 27/01/2021.
//

#include "Point.h"

Point::Point()
: m_x(0.0), m_y(0.0)
{
}

Point::Point(double _cX, double _cY, double _cZ)
: m_x(_cX), m_y(_cY), m_z(_cZ)
{
}


Point Point::crossProd(const Point& v1, const Point& v2)
{
	return Point((v1.m_y*v2.m_z - v1.m_z*v2.m_y), (v1.m_z * v2.m_x - v1.m_x * v2.m_z), (v1.m_x * v2.m_y - v1.m_y * v2.m_x));
}

void Point::normalize()
{
	double norm = getNorm();
	if (norm > 0 )
	{
		m_x /= norm;
		m_y /= norm;
		m_z /= norm;
	}

}

