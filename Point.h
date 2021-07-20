//
// Created by hilto on 27/01/2021.
//
#ifndef MAIN_POINT_H
#define MAIN_POINT_H

#include <vector>


class Vertex;
class Point {
public:
    double m_x;
    double m_y;
    double m_z;
    
public:
    Point();
    Point(double _cX, double _cY, double _cZ);
    // Addition
    friend inline Point operator+(Point p1,Point p2) {
        return Point(p1.m_x + p2.m_x, p1.m_y + p2.m_y, p1.m_z + p2.m_z);
    };
    // Division
    friend inline Point operator/(Point p, double a) {
        if (a == 0.0)
          return Point(0.0, 0.0, 0.0);

        return Point(p.m_x / a, p.m_y / a, p.m_z / a);
    }
    friend inline Point operator*(const double& a, const Point& p) {
        return Point(a * p.m_x, a * p.m_y, a * p.m_z);
    }
    friend inline double operator*(const Point& p1, const Point& p2) {
        return double (p1.m_x*p2.m_x + p1.m_y * p2.m_y + p1.m_z * p2.m_z);
    }
    friend inline Point operator-(const Point& p1,const Point& p2) {
        return Point(p1.m_x - p2.m_x, p1.m_y - p2.m_y, p1.m_z - p2.m_z);
    }

    inline double getNorm() const {
        return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
    }

    inline double operator[](int& id) {
        switch (id)
        {
        case 0:
            return m_x;
        break;
        case 1:
            return m_y;
        break;
        case 2:
            return m_z;
        break;
        }
    }
    static Point crossProd(const Point& v1, const Point& v2);
    void normalize();

};



#endif //MAIN_POINT_H
