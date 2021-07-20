#include "RTable.h"
#include <omp.h>
#include <iostream>

RTable::RTable()
{


}

RTable::~RTable()
{
    delete[] m_RecSxFull;
    delete[] m_RecRy;
    delete[] m_RecRyFull;
    delete[] m_RecSx;
    //delete[] m_RecSxFull;

}

std::complex<double> RTable::getPos(const int & n, const int & m)
{

 /*   const  int pos = getIdNM(n, m);
    if (pos < 0)
    {
        return std::complex<double>(0, 0);
    }
    return m_Ry[pos];*/


    if (m >= 0)
    {
        const  int pos = getIdNM(n, m);
        if (pos < 0)
        {
            return std::complex<double>(0, 0);
        }
        return m_Ry[pos];
    }
    else
    {

        const  int pos = getIdNM(n, -m);
        if (pos < 0)
        {
            return std::complex<double>(0, 0);
        }
        if (-m % 2 == 0)
        {
            //return (std::conj(m_Ry[pos]));
            return std::complex<double>(0, 0);

        }
        else
        {
            //return (-std::conj(m_Ry[pos]));
            return std::complex<double>(0, 0);

        }
        

    }
   

   
   
    return std::complex<double> ();
}

void RTable::setSizeRTable(const int& N)
{
    m_N = N;
    m_Ry.resize(28);

}

void RTable::setSizeRecTableR(const int& N)
{
    m_NR = 2*N + 1;
    m_totNR = (m_NR + 2) * (m_NR + 1) * 0.5;
    m_RecRy = new std::complex<double>[m_totNR];
    m_RecRyFull = new std::complex<double>[(m_NR+1)* (m_NR+1)];

}

void RTable::setSizeRecTableS(const int& N)
{
    m_NS = N;
    m_totNS = (N + 2) * (N + 1) * 0.5;
    m_RecSx = new std::complex<double>[m_totNS];
    m_RecSxFull = new std::complex<double>[(N + 1) * (N + 1)];

}

std::vector<std::complex<double>> RTable::evaluateTableR(const Point& y)
{

  double y3 = y.m_z;
  std::complex<double> yb(y.m_x, y.m_y);
  //std::complex<double> conjY(y.m_x, -y.m_y);

  double AbsSquaredY = (y.m_x * y.m_x + y.m_y * y.m_y);


  std::complex<double> powYb2 = yb * yb;
  std::complex<double> powYb3 = yb * powYb2;
  std::complex<double> powYb4 = yb * powYb3;
  std::complex<double> powYb5 = yb * powYb4;
  std::complex<double> powYb6 = yb * powYb5;
  std::complex<double> powYb7 = yb * powYb6;
  std::complex<double> powYb8 = yb * powYb7;
  std::complex<double> powYb9 = yb * powYb8;
  std::complex<double> powYb10 = yb * powYb9;


  double powY32 = y3 * y3;
  double powY33 = y3 * powY32;
  double powY34 = y3 * powY33;
  double powY35 = y3 * powY34;
  double powY36 = y3 * powY35;
  double powY37 = y3 * powY36;
  double powY38 = y3 * powY37;
  double powY39 = y3 * powY38;
  double powY310 = y3 * powY39;



  double powAbsSquaredY2 = AbsSquaredY * AbsSquaredY;
  double powAbsSquaredY3 = AbsSquaredY * powAbsSquaredY2;
  double powAbsSquaredY4 = AbsSquaredY * powAbsSquaredY3;
  double powAbsSquaredY5 = AbsSquaredY * powAbsSquaredY4;


  //m_Ry = { 1, y3, -yb / 2., (-AbsSquaredY + 2 * powY32) / 4., -(y3 * yb) / 2., powYb2 / 8., -(AbsSquaredY * y3) / 4. + powY33 / 6., ((AbsSquaredY - 4 * powY32) * yb) / 16., (y3 * powYb2) / 8.,
  //    -powYb3 / 48., (3 * powAbsSquaredY2 + 8. * (powY34 - 3 * powY32 * powYb2)) / 192., -(y3 * (-3 * AbsSquaredY + 4 * powY32) * yb) / 48.,
  //    -((AbsSquaredY - 6 * powY32) * powYb2) / 96., -(y3 * powYb3) / 48., powYb4 / 384., (15 * powAbsSquaredY2 * y3 + 8. * (powY35 - 5 * powY33 * powYb2)) / 960.,
  //    -((powAbsSquaredY2 - 12 * AbsSquaredY * powY32 + 8 * powY34) * yb) / 384., (y3 * (-AbsSquaredY + 2 * powY32) * powYb2) / 96.,
  //    ((AbsSquaredY - 8 * powY32) * powYb3) / 768., (y3 * powYb4) / 384., -powYb5 / 3840.,
  //    (-5 * powAbsSquaredY3 + 90 * powAbsSquaredY2 * powY32 - 120 * AbsSquaredY * powY34 + 16 * powY36) / 11520.,
  //    -(y3 * (5 * powAbsSquaredY2 - 20 * AbsSquaredY * powY32 + 8 * powY34) * yb) / 1920., (powYb2 * (powAbsSquaredY2 + 16. * (powY34 - powY32 * powYb2))) / 3072.,
  //    -(y3 * (-3 * AbsSquaredY + 8 * powY32) * powYb3) / 2304., -((AbsSquaredY - 10 * powY32) * powYb4) / 7680., -(y3 * powYb5) / 3840., powYb6 / 46080. };

  //std::vector<std::complex<double>> out = { 1., y3, -0.5 * yb, 0.25 * (-1. * AbsSquaredY + 2. * powY32), -0.5 * y3 * yb, 0.125 * powYb2, 0.16666666666666666 * powY33 - 0.25 * AbsSquaredY * y3, 0.0625 * (AbsSquaredY - 4. * powY32) * yb,
  //   0.125 * powYb2 * y3, -0.020833333333333332 * powYb3, 0.005208333333333333 * (3. * powAbsSquaredY2 + 8. * (powY34 - 3. * powY32 * powYb2)),
  //   -0.020833333333333332 * (-3. * AbsSquaredY + 4. * powY32) * y3 * yb, -0.010416666666666666 * (AbsSquaredY - 6. * powY32) * powYb2, -0.020833333333333332 * powYb3 * y3, 0.0026041666666666665 * powYb4,
  //   0.0010416666666666667 * (8. * (powY35 - 5. * powY33 * powYb2) + 15. * powAbsSquaredY2 * y3), -0.0026041666666666665 * (powAbsSquaredY2 - 12. * AbsSquaredY * powY32 + 8. * powY34) * yb,
  //   0.010416666666666666 * (-1. * AbsSquaredY + 2. * powY32) * powYb2 * y3, 0.0013020833333333333 * (AbsSquaredY - 8. * powY32) * powYb3, 0.0026041666666666665 * powYb4 * y3,
  //   -0.00026041666666666666 * powYb5, 0.00008680555555555556 * (-5. * powAbsSquaredY3 + 90. * powAbsSquaredY2 * powY32 - 120. * AbsSquaredY * powY34 + 16. * powY36),
  //   -0.0005208333333333333 * (5. * powAbsSquaredY2 - 20. * AbsSquaredY * powY32 + 8. * powY34) * y3 * yb, 0.0003255208333333333 * powYb2 * (powAbsSquaredY2 + 16. * (powY34 - 1. * powY32 * powYb2)),
  //   -0.00043402777777777775 * (-3. * AbsSquaredY + 8. * powY32) * powYb3 * y3, -0.00013020833333333333 * (AbsSquaredY - 10. * powY32) * powYb4, -0.00026041666666666666 * powYb5 * y3,
  //   0.00002170138888888889 * powYb6 };
  std::vector<std::complex<double>> out;
  if (m_N <= 4)
  {
    out = {1., y3, -0.5 * yb, 0.25 * (-1. * AbsSquaredY + 2. * powY32), -0.5 * y3 * yb, 0.125 * powYb2, -0.25 * AbsSquaredY * y3 + 0.16666666666666666 * powY33,
      0.0625 * (AbsSquaredY - 4. * powY32) * yb, 0.125 * y3 * powYb2, -0.020833333333333332 * powYb3,
      0.005208333333333333 * (3. * powAbsSquaredY2 - 24. * AbsSquaredY * powY32 + 8. * powY34), -0.020833333333333332 * y3 * (-3. * AbsSquaredY + 4. * powY32) * yb,
      -0.010416666666666666 * (AbsSquaredY - 6. * powY32) * powYb2, -0.020833333333333332 * y3 * powYb3, 0.0026041666666666665 * powYb4  };

  }
  else
  {
    out = { 1., y3, -0.5 * yb, 0.25 * (-1. * AbsSquaredY + 2. * powY32), -0.5 * y3 * yb, 0.125 * powYb2, -0.25 * AbsSquaredY * y3 + 0.16666666666666666 * powY33,
  0.0625 * (AbsSquaredY - 4. * powY32) * yb, 0.125 * y3 * powYb2, -0.020833333333333332 * powYb3,
  0.005208333333333333 * (3. * powAbsSquaredY2 - 24. * AbsSquaredY * powY32 + 8. * powY34), -0.020833333333333332 * y3 * (-3. * AbsSquaredY + 4. * powY32) * yb,
  -0.010416666666666666 * (AbsSquaredY - 6. * powY32) * powYb2, -0.020833333333333332 * y3 * powYb3, 0.0026041666666666665 * powYb4,
  0.0010416666666666667 * y3 * (15. * powAbsSquaredY2 - 40. * AbsSquaredY * powY32 + 8. * powY34),
  -0.0026041666666666665 * (powAbsSquaredY2 - 12. * AbsSquaredY * powY32 + 8. * powY34) * yb, 0.010416666666666666 * y3 * (-1. * AbsSquaredY + 2. * powY32) * powYb2,
  0.0013020833333333333 * (AbsSquaredY - 8. * powY32) * powYb3, 0.0026041666666666665 * y3 * powYb4, -0.00026041666666666666 * powYb5,
  0.00008680555555555556 * (-5. * powAbsSquaredY3 + 90. * powAbsSquaredY2 * powY32 - 120. * AbsSquaredY * powY34 + 16. * powY36),
  -0.0005208333333333333 * y3 * (5. * powAbsSquaredY2 - 20. * AbsSquaredY * powY32 + 8. * powY34) * yb,
  0.0003255208333333333 * (powAbsSquaredY2 - 16. * AbsSquaredY * powY32 + 16. * powY34) * powYb2,
  -0.00043402777777777775 * y3 * (-3. * AbsSquaredY + 8. * powY32) * powYb3, -0.00013020833333333333 * (AbsSquaredY - 10. * powY32) * powYb4,
  -0.00026041666666666666 * y3 * powYb5, 0.00002170138888888889 * powYb6,
  0.00001240079365079365 * y3 * (-35. * powAbsSquaredY3 + 210. * powAbsSquaredY2 * powY32 - 168. * AbsSquaredY * powY34 + 16. * powY36),
  0.000010850694444444445 * (5. * powAbsSquaredY3 - 120. * powAbsSquaredY2 * powY32 + 240. * AbsSquaredY * powY34 - 64. * powY36) * yb,
  0.00002170138888888889 * y3 * (15. * powAbsSquaredY2 - 80. * AbsSquaredY * powY32 + 48. * powY34) * powYb2,
  -0.000010850694444444445 * (3. * powAbsSquaredY2 - 60. * AbsSquaredY * powY32 + 80. * powY34) * powYb3,
  0.00004340277777777778 * y3 * (-3. * AbsSquaredY + 10. * powY32) * powYb4, 0.000010850694444444445 * (AbsSquaredY - 12. * powY32) * powYb5,
  0.00002170138888888889 * y3 * powYb6, -1.5500992063492063e-6 * powYb7,
  1.937624007936508e-7 * (35. * powAbsSquaredY4 - 1120. * powAbsSquaredY3 * powY32 + 3360. * powAbsSquaredY2 * powY34 - 1792. * AbsSquaredY * powY36 +
    128. * powY38), -1.5500992063492063e-6 * y3 * (-35. * powAbsSquaredY3 + 280. * powAbsSquaredY2 * powY32 - 336. * AbsSquaredY * powY34 + 64. * powY36) * yb,
  -5.4253472222222224e-6 * (powAbsSquaredY3 - 30. * powAbsSquaredY2 * powY32 + 80. * AbsSquaredY * powY34 - 32. * powY36) * powYb2,
  -0.000010850694444444445 * y3 * (3. * powAbsSquaredY2 - 20. * AbsSquaredY * powY32 + 16. * powY34) * powYb3,
  2.7126736111111112e-6 * (powAbsSquaredY2 - 24. * AbsSquaredY * powY32 + 40. * powY34) * powYb4,
  -0.000010850694444444445 * y3 * (-1. * AbsSquaredY + 4. * powY32) * powYb5, -7.750496031746032e-7 * (AbsSquaredY - 14. * powY32) * powYb6,
  -1.5500992063492063e-6 * y3 * powYb7, 9.68812003968254e-8 * powYb8,
  2.152915564373898e-8 * y3 * (315. * powAbsSquaredY4 - 3360. * powAbsSquaredY3 * powY32 + 6048. * powAbsSquaredY2 * powY34 - 2304. * AbsSquaredY * powY36 +
    128. * powY38), -9.68812003968254e-8 * (7. * powAbsSquaredY4 - 280. * powAbsSquaredY3 * powY32 + 1120. * powAbsSquaredY2 * powY34 -
      896. * AbsSquaredY * powY36 + 128. * powY38) * yb, 7.750496031746032e-7 * y3 *
  (-7. * powAbsSquaredY3 + 70. * powAbsSquaredY2 * powY32 - 112. * AbsSquaredY * powY34 + 32. * powY36) * powYb2,
  4.5211226851851854e-7 * (powAbsSquaredY3 - 36. * powAbsSquaredY2 * powY32 + 120. * AbsSquaredY * powY34 - 64. * powY36) * powYb3,
  2.7126736111111112e-6 * y3 * (powAbsSquaredY2 - 8. * AbsSquaredY * powY32 + 8. * powY34) * powYb4,
  -1.937624007936508e-7 * (powAbsSquaredY2 - 28. * AbsSquaredY * powY32 + 56. * powY34) * powYb5,
  2.5834986772486774e-7 * y3 * (-3. * AbsSquaredY + 14. * powY32) * powYb6, 4.84406001984127e-8 * (AbsSquaredY - 16. * powY32) * powYb7, 9.68812003968254e-8 * y3 * powYb8,
  -5.382288910934745e-9 * powYb9, 1.0764577821869488e-9 * (-63. * powAbsSquaredY5 + 3150. * powAbsSquaredY4 * powY32 - 16800. * powAbsSquaredY3 * powY34 +
    20160. * powAbsSquaredY2 * powY36 - 5760. * AbsSquaredY * powY38 + 256. * powY310),
  -1.076457782186949e-8 * y3 * (63. * powAbsSquaredY4 - 840. * powAbsSquaredY3 * powY32 + 2016. * powAbsSquaredY2 * powY34 - 1152. * AbsSquaredY * powY36 +
    128. * powY38) * yb, 8.073433366402117e-9 * (7. * powAbsSquaredY4 - 336. * powAbsSquaredY3 * powY32 + 1680. * powAbsSquaredY2 * powY34 -
      1792. * AbsSquaredY * powY36 + 384. * powY38) * powYb2, -6.458746693121694e-8 * y3 *
  (-7. * powAbsSquaredY3 + 84. * powAbsSquaredY2 * powY32 - 168. * AbsSquaredY * powY34 + 64. * powY36) * powYb3,
  -3.229373346560847e-8 * (powAbsSquaredY3 - 42. * powAbsSquaredY2 * powY32 + 168. * AbsSquaredY * powY34 - 112. * powY36) * powYb4,
  -1.2917493386243387e-8 * y3 * (15. * powAbsSquaredY2 - 140. * AbsSquaredY * powY32 + 168. * powY34) * powYb5,
  4.0367166832010585e-9 * (3. * powAbsSquaredY2 - 96. * AbsSquaredY * powY32 + 224. * powY34) * powYb6,
  -1.6146866732804234e-8 * y3 * (-3. * AbsSquaredY + 16. * powY32) * powYb7, -2.6911444554673723e-9 * (AbsSquaredY - 18. * powY32) * powYb8,
  -5.382288910934745e-9 * y3 * powYb9, 2.691144455467372e-10 * powYb10 };

  }



    return out;
  

}

std::vector<std::complex<double>> RTable::evaluateRecursiveTableR(const Point& y)
{
    std::vector<std::complex<double>> out1(m_totNR);
    std::vector<std::complex<double>> out2((m_NR + 1) * (m_NR + 1));
    
    const double y3 = y.m_z;
    const std::complex<double> yb(y.m_x, y.m_y);
    const double r2 = y * y;
    const int N = m_NR ;
    out1[0] = 1;
    int idBef, idNext;
    double Const , invI;
    idBef = 0;
    // first diagonal 
    for (int i = 1; i < N + 1; i++)
    {
        idNext = idBef + i + 1; // (n+1)*n*0.5 + n
        invI = 1. / i;
        out1[idNext] = -1.*yb * 0.5 * invI * out1[idBef];
        idBef = idNext;
    }
    // second diagonal 
    // m = i -1
    idBef = 0;
    for (int i = 1; i < N + 1; i++)
    {
        idNext = idBef + i;
        Const = 1. / (2*i -1 ) ;
        out1[idNext] = Const * (2*i - 1) * y3 * out1[idBef];
        idBef = idBef + i + 1;
    }
    int constIdNM = 1;
    int constIdNMinusOneM = 0;
    int idNM , idNMinusOneM;
    for (int m = 0; m < N ; m++)
    {
        idNM = constIdNM;
        idNMinusOneM = constIdNMinusOneM;
        for (int i = m + 1; i < N; i++)
        {
            idNext = idNM + (i+1) ;
            Const = 1./((i + m + 1) * (i + 1 - m));
            out1[idNext] = Const*1.*(((2 * i) + 1) * y3 * out1[idNM] - r2 * out1[idNMinusOneM]);
            idNMinusOneM = idNM;
            idNM = idNext;
        }
        constIdNM = constIdNM + 3 + m;
        constIdNMinusOneM = constIdNMinusOneM + 2 + m;
    }
    // BUILD FULL TABLE
    int count, posMean, posMPos,posMNeg;
    count = 0;
    for (int n = 0; n < m_NR + 1; n++)
    {
        posMean = n * n + n; // this is the position of middle elements on pascal triangle
        out2[posMean] = out1[count];
        count++;
        for (int m = 1; m < n + 1; m++)
        {
            posMPos = posMean + m;  // positive and negative indexs on pascal triangle
            posMNeg = posMean - m;
            std::complex<double>* Ry = &out1[count];
            out2[posMPos] = *Ry;
            short int minus1 = (m % 2) == 0 ? 1 : -1;
            out2[posMNeg] = (minus1 * 1. * std::conj(*Ry));
            count++;
        }
    }
    return out2;
}

std::vector<std::complex<double>>  RTable::evaluateRecursiveTableS(const Point& x)
{
  std::vector<std::complex<double>> out((m_NS + 1) * (m_NS + 1));
  std::vector<std::complex<double>> out1 ( m_totNS);
   

    const double x3 = x.m_z;
    const std::complex<double> xb(x.m_x, -x.m_y);
    const double r2 = x * x;
    const double ir2 = 1 / (x * x);
    const int N = m_NS;
    out1[0] = 1/ (x.getNorm());
    //
    int idBef, idNext , Const;
    // first diagonal
    idBef = 0;
    for (int i = 1; i < N + 1; i++)
    {
        idNext = idBef + i + 1; // (n+1)*n*0.5 + n
        out1[idNext] = (2 * i - 1) * 1. * xb * ir2 * -out1[idBef];
        idBef = idNext;
    }
    // second diagonal 
    // m = i -1
    idBef = 0;
    for (int i = 1; i < N + 1; i++)
    {
        idNext = idBef + i;
        out1[idNext] = ir2 * (2 * i - 1) * x3 * out1[idBef];
        idBef = idBef + i + 1;
    }
    int constIdNM = 1;
    int constIdNMinusOneM = 0;
    int idNM, idNMinusOneM;
    for (int m = 0; m < N; m++)
    {
        idNM = constIdNM;
        idNMinusOneM = constIdNMinusOneM;
        for (int i = m + 1; i < N; i++)
        {
            idNext = idNM + (i + 1);
            int Const = (i+m)*(i-m);
            out1[idNext] = ir2 * ( ( 2 * i + 1)  * x3 * out1[idNM] - Const * 1. * out1[idNMinusOneM]);
            idNMinusOneM = idNM;
            idNM = idNext;
        }
        constIdNM = constIdNM + 3 + m;
        constIdNMinusOneM = constIdNMinusOneM + 2 + m;
    }
    // BUILD FULL TABLE
    int count, posMean, posMPos, posMNeg;
    count = 0;
    for (int n = 0; n < m_NS + 1; n++)
    {
        posMean = n * n + n; // this is the position of middle elements on pascal triangle
        out[posMean] = out1[count];
        count++;
        for (int m = 1; m < n + 1; m++)
        {
            posMPos = posMean + m;  // positive and negative indexs on pascal triangle
            posMNeg = posMean - m;
            std::complex<double>* Sx = &out1[count];
            out[posMPos] = *Sx;
            short int minus1 = (m % 2) == 0 ? 1 : -1;
            out[posMNeg] = (minus1 * 1. * std::conj(*Sx));
            count++;
        }
    }
    return out;
}

std::complex<double> RTable::getPosCdR(const int n, const int m, std::vector<std::complex<double>> & _R)
{
/*    static std::vector<std::complex<double>>* R;
    if (!(_R == nullptr))
    {
      R = _R;
      return 0.0;
    }   */ 

    if (m > n - 2 )
      return 0.0;
    
    const int new_n = n - 1;
    const int new_m = m + 1;
    int pos = getIdNM(new_n, abs(new_m));
    if (new_m >= 0)
    {
      return _R[pos];
    }
    else
    {
      short int minus1 = (new_m % 2) == 0 ? 1 : -1;
      return (minus1 * 1. * std::conj(_R[pos]));
    }
    
}

std::complex<double> RTable::getProjdR(const int n, const int m, std::vector<std::complex<double>>& R)
{
  if (std::abs(m-1) > (n-1) )
    return 0.0;

  const int new_n = n - 1;
  const int new_m = m - 1 ;
  int pos = getIdNM(new_n, abs(new_m));
  if (new_m >= 0 )
  {
    return R[pos];
  }
  else
  {
    short int minus1 = (new_m % 2) == 0 ? 1 : -1;
    return (minus1 * 1. * std::conj(R[pos]));
  }
}

int RTable::getPosZdR(const int n, const int m)
{
  if (m < -n + 1 || m > n - 1)
  {
    return -1;
  }
  
  return getIdNM(n - 1, m);


}

const int RTable::getIdNM(const int& N, const int& M)
{
    //const int MTotal = 2 * N + 1;
    //const int MVirtual = M + (MTotal - 1) / 2;

    if (M < 0)
    {
        return -1;
    }

    // convert the n,m index to a position at a row vector;
    const int totN = (N+1)*(N)*0.5;

    return (totN  + M);
}
