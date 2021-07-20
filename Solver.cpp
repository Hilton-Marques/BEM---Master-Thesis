#include <cmath>
#include "Solver.h"
#include "GaussQuadrature.h"
#include "RTable.h"
#include <omp.h>
#include <iostream>

extern const double M_PI;
extern const double Const;

Solver::Solver(const int &N, const int &nGFMM, const int &nNodes): m_N(N), m_nGFMM(nGFMM), m_nNodes(nNodes)
{
    m_NTot = (N + 1) * (N + 1);
    m_RTable.setSizeRTable(N);
    m_RTable.setSizeRecTableR(N);
    m_RTable.setSizeRecTableS(N);
    m_gauss.init(m_nGBEM);
}

Solver::~Solver()
{
    delete &m_Hd;
    delete &m_Gt;
}

void Solver::CalculateME(Face* element)
{
    //element->m_fatherElement->setFMMSizes(m_NTot);
    Face* fatherElement = element->m_fatherElement;
    element->setFMMSizes(m_NTot);
    GaussQuad gauss(m_nGFMM);
    Point yc = element->m_fatherElement->getYc();

    Point* normal = &element->m_normal;
    const std::complex<double> CNormal(normal->m_x, -normal->m_y);
    const double ZNormal = normal->m_z;
    double* jacobian = & element->m_jacobian;


    for (int i = 0; i < 3; i++)
    {
        Vertex* field = element->m_points[i];
        if (field->m_bd)
        {
            const double u = field->m_u; // temperature

            for (int j = 0; j < gauss.m_totN; j++)
            {
                double* Xi = &gauss.m_gaussPointsXi[j];
                double* Eta = &gauss.m_gaussPointsEta[j];
                double* W = &gauss.m_gaussWeights[j];
                double shapeF[3] = { 1. - *Xi - *Eta , *Xi, *Eta };
                Point xksi = Multi(shapeF, element->m_points);
                Point r = xksi - yc;
                //RTable Ry(m_N, r);
                double Const2 = u * shapeF[i] * (*W);
                int zdRpos;
                std::complex<double> adR, bdR, zdR;
                std::complex<double> xdR, ydR;
                std::complex<double> I(0, 1);
                int count = 0;  // carry the index of array Ry
                std::vector<std::complex<double>> R = m_RTable.evaluateTableR(r);
               
               // m_RTable.getPosCdR(-1, -1, R);
                for (int n = 0; n < m_N + 1; n++)
                {

                  int posMean = n * n + n; // this is the position of middle elements on pascal triangle

                  zdRpos = m_RTable.getPosZdR(n, 0);
                  zdR = zdRpos != -1 ? R[zdRpos] : 0.0;

                  bdR = m_RTable.getPosCdR(n, 0, R);
                  adR = m_RTable.getProjdR(n, 0, R);

                  xdR = -0.5 * (adR - bdR);
                  ydR = -0.5 * I*(adR + bdR);

                  //element->m_MEH[posMean] += +Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);
                  fatherElement->m_MEH[posMean] +=  + Const2 * (xdR*normal->m_x + ydR*normal->m_y + zdR*normal->m_z);
                   

                  for (int m = 1; m < n + 1; m++)
                  {                                       
                    int posMPos = posMean + m;  // positive and negative indexs on pascal triangle
                    int posMNeg = posMean - m;                  
                    
                    zdRpos = m_RTable.getPosZdR(n, m);

                    short int minus1 = (m % 2) == 0 ? 1 : -1;

                    zdR = zdRpos != -1 ? minus1 *1.0* std::conj(R[zdRpos]) : 0.0;

                    
                    bdR = m_RTable.getPosCdR(n, -m , R);
                    
                    adR = m_RTable.getProjdR(n, -m, R);


                    xdR = -0.5 * (adR - bdR);
                    ydR = -0.5 * I*(adR + bdR);

                    //element->m_MEH[posMNeg] += Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);
                    fatherElement->m_MEH[posMNeg] += Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);


                    zdR = zdRpos != -1 ? R[zdRpos] : 0.0;

                    bdR = m_RTable.getPosCdR(n, m , R);
                    adR = m_RTable.getProjdR(n, m, R);


                    xdR = -0.5 * (adR - bdR);
                    ydR = -0.5 * I *(adR + bdR);

                    //element->m_MEH[posMPos] += Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);
                    fatherElement->m_MEH[posMPos] +=  Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);



                    count++;
                  }
                }

        
            } 


        }
        if (element->m_bd)
        {       
            const double q = element->m_q[i]; // gradient
            const double Const = (*jacobian) * q;

                for (int j = 0; j < gauss.m_totN; j++)
                {
                    double* Xi = &gauss.m_gaussPointsXi[j];
                    double* Eta = &gauss.m_gaussPointsEta[j];
                    double* W = &gauss.m_gaussWeights[j];
                    double shapeF[3] = {  1. - *Xi - *Eta, *Xi, *Eta };
                    Point xksi = Multi(shapeF, element->m_points);
                    const Point r = xksi - yc;
                    double Const2 = *W * Const * shapeF[i];
                    int count = 0 ;  // carry the index of array Ry
                    std::vector<std::complex<double>> R = m_RTable.evaluateTableR(r);
                    for (int n = 0; n < m_N + 1; n++) 
                    {
                        int posMean = n * n + n; // this is the position of middle elements on pascal triangle

                       // element->m_MEG[posMean] += Const2 * (R[count]);
                        
                        fatherElement->m_MEG[posMean] +=  Const2 *( R [count]);
                        
                        count++;
                        for (int m = 1; m < n + 1; m++)
                        {
                           int posMPos = posMean + m;  // positive and negative indexs on pascal triangle
                           int posMNeg = posMean - m;
                           //element->m_MEG[posMPos] += Const2 * (R[count]);
                           fatherElement->m_MEG[posMPos] +=  Const2 * (R[count]);

                           short int minus1 = (m % 2) == 0 ? 1 : -1; 
                           //element->m_MEG[posMNeg] += Const2 * (minus1 * 1. * std::conj((R[count])));
                           fatherElement->m_MEG[posMNeg] +=  Const2 * (minus1*1.* std::conj((R[count])) );

                           count++;
                        }
                    }

                }
            

            
        }
    }



}

void Solver::CalculateNearInt(std::vector<Vertex*> sourcesNodes, Face* element)
{
    
    Point* normal = &element->m_normal;
    double* jacobian = &element->m_jacobian;
    const double ConstG = Const  * (*jacobian);

    
    for (int i = 0; i < 3; i++)
    {
      Vertex* pt = element->m_points[i];

      if (pt->m_bd)
      {
        const double u = pt->m_u;
        const double Const2 = Const * u;

        for (Vertex* source : sourcesNodes)
        {
          // Calculate diagonal of H Matrix

          if (source->m_id == pt->m_id)
          {
            if (!source->markTemp)
            
            {
              source->markTemp = true;
              #pragma omp atomic
              m_Hd[{source->m_id, 0}] +=  Const2* pt->calculateSolidAngle();
             
            }
          }
          else
          {
            // Integrate Hu
            double H = 0.0;
            for (int j = 0; j < m_gauss.m_totN; j++)
            {
              double* Xi = &m_gauss.m_gaussPointsXi[j];
              double* Eta = &m_gauss.m_gaussPointsEta[j];
              double* W = &m_gauss.m_gaussWeights[j];
              double shapeF[3] = { 1. - *Xi - *Eta, *Xi, *Eta, };
              Point xksi = Multi(shapeF, element->m_points);
              Point r = xksi - source->m_coord;
              double rNorm = r.getNorm();
              double rNorm3 = rNorm * rNorm * rNorm;
              H = H + (- (r * (*normal) ) / rNorm3 ) * (*W) * shapeF[i];
            }

            #pragma omp atomic
            m_Hd[{source->m_id, 0}] +=  Const2 * H;

          }
        }

      }
      if (element->m_bd)
      {
        const double q = element->m_q[i];
        for (Vertex* source : sourcesNodes)
        {          

          // Integrate Gq
          double G = 0.0;
          for (int j = 0; j < m_gauss.m_totN; j++)
          {
            double* Xi = &m_gauss.m_gaussPointsXi[j];
            double* Eta = &m_gauss.m_gaussPointsEta[j];
            double* W = &m_gauss.m_gaussWeights[j];
            double shapeF[3] = { 1. - *Xi - *Eta , *Xi, *Eta };
            Point xksi = Multi(shapeF, element->m_points);
            Point r = xksi - source->m_coord;
            double rNorm = r.getNorm();
            G = G + (1 / rNorm) * (*W) * shapeF[i];
          }
          #pragma omp atomic
          m_Gt[{source->m_id, 0}] += q *ConstG * G;
        }

      }
    }


}

void Solver::CalculateFarInt(std::vector<Vertex*> sourcesNodes, Face* element)
{
    //std::vector<Face*> grangranChildrenElement;
    //element->getLeafElements(&grangranChildrenElement);
    Point yc = element->getYc();
    std::vector<std::complex<double>> Sb;
    for (Vertex* source : sourcesNodes) 
    {
      //for (Face* gran : grangranChildrenElement)
      //{
      //  source->m_closeElements.push_back(gran);
      //}
      Point x = source->m_coord - yc;
      Sb = m_RTable.evaluateRecursiveTableS(x);
      //double valueG = Const * Dot(element->m_MEG, &Sb);
      //double valueH = Const * Dot(element->m_MEH, &Sb);
#pragma omp atomic
      m_Gt[{source->m_id, 0}] += Const * Dot(element->m_MEG, &Sb);
#pragma omp atomic
      m_Hd[{source->m_id, 0}] += Const * Dot(element->m_MEH, &Sb);




      
    }
}

void Solver::CalculateMMT(Face* element)
{    
    Point ycL = element->getYc();
    for (Face* child : element->m_childrenElement)
    {
        Point yc = child->getYc();
        Point r = yc - ycL;
        int id;
        std::vector<std::complex<double>> R = m_RTable.evaluateRecursiveTableR(r);
        int count = 0;
        for (int n = 0; n < m_N + 1; n++)
        {
            for (int m = -n; m < n+1; m++)
            {
                for (int nc = 0; nc < m_N + 1; nc++)
                {
                    for (int mc = -nc; mc < nc + 1; mc++)
                    {
                        id = getPos(n - nc, m - mc);
                        if ( id > -1)
                        {
                            element->m_MEG[count] = element->m_MEG[count] + R[getPos(nc, mc)] * child->m_MEG[id];
                            element->m_MEH[count] = element->m_MEH[count] + R[getPos(nc, mc)] * child->m_MEH[id];
                            
                        }
                    }
                }
                count++;
            }
        }
        child->reset();
    }
}



double Solver::Dot(std::complex<double> R[], std::vector<std::complex<double>>* Sb)
{
  std::complex <double> dot = 0.0;
  for (int i = 0; i < m_NTot; i++)
  {
    dot = dot + (R[i] * (*Sb)[i]);
  }
  return std::real(dot);
}

int Solver::getPos(const int& N, const int& M)
{

        const int mTot = 2 * N + 1;
        const int mVirtual = M + (mTot - 1) * 0.5 ;
        if (mVirtual >= mTot || mVirtual < 0)
        {
            return -1;
        }
        const int nTot = (N ) * (N );
        return (nTot + mVirtual);

 
}




std::complex<double> operator*(const std::vector<std::complex<double>>& dR,const Point& x)
{
    return std::complex<double>(dR[0]*x.m_x + dR[1]*x.m_y + dR[2]*x.m_z);
}

void Solver::prepareElementMatrix(std::vector<Vertex*> sourceNodes , Face* element)
{
  Point* normal = &element->m_normal;
  double* jacobian = &element->m_jacobian;
  Matrix GG(sourceNodes.size(), 3);
  Matrix HH(sourceNodes.size(), 3);

  for (int j = 0; j < 3; j++)
  {
    Vertex* pt = element->m_points[j];

    if (pt->m_bd)
    {
      const double Const2 = Const;

     for (size_t i = 0; i < sourceNodes.size(); i++)
       {
       Vertex* source = sourceNodes[i];
        if (source->m_id == pt->m_id)
        {
          if (!source->markTemp)
          {
            HH[{i,j}] = Const2 * pt->calculateSolidAngle();
            source->markTemp = true;
          }
          // Calculate diagonal of H Matrix
        }
        else
        {
          // Integrate Hu
          double H = 0.0;
          for (int k = 0; k < m_gauss.m_totN; k++)
          {
            double* Xi = &m_gauss.m_gaussPointsXi[k];
            double* Eta = &m_gauss.m_gaussPointsEta[k];
            double* W = &m_gauss.m_gaussWeights[k];
            double shapeF[3] = { 1. - *Xi - *Eta, *Xi, *Eta, };
            Point xksi = Multi(shapeF, element->m_points);
            Point r = xksi - source->m_coord;
            double rNorm = r.getNorm();
            double rNorm3 = rNorm * rNorm * rNorm;
            H = H + (-(r * (*normal)) / rNorm3) * (*W) * shapeF[j];
          }
          HH[{i, j}] = Const2 * H;
        }
      }
    }
    if (element->m_bd)
    {
      const double Const2 = Const * (*jacobian);
      for (size_t i = 0; i < sourceNodes.size(); i++)
      {
        Vertex* source = sourceNodes[i];
        // Integrate Gq
        double G = 0.0;
        for (int k = 0; k < m_gauss.m_totN; k++)
        {
          double* Xi = &m_gauss.m_gaussPointsXi[k];
          double* Eta = &m_gauss.m_gaussPointsEta[k];
          double* W = &m_gauss.m_gaussWeights[k];
          double shapeF[3] = { 1. - *Xi - *Eta , *Xi, *Eta };
          Point xksi = Multi(shapeF, element->m_points);
          Point r = xksi - source->m_coord;
          double rNorm = r.getNorm();
          G = G + (1 / rNorm) * (*W) * shapeF[j];
        }
        GG[{i, j}] =  Const2 * G;
      }

    }
  }
  element->m_H = HH;
  element->m_G = GG;

}
void Solver::CalculateNearIntMatrix(std::vector<Vertex*> sourceNodes , Face* element)
{
  Matrix d(3, 1);
  Matrix t(3, 1);
  for (int i = 0; i < 3; i++)
  {
    Vertex* pt = element->m_points[i];
    if (pt->m_bd)
    {
      d[{i, 0}] = pt->m_u;
    }
    if (element->m_bd)
    {
      t[{i, 0}] = element->m_q[i];
    }
  }
  for (int i = 0; i < sourceNodes.size(); i++)
  {
    Vertex* pt = sourceNodes[i];
    Matrix rowH = element->m_H.getSubMatrix(i + 1, i + 1, 1, 3);
    Matrix rowG = element->m_G.getSubMatrix(i+1, i+1, 1, 3);
    #pragma omp atomic
    m_Hd[{pt->m_id, 0}] +=  rowH.dot(d);
    #pragma omp atomic
    m_Gt[{pt->m_id, 0}] += rowG.dot(t);
  }

}
void Solver::CalculateMEForGMRES(Face* element)
{
  GaussQuad gauss(m_nGFMM);
  Point yc = element->m_fatherElement->getYc();
  Point* normal = &element->m_normal;
  const std::complex<double> CNormal(normal->m_x, -normal->m_y);
  const double ZNormal = normal->m_z;
  double* jacobian = &element->m_jacobian;


  for (int i = 0; i < 3; i++)
  {
    Vertex* field = element->m_points[i];
    std::complex<double>* MEH_i = element->m_MEH_p[i];
    std::complex<double>* MEG_i = element->m_MEG_p[i];

    if (field->m_bd)
    {
      const double u = field->m_u; // temperature

      for (int j = 0; j < gauss.m_totN; j++)
      {
        double* Xi = &gauss.m_gaussPointsXi[j];
        double* Eta = &gauss.m_gaussPointsEta[j];
        double* W = &gauss.m_gaussWeights[j];
        double shapeF[3] = { 1. - *Xi - *Eta , *Xi, *Eta };
        Point xksi = Multi(shapeF, element->m_points);
        Point r = xksi - yc;
        //RTable Ry(m_N, r);
        double Const2 = u * shapeF[i] * (*W);
        int zdRpos;
        std::complex<double> adR, bdR, zdR;
        std::complex<double> xdR, ydR;
        std::complex<double> I(0, 1);
        int count = 0;  // carry the index of array Ry
        std::vector<std::complex<double>> R = m_RTable.evaluateTableR(r);

        // m_RTable.getPosCdR(-1, -1, R);
        for (int n = 0; n < m_N + 1; n++)
        {

          int posMean = n * n + n; // this is the position of middle elements on pascal triangle

          zdRpos = m_RTable.getPosZdR(n, 0);
          zdR = zdRpos != -1 ? R[zdRpos] : 0.0;

          bdR = m_RTable.getPosCdR(n, 0, R);
          adR = m_RTable.getProjdR(n, 0, R);

          xdR = -0.5 * (adR - bdR);
          ydR = -0.5 * I * (adR + bdR);

          MEH_i[posMean] += +Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);



          for (int m = 1; m < n + 1; m++)
          {
            int posMPos = posMean + m;  // positive and negative indexs on pascal triangle
            int posMNeg = posMean - m;

            zdRpos = m_RTable.getPosZdR(n, m);

            short int minus1 = (m % 2) == 0 ? 1 : -1;

            zdR = zdRpos != -1 ? minus1 * 1.0 * std::conj(R[zdRpos]) : 0.0;


            bdR = m_RTable.getPosCdR(n, -m, R);

            adR = m_RTable.getProjdR(n, -m, R);


            xdR = -0.5 * (adR - bdR);
            ydR = -0.5 * I * (adR + bdR);

            MEH_i[posMNeg] += Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);


            zdR = zdRpos != -1 ? R[zdRpos] : 0.0;

            bdR = m_RTable.getPosCdR(n, m, R);
            adR = m_RTable.getProjdR(n, m, R);


            xdR = -0.5 * (adR - bdR);
            ydR = -0.5 * I * (adR + bdR);

            MEH_i[posMPos] += Const2 * (xdR * normal->m_x + ydR * normal->m_y + zdR * normal->m_z);



            count++;
          }
        }


      }


    }
    if (element->m_bd)
    {
      const double q = element->m_q[i]; // gradient
      const double Const = (*jacobian) * q;

      for (int j = 0; j < gauss.m_totN; j++)
      {
        double* Xi = &gauss.m_gaussPointsXi[j];
        double* Eta = &gauss.m_gaussPointsEta[j];
        double* W = &gauss.m_gaussWeights[j];
        double shapeF[3] = { 1. - *Xi - *Eta, *Xi, *Eta };
        Point xksi = Multi(shapeF, element->m_points);
        const Point r = xksi - yc;
        double Const2 = *W * Const * shapeF[i];
        int count = 0;  // carry the index of array Ry
        std::vector<std::complex<double>> R = m_RTable.evaluateTableR(r);
        for (int n = 0; n < m_N + 1; n++)
        {
          int posMean = n * n + n; // this is the position of middle elements on pascal triangle

          MEG_i[posMean] = element->m_MEG[posMean] + Const2 * (R[count]);
          count++;
          for (int m = 1; m < n + 1; m++)
          {
            int posMPos = posMean + m;  // positive and negative indexs on pascal triangle
            int posMNeg = posMean - m;
            MEG_i[posMPos] = element->m_MEG[posMPos] + Const2 * (R[count]);

            short int minus1 = (m % 2) == 0 ? 1 : -1;
            MEG_i[posMNeg] = element->m_MEG[posMNeg] + Const2 * (minus1 * 1. * std::conj((R[count])));
            count++;
          }
        }

      }



    }
  }



}

void Solver::CalculateNearIntExact(std::vector<Vertex*> sourcesNodes, Face* element)
{
  double jacobian = element->m_jacobian;
  double jac2 = element->m_jacobian * element->m_jacobian;
  double tol = 1.0e-10;
  Point u = (element->m_points[0]->m_coord - element->m_points[2]->m_coord);
  Point v = (element->m_points[1]->m_coord - element->m_points[2]->m_coord);

  for (Vertex* source : sourcesNodes)
  {
    
    // Calculate Integral
    Point x0 = source->m_coord - element->m_points[2]->m_coord;
    double X00 = x0 * x0;
    Matrix A(3, 2);
    Matrix M(3, 3);
    for (int i = 0; i < 3; i++)
    {
      A[{ i, 0 }] = u[i]; A[{i, 1}] = v[i];
      M[{i, 0}] = u[i]; M[{i, 1}] = v[i];  M[{i, 2}] = x0[i];
    }
    Matrix At = A.transpose();
    Matrix X = At * A;
    double X01 = u * x0; double X02 = v * x0;
    double ab1 = (X[{1, 1}] * X01 - X[{0, 1}] * X02) / jac2;
    double ab2 = (X[{0, 0}] * X02 - X[{0, 1}] * X01) / jac2;
    double d11 = M[{1, 1}] * M[{2, 2}] - M[{1, 2}] * M[{2, 1}];
    double d22 = M[{1, 0}] * M[{2, 2}] - M[{1, 2}] * M[{2, 0}];
    double d33 = M[{1, 0}] * M[{2, 1}] - M[{1, 1}] * M[{2, 0}];
    double d00 = M[{0, 0}] * d11 - M[{0, 1}] * d22 + M[{0, 2}] * d33;
    double d = X00 - X01 * ab1 - X02 * ab2;
    double Atil[3] = { X[{1,1}], X[{0,0}], X[{0,0}] - 2 * X[{0,1}] + X[{1,1}] };
    double Btil[3] = { 2 * (X02 - X[{1,1}]), -2 * X01, 2 * (X[{0,1}] - X[{0,0}] - X02 + X01) };
    double Ctil[3] = { X[{1,1}] - 2 * X02 + X00, X00, X[{0,0}] - 2 * X01 + X00 };
    Matrix Ns1(3, 3, { ab1,-ab1,0, ab2,1 - ab2,-1, 1 - ab1 - ab2,ab1 + ab2 - 1,1 });
    Matrix Ns2(3, 3, { ab1,-ab1,1, ab2,-ab2,0, 1 - ab1 - ab2,ab1 + ab2,-1 });
    Matrix Ns3(3, 3, { ab1,1 - ab1,-1, ab2,-ab2,1, 1 - ab1 - ab2,ab1 + ab2 - 1,0 });
    std::vector<Matrix> Ns = { Ns1, Ns2, Ns3 };
    double jacL[3] = { ab1, ab2, 1 - ab1 - ab2 };

    Matrix outH(3, 1);  // H integral
    Matrix outG(3, 1);  // G integral

    if (std::abs(d00) > tol)
    {
      for (int i = 0; i < 3; i++)
      {
        double jacI = jacL[i];
        if (std::abs(jacI) > tol)
        {
          Matrix IH(3, 1);
          Matrix IG(3, 1);

          double a = Atil[i];
          double b = Btil[i];
          double c = Ctil[i];
          std::complex<double> sqrtA = sqrt(a);
          std::complex<double> sqrtD = sqrt(d);
          std::complex<double> sqrtC = sqrt(c);
          std::complex<double> sqrtAC = sqrt(a * c);
          std::complex<double> sqrtAD = sqrtA * sqrtD;
          std::complex<double> sqrtCD = sqrt(c * d);
          std::complex<double> sqrtABC = sqrt(a + b + c);
          std::complex<double> sqrtAABC = sqrtA * sqrtABC;
          std::complex<double> sqrtDABC = sqrtD * sqrtABC;
          std::complex<double> b24ac = b * b - 4. * a * c;
          std::complex<double> sqrtb24ac4ad = sqrt((b24ac + 4. * a * d));
          std::complex<double> sqrtCMD = sqrt(c - d);
          std::complex<double> sqrtABCMD = sqrt(a + b + c - d);
          std::complex<double> logB2SqrtAC = log(b + 2. * sqrtAC);
          std::complex<double> log2ABSQRTABC = log(2. * a + b + 2. * sqrtAABC);
          std::complex<double> logABCABCMD = log(sqrtABC + sqrtABCMD);
          std::complex<double> logSqrtCSqrtCMD = log(sqrtC + sqrtCMD);
          std::complex<double> logD = log(d);
          std::complex<double> B24ADMC = b24ac + 4. * a * d;
          std::complex<double> sqrtB24ADMC = sqrt(B24ADMC);


          IH[{0, 0}] = std::real((-log(4. * a * (c + sqrtCD) + b * (-b + sqrtb24ac4ad))
            + log(4. * a * (c + sqrtCD) - b * (b + sqrtb24ac4ad))
            + -log(-b * (b + sqrtb24ac4ad) - 2. * a * (-2. * c - 2. * sqrtABC *
              sqrtD + sqrtb24ac4ad)) + log(b * (-1. * b + sqrtb24ac4ad) + 2. * a * (2. * c + 2. * sqrtDABC
                + sqrtb24ac4ad))) / (sqrtD * sqrtB24ADMC));
          IH[{1, 0}] = std::real((4. * sqrtA * (logB2SqrtAC - log2ABSQRTABC) -
            2. * b * logSqrtCSqrtCMD / sqrtCMD +
            2. * (2. * a + b) * logABCABCMD / sqrtABCMD +
            logD * (b * (1. / sqrtCMD - 1. / sqrtABCMD) - (2. * a / sqrtABCMD))) /
            -B24ADMC);
          IH[{2, 0}] = std::real(((2. * b * (log2ABSQRTABC - logB2SqrtAC) / sqrtA) +
            4. * sqrtCMD * logSqrtCSqrtCMD - 2. * sqrtCMD * logD -
            ((b + 2. * c - 2. * d) * (2. * logABCABCMD - logD) / sqrtABCMD)) / -B24ADMC);

          outH += Ns[i] * IH * -jacI;

          IG[{0, 0}] = std::real((-sqrtb24ac4ad * logB2SqrtAC + (sqrtb24ac4ad)*log2ABSQRTABC
            - sqrtAD * log(-b24ac + 4. * a * sqrtCD - b * (sqrtb24ac4ad)) +
            sqrtAD * log(-b24ac + 4. * a * sqrtDABC - 2. * a * (sqrtb24ac4ad)-b * (sqrtb24ac4ad)) +
            sqrtAD * log(-b24ac + 4. * a * sqrtCD + b * (sqrtb24ac4ad)) -
            sqrtAD * log(-b24ac + 4. * a * sqrtDABC + 2. * a * (sqrtb24ac4ad)+b * (sqrtb24ac4ad))) / (sqrtA * sqrtb24ac4ad));

          IG[{1, 0}] = std::real((b24ac * sqrtCMD * sqrtABCMD *
            logB2SqrtAC - b24ac * sqrtCMD * sqrtABCMD * log2ABSQRTABC +
            sqrtA * d * (2. * b * sqrtABCMD * logSqrtCSqrtCMD -
              2. * (2. * a + b) * sqrtCMD * logABCABCMD +
              (2. * a * sqrtCMD + b * sqrtCMD - b * sqrtABCMD) * logD))
            / (2. * sqrtA * (-b24ac - 4. * a * d) * sqrtCMD * sqrtABCMD));

          IG[{2, 0}] = std::real((-b * b24ac * sqrtCMD * sqrtABCMD *
            logB2SqrtAC - b * b24ac * sqrtCMD * sqrtABCMD * log2ABSQRTABC +
            2. * sqrtA * ((sqrtC - sqrtABC) * sqrtABCMD
              * sqrtb24ac4ad * sqrtb24ac4ad + a * d * (-4. * sqrtCMD * sqrtABCMD * logSqrtCSqrtCMD + (b + 2. * c + -2. * d)
                * (2. * logABCABCMD - logD) + 2. * sqrtCMD * sqrtABCMD * logD))) /
            (4. * sqrtA * a * (-b24ac - 4. * a * d) * sqrtABCMD));

          outG += Ns[i] * IG * jacI;

        }

      }
    }
    else
    {
      for (int i = 0; i < 3; i++)
      {
        double jacI = jacL[i];
        if (std::abs(jacI) > tol)
        {
          Matrix I(3, 1);
          double a = Atil[i];
          double b = Btil[i];
          double c = Ctil[i];
          std::complex<double> sqrtA = sqrt(a);
          std::complex<double> sqrtAC = sqrt(a * c);
          std::complex<double> sqrtABC = sqrt(a + b + c);
          std::complex<double> sqrtAABC = sqrtA * sqrtABC;
          std::complex<double> logB2SqrtAC = log(b + 2. * sqrtAC);
          std::complex<double> log2ABSQRTABC = log(2. * a + b + 2. * sqrtAABC);

          I[{0, 0}] = std::real((log2ABSQRTABC - logB2SqrtAC) / sqrtA);
          I[{1, 0}] = I[{0, 0}] * 0.5;
          I[{2, 0 }] = std::real((-2. * sqrtAC + 2. * sqrtAABC - b * sqrtA * I[{0, 0}]) / (4. * sqrtA * a));

          outG += Ns[i] * I * jacI;
        }

      }
    }
   outH = outH * d00;
   const double ConstG = Const  * (jacobian);


      // MatrixVector multiplication
      for (int i = 0; i < 3; i++)
      {
        Vertex* pt = element->m_points[i];
        if (pt->m_bd)
        {
          const double u = pt->m_u;
          const double Const2 = Const * u;

          if (source->m_id == pt->m_id)
          {
            // Calculate diagonal of H Matrix
            if (!source->markTemp)
            {
              source->markTemp = true;
              #pragma omp atomic
              m_Hd[{source->m_id, 0}] += Const2 * pt->calculateSolidAngle();

            }
          }
          else
          {
            #pragma omp atomic
            m_Hd[{source->m_id, 0}] += Const2 * outH[{i, 0}];

          }
        }
        if (element->m_bd)
        {
          const double q = element->m_q[i];
          #pragma omp atomic
          m_Gt[{source->m_id, 0}] += q * ConstG * outG[{i, 0}];
        }
      }
    }
  

}
