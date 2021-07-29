//
// Created by hilto on 27/01/2021.
//

#include <algorithm> // find function
#include <iostream>
#include "Face.h"
#include <set> // erase duplicates


Face::Face() {

}

Face::Face(int cHedInc, double q[] , short int level) {
    m_hedInc = cHedInc;
    m_nL = level;
    m_heds.resize(3); // cant call a method inside class declaration
    if (!(q[0] == NAN))
    {
      m_q[0] = q[0];
      m_q[1] = q[1];
      m_q[2] = q[2];
      m_bd = true;
    }

    //getIdcPoints(heds);

}


std::vector<Face*> Face::getAdjacentElements(bool flag)
{
    //std::vector<int> elIds;
    std::vector<Face*> adjacentElements;
    //adjacentElements.reserve(15);
    for (int i = 0; i < 3; i++)
    {
        HalfEdge* hed = m_heds[i];
        m_adjacentElementsEdges[i] = hed->m_edge->getTwin(hed->m_id)->m_el;
        int p2Target = hed->m_inc[1];
        // get next of twin
        hed = hed->m_edge->getTwin(hed->m_id)->m_heNext;
        while (true)
        {
            int p2 = hed->m_inc[1];
            if (p2 == p2Target)
            {
                break;
            }
            //int newId = hed->m_elId;
            Face* element = hed->m_el;
            adjacentElements.push_back(element);
      /*      if (std::find(adjacentElements.begin(), adjacentElements.end(), element) == adjacentElements.end())
            {
                adjacentElements.push_back(element);               
            }*/
            hed = hed->m_edge->getTwin(hed->m_id)->m_heNext;
        }
   
    }
    //m_adjacentFaces = adjacentElements;
    adjacentElements.push_back(this);
    // erase duplicates
    std::set<Face*> s(adjacentElements.begin(), adjacentElements.end());
    adjacentElements.assign(s.begin(), s.end());
    return adjacentElements;
}

Point Face::getNormalUnitary()
{
    return Point(m_normal / m_jacobian);
}

void Face::setNormal()
{
    Point v1 = m_points[1]->m_coord - m_points[0]->m_coord;
    Point v2 = m_points[2]->m_coord - m_points[0]->m_coord;
    m_normal = Point::crossProd(v2, v1);
    m_jacobian = m_normal.getNorm();
}

void Face::setFMMSizes(const int& N)
{

    m_MEG = new std::complex<double>[N];
    m_MEH = new std::complex<double>[N];
    m_N = N;

}

Point Face::getYc()
{
    // Calculate centroid
    return Point((m_points[0]->m_coord + m_points[1]->m_coord + m_points[2]->m_coord)/3.0);
}

void Face::getLeafElements(std::vector<Face*> * leafElements )
{
  if (m_isLeaf)
  {
    leafElements->push_back(this);
    return;
  }
  else 
  {
    for (Face* child : m_childrenElement)
    {

        child->getLeafElements(leafElements);

    }

  }

    
}
void Face::getLevelElements(std::vector<Face*>* leafElements, short int level)
{
  if (m_nL == level)
  {
    leafElements->push_back(this);
    return;
  }
  else
  {
    for (Face* child : m_childrenElement)
    {

      child->getLevelElements(leafElements, level);

    }

  }


}

void Face::accumulate()
{
  for (Face* child : m_childrenElement)
  {
    for (int k = 0; k < m_N; k++)
    {
      m_MEG[k] += child->m_MEG[k];
      m_MEH[k] = m_MEH[k] + child->m_MEH[k];
    }
    child->reset();
  }
}

void Face::reset()
{
  //delete[] m_MEG;
  //delete[] m_MEH;
  for (int i = 0; i < m_N; i++)
  {
    m_MEG[i] = 0.0;
    m_MEH[i] = 0.0;
  }
}

std::vector<Vertex*> Face::getBoundaryPoints()
{
  std::vector<Vertex*> boundaryPoints;
  std::vector<Edge*> boundaryEdges;

  std::vector<Face*> adjacentFaces = getAdjacentElements();


  adjacentFaces.push_back(this);
  // find first Element
  HalfEdge* firstHed = nullptr;
  HalfEdge* secondHed = nullptr;
  HalfEdge* hedTwin = nullptr;
  Face* element;
  bool flag = false;
  bool mark;
  for (Face* adj : adjacentFaces)
  {
    for (HalfEdge* hed : adj->m_heds)
    {
      hedTwin = hed->m_edge->getTwin(hed->m_id);
      element = hedTwin->m_el;
      mark = !(std::find(adjacentFaces.begin(), adjacentFaces.end(), element) == adjacentFaces.end());
      if (!(mark))
      {
        firstHed = hed;
        flag = true;
        break;
      }
    }
    if (flag)
      break;
  }
  if (!flag)
  {
    return boundaryPoints;
  }
  Vertex* v0;
  Vertex* v1 = nullptr;
  v0 = firstHed->getP0();
  v0->markTemp = true;
  boundaryPoints.push_back(v0);
  boundaryEdges.push_back(firstHed->m_edge);
  while (v0 != v1)
  {
    secondHed = firstHed->m_heNext;
    element = secondHed->m_el;
    mark = !(std::find (adjacentFaces.begin(), adjacentFaces.end(), element ) == adjacentFaces.end());
    while (mark)
    {
      firstHed = secondHed;
      hedTwin = firstHed->m_edge->getTwin(firstHed->m_id);
      secondHed = hedTwin->m_heNext;
      element = secondHed->m_el;
      mark = !(std::find(adjacentFaces.begin(), adjacentFaces.end(), element) == adjacentFaces.end());
    }
    v1 = firstHed->getP0();
    boundaryPoints.push_back(v1);
    v1->markTemp = true;
    boundaryEdges.push_back(firstHed->m_edge);
  }
  // get leaf nodes
  for (Edge* edge : boundaryEdges)
  {
    edge->getLeafNodes(&boundaryPoints);
  }





  return boundaryPoints;
}

void Face::showMEG()
{
  std::cout << m_hedInc << "\n";
  std::cout << "MEG" << "\n";
  std::cout.precision(15);
  for (int i = 0; i < m_N; i++)
  {
    std::cout << std::imag(m_MEG[i]) << "\n";

  }
}

void Face::showMEH()
{
  std::cout << m_hedInc << "\n";
  std::cout << "MEH" << "\n";
  std::cout.precision(15);
  for (int i = 0; i < m_N; i++)
  {
    std::cout << std::real(m_MEH[i]) << "\n";

  }
}

void Face::set_FMM_ME_sizes()
{
  m_MEG_p.resize(3);
  m_MEH_p.resize(3);
  for (int i = 0; i < 3; i++)
  {
    m_MEH_p[i] = new std::complex<double>[m_N];
    m_MEG_p[i] = new std::complex<double>[m_N];
  }
}

Matrix Face::exactIntG(Vertex* source)
{
  double jac2 = m_jacobian * m_jacobian;
  double tol = 1.0e-10;
  Point u = (m_points[0]->m_coord - m_points[2]->m_coord);
  Point v = (m_points[1]->m_coord - m_points[2]->m_coord);
  Point x0 = source->m_coord - m_points[2]->m_coord;
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
  Matrix out(3, 1);
  if (std::abs(d00) > tol)
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

        I[{0, 0}] = std::real((-sqrtb24ac4ad * logB2SqrtAC + (sqrtb24ac4ad)*log2ABSQRTABC
          - sqrtAD * log(-b24ac + 4. * a * sqrtCD - b * (sqrtb24ac4ad)) +
          sqrtAD * log(-b24ac + 4. * a * sqrtDABC - 2. * a * (sqrtb24ac4ad)-b * (sqrtb24ac4ad)) +
          sqrtAD * log(-b24ac + 4. * a * sqrtCD + b * (sqrtb24ac4ad)) -
          sqrtAD * log(-b24ac + 4. * a * sqrtDABC + 2. * a * (sqrtb24ac4ad)+b * (sqrtb24ac4ad))) / (sqrtA * sqrtb24ac4ad));

        I[{1, 0}] = std::real((b24ac * sqrtCMD * sqrtABCMD *
          logB2SqrtAC - b24ac * sqrtCMD * sqrtABCMD * log2ABSQRTABC +
          sqrtA * d * (2. * b * sqrtABCMD * logSqrtCSqrtCMD -
            2. * (2. * a + b) * sqrtCMD * logABCABCMD +
            (2. * a * sqrtCMD + b * sqrtCMD - b * sqrtABCMD) * logD))
          / (2. * sqrtA * (-b24ac - 4. * a * d) * sqrtCMD * sqrtABCMD));

        I[{2, 0}] = std::real((-b * b24ac * sqrtCMD * sqrtABCMD *
          logB2SqrtAC - b * b24ac * sqrtCMD * sqrtABCMD * log2ABSQRTABC +
          2. * sqrtA * ((sqrtC - sqrtABC) * sqrtABCMD
            * sqrtb24ac4ad * sqrtb24ac4ad + a * d * (-4. * sqrtCMD * sqrtABCMD * logSqrtCSqrtCMD + (b + 2. * c + -2. * d)
              * (2. * logABCABCMD - logD) + 2. * sqrtCMD * sqrtABCMD * logD))) /
          (4. * sqrtA * a * (-b24ac - 4. * a * d) * sqrtABCMD));
       

        out += Ns[i] * I * jacI;

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

        I[{0, 0}] = std::real((log2ABSQRTABC - logB2SqrtAC) / sqrtA);
        I[{1, 0}] = I[{0, 0}] * 0.5;
        I[{2, 0 }] = std::real((-2. * sqrtAC + 2. * sqrtAABC - b * sqrtA * I[{0, 0}]) / (4. * sqrtA * a));

        out += Ns[i] * I * jacI;
      }
    }
  }

  return out;
}

Matrix Face::exactIntH(Vertex* source)
{
  double jac2 = m_jacobian * m_jacobian;
  double tol = 1.0e-10;
  Point u = (m_points[0]->m_coord - m_points[2]->m_coord);
  Point v = (m_points[1]->m_coord - m_points[2]->m_coord);
  Point x0 = source->m_coord - m_points[2]->m_coord;
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
  Matrix out(3, 1);
  if (std::abs(d00) > tol)
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

        I[{0, 0}] = std::real((-log(4. * a * (c + sqrtCD) + b * (-b + sqrtb24ac4ad))
          + log(4. * a * (c + sqrtCD) - b * (b + sqrtb24ac4ad))
          + -log(-b * (b + sqrtb24ac4ad) - 2. * a * (-2. * c - 2. * sqrtABC *
            sqrtD + sqrtb24ac4ad)) + log(b * (-1. * b + sqrtb24ac4ad) + 2. * a * (2. * c + 2. * sqrtDABC
              + sqrtb24ac4ad))) / (sqrtD * sqrtB24ADMC));
        I[{1, 0}] = std::real((4. * sqrtA * (logB2SqrtAC - log2ABSQRTABC) -
          2. * b * logSqrtCSqrtCMD / sqrtCMD +
          2. * (2. * a + b) * logABCABCMD / sqrtABCMD +
          logD * (b * (1. / sqrtCMD - 1. / sqrtABCMD) - (2. * a / sqrtABCMD))) /
          -B24ADMC);
        I[{2, 0}] = std::real(((2. * b * (log2ABSQRTABC - logB2SqrtAC) / sqrtA) +
          4. * sqrtCMD * logSqrtCSqrtCMD - 2. * sqrtCMD * logD -
          ((b + 2. * c - 2. * d) * (2. * logABCABCMD - logD) / sqrtABCMD)) / -B24ADMC);

        out += Ns[i] * I * -jacI;

      }
    }
  }
  out = out * d00;

  return out;
}

std::vector<Face*> Face::getChildrenElement()
{
  return std::vector<Face*> {m_childrenElement[0], m_childrenElement[1], m_childrenElement[2], m_childrenElement[3]};
}

double Face::getRadius()
{
    Point L1 = m_points[1]->m_coord - m_points[0]->m_coord;
    Point L2 = m_points[2]->m_coord - m_points[0]->m_coord;
    Point L3 = m_points[2]->m_coord - m_points[1]->m_coord;
    double L12 = L1 * L1;
    double L22 = L2 * L2;
    double L32 = L3 * L3;
    if (L12 >= L22)
    {
        if (L12 >= L32)
        {
            return L12;
        }
        else
        {
            return L32;
        }
    }
    else
    {
        if (L22 >= L32)
        {
            return L22;
        }
        else
        {
            return L32;
        }
    }
}

Face::~Face()
{
  // reset is doing the job
    //delete[] m_MEG;
    //delete[] m_MEH;
}






