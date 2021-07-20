#include "gmres.h"
#include<iostream>


GMRES::GMRES(FMM* fmm)
{
  m_fmm = fmm;
}

void GMRES::Update(Matrix& x, int k, Matrix& h, Matrix & s, Matrix& Q)
{
  Matrix y = s.getSubMatrix(1,k+1,1,1);
  int n = Q.rows();
  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y[{i, 0}] /= h[{i, i}];
    for (int j = i - 1; j >= 0; j--)
      y[{j, 0}] -= h[{j, i}] * y[{i, 0}];
  }

  Matrix v = Q.getSubMatrix(1, n, 1, k+1);
    x += v * y;
}

Matrix GMRES::Solver(Matrix &b,int m,  int max_iter, double tol)
{
  int i, j = 1, k;
  int n = b.rows();
  double bNorm = b.norm();
  double rNorm;
  Matrix H(m + 1, m);
  Matrix s(n, 1);
  Matrix cs(m + 1,1);
  Matrix sn(m + 1,1);
  Matrix e1(n , 1);
  Matrix Q(n, m + 1);
  e1[{0, 0}] = 1.0;

  Matrix x(n, 1);
  x = b;
 
  //Matrix x = b;
  Matrix r = b;
 
  Matrix v = b;
  Matrix res = b;
  Matrix wH = b;
  Matrix w(m + 1,1);

  while (j <= max_iter) {

    res = m_fmm->matrixVectorMulti(x);


    r = (b - res); 
    rNorm = r.norm();
    Matrix rNormr = r * (1 / rNorm);
    Q.setSubMatrix(1,n,1,1, rNormr);    // ??? r / beta
    s= e1*rNorm;
    for (i = 0; i < m ; i++) 
    {
      v = Q.getSubMatrix(1, n, i+1, i+1);

      w = m_fmm->matrixVectorMulti(v);
      
      for (k = 0; k <= i; k++) {
        v = Q.getSubMatrix(1, n, k+1, k+1);
        H [{ k, i }] = w.dot(v);
        w -=  v * H[{k, i}];
      }
      H[{ i + 1, i }] = w.norm();
      wH = w * (1.0 / H[{i + 1, i}]);
      Q.setSubMatrix(1, n, i + 2, i + 2, wH); // ??? w / H(i+1, i)

      for (k = 0; k < i; k++)
        ApplyPlaneRotation(H[{k, i}], H[{k + 1, i}], cs[{k, 0}], sn[{k, 0}]);

      GeneratePlaneRotation(H[{i, i}], H[{i + 1, i}], cs[{i, 0}], sn[{i, 0}] );
      ApplyPlaneRotation(H[{i, i}], H[{i + 1, i}], cs[{i, 0}], sn[{i, 0}]);
      ApplyPlaneRotation(s[{i, 0}], s[{i + 1, 0}], cs[{i, 0}], sn[{i, 0}]);
      double residBef = m_resid;
      m_resid = std::abs(s[{i + 1,0}]) / b.norm();
      double err = std::abs(m_resid - residBef) / residBef;

      if (m_resid < tol || err < tol  ) {
        Update(x, i, H, s, Q);
        return x;
      }

    }
    Update(x, m - 1, H, s, Q);
    res = m_fmm->matrixVectorMulti(x);
    r = b - res;
    rNorm = r.norm();
    s[{i + 1, 0}] = rNorm;
    m_resid = std::abs(s[{i + 1, 0}]) / b.norm();
    if (m_resid < tol) {
      return x;
    }
    j++;
  }
  return x;
}

void GMRES::GeneratePlaneRotation(double& dx, double& dy, double& cs, double& sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  }
  else if (abs(dy) > abs(dx)) {
    double temp = dx / dy;
    sn = 1.0 / sqrt(1.0 + temp * temp);
    cs = temp * sn;
  }
  else {
    double temp = dy / dx;
    cs = 1.0 / sqrt(1.0 + temp * temp);
    sn = temp * cs;
  }
}

void GMRES::ApplyPlaneRotation(double& dx, double& dy, double& cs, double& sn)
{
  double temp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

Matrix GMRES::SolveM(Matrix& A, Matrix& b, int m, int max_iter, double tol)
{
  double resid;
  int i, j = 1, k;
  int n = b.rows();
  double bNorm = b.norm();
  double rNorm;
  Matrix H(m + 1, m);
  Matrix s(m + 1, 1);
  Matrix cs(m + 1, 1);
  Matrix sn(m + 1, 1);
  Matrix e1(m+1, 1);
  Matrix Q(n, m + 1);
  e1[{0, 0}] = 1.0;
  Matrix x = b;
  Matrix r = b;

  Matrix v = b;
  Matrix wH = b;
  Matrix w(m + 1, 1);

  while (j <= max_iter) {

    r = (b - A*x);
    rNorm = r.norm();
    Matrix rNormr = r * (1 / rNorm);
    Q.setSubMatrix(1, n, 1, 1, rNormr);    // ??? r / beta
    s = e1 * rNorm;

    for (i = 0; i < m; i++)
    {
      v = Q.getSubMatrix(1, n, i + 1, i + 1);
      w = A*v;
      for (k = 0; k <= i; k++) {
        v = Q.getSubMatrix(1, n, k + 1, k + 1);
        H[{ k, i }] = w.dot(v);
        w -= v * H[{k, i}];
      }
      H[{ i + 1, i }] = w.norm();
      wH = w * (1.0 / H[{i + 1, i}]);
      Q.setSubMatrix(1, n, i + 2, i + 2, wH); // ??? w / H(i+1, i)

      for (k = 0; k < i; k++)
        ApplyPlaneRotation(H[{k, i}], H[{k + 1, i}], cs[{k, 0}], sn[{k, 0}]);

      GeneratePlaneRotation(H[{i, i}], H[{i + 1, i}], cs[{i, 0}], sn[{i, 0}]);
      ApplyPlaneRotation(H[{i, i}], H[{i + 1, i}], cs[{i, 0}], sn[{i, 0}]);
      ApplyPlaneRotation(s[{i, 0}], s[{i + 1, 0}], cs[{i, 0}], sn[{i, 0}]);

      resid = std::abs(s[{i + 1, 0}]) / b.norm();
      if (resid < tol) {
        Update(x, i, H, s, Q);
        tol = resid;
        max_iter = j;
        return x;
      }

    }
    Update(x, m - 1, H, s, Q);
    r = b - A*x;
    rNorm = r.norm();
    s[{i + 1, 0}] = rNorm;
    resid = std::abs(s[{i + 1, 0}]) / b.norm();
    if (resid < tol) {
      tol = resid;
      max_iter = j;
      return x;
    }
    j++;

  }
  return x;
}

Matrix GMRES::BiCG(Matrix& b, int m, int max_iter, double tol)
{
  double resid;
  double rho_1, rho_2, alpha, beta, omega;
  Matrix p{ b };
  Matrix phat{ b };
  Matrix s{ b };
  Matrix shat{ b };
  Matrix t{ b };
  Matrix v{ b };

  Matrix x = b;
  Matrix res = b;
  double normb = b.norm();
  res = m_fmm->matrixVectorMulti(x);
  Matrix r = b - res;
  Matrix rtilde = r;

  if (normb == 0.0)
    normb = 1;


  for (int i = 1; i <= max_iter; i++) {
    rho_1 = rtilde.dot(r);
    if (rho_1 == 0) {
      tol = r.norm()/ normb;
      //return 2;
    }
    if (i == 1)
      p = r;
    else {
      beta = (rho_1 / rho_2) * (alpha / omega);
      p = r +  (p -  v*omega)*beta;
    }
    phat = p;
    v = m_fmm->matrixVectorMulti(phat);
    
    alpha = rho_1 / rtilde.dot(v);
    s = r - v*alpha;
    if ((resid = s.norm() / normb) < tol) {
      x +=  phat * alpha;
      tol = resid;
      return x;
    }
    shat = s;
    m_fmm->matrixVectorMulti(shat);
    omega = t.dot(s) / t.dot(t);
    x += phat*alpha + shat*omega;
    r = s - t*omega;

    rho_2 = rho_1;
    if ((resid = r.norm() / normb) < tol) {
      tol = resid;
      max_iter = i;
      return x;
    }
    if (omega == 0) {
      tol = r.norm() / normb;
      return x;
    }
  }

  tol = resid;
  return x;
}

Matrix GMRES::BiCGM(Matrix& A, Matrix& b, int m, int max_iter, double tol)
{
  double resid;
  double rho_1, rho_2, alpha, beta, omega;
  Matrix p{ b };
  Matrix phat{ b };
  Matrix s{ b };
  Matrix shat{ b };
  Matrix t{ b };
  Matrix v{ b };

  Matrix x = b;
  double normb = b.norm();
  Matrix r = b - A*(x);
  Matrix rtilde = r;

  if (normb == 0.0)
    normb = 1;


  for (int i = 1; i <= max_iter; i++) {
    rho_1 = rtilde.dot(r);
    if (rho_1 == 0) {
      tol = r.norm() / normb;
      return x;
    }
    if (i == 1)
      p = r;
    else {
      beta = (rho_1 / rho_2) * (alpha / omega);
      p = r + (p - v * omega) * beta;
    }
    phat = p;
    v = A * (phat);
    alpha = rho_1 / rtilde.dot(v);
    s = r - v * alpha;
    if ((resid = s.norm() / normb) < tol) {
      x += phat * alpha;
      tol = resid;
      return x;
    }
    shat = s;
    t = A*(shat);
    omega = t.dot(s) / t.dot(t);
    x += phat * alpha + shat * omega;
    r = s - t * omega;

    rho_2 = rho_1;
    if ((resid = r.norm() / normb) < tol) {
      tol = resid;
      max_iter = i;
      return x;
    }
    if (omega == 0) {
      tol = r.norm() / normb;
      return x;
    }
  }

  tol = resid;
  return x;
}
