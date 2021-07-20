#include <cmath>
#include <stdexcept>
#include <iostream>
#include "matrix.h"

Matrix::Matrix(int r, int c) : mNumRows(r), mNumColumns(c) {
	mData = std::vector<double>(r * c, 0.0);
}

Matrix::Matrix(int dim) : mNumRows(dim), mNumColumns(dim) {
	mData = std::vector<double>(dim * dim, 0.0);
}

Matrix::Matrix(int r, int c, std::vector<double> data) :
	mNumRows(r), mNumColumns(c), mData(data) {

}

Matrix::~Matrix() {
	mData.clear();
}

double& Matrix::operator[](std::pair<int, int> ij) {
	const int index = ij.first * mNumColumns + ij.second;
	return mData.at(index);
}

Matrix Matrix::operator*(Matrix& other) {
	Matrix res{ mNumRows, other.mNumColumns };
	for (int r = 0; r < res.mNumRows; ++r) {
		for (int c = 0; c < res.mNumColumns; ++c) {
			for (int n = 0; n < this->mNumColumns; ++n) {
				const double rowVal = this->operator[]({ r, n });
				const double colVal = other[{n, c}];
				const double prod = rowVal * colVal;
				res[{r, c}] += prod;
			}
		}
	}
	return res;
}

Matrix Matrix::operator*(const double& other) const {
	Matrix copy{ *this };
	for (int i = 0; i < copy.mNumRows; ++i) {
		for (int j = 0; j < copy.mNumColumns; ++j) {
			copy[{i, j}] *= other;
		}
	}
	return copy;
}

Matrix Matrix::operator+(Matrix other) const {
	Matrix copy{ *this };
	for (int i = 0; i < copy.mNumRows; ++i) {
		for (int j = 0; j < copy.mNumColumns; ++j) {
			copy[{i, j}] += other[{i, j}];
		}
	}
	return copy;
}

void Matrix::operator+=(Matrix other) {
	for (int i = 0; i < mNumRows; ++i) {
		for (int j = 0; j < mNumColumns; ++j) {
			this->operator[]({ i, j }) += other[{i, j}];
		}
	}
}

void Matrix::operator-=(Matrix other) {
	for (int i = 0; i < mNumRows; ++i) {
		for (int j = 0; j < mNumColumns; ++j) {
			this->operator[]({ i, j }) -= other[{i, j}];
		}
	}
}

Matrix Matrix::operator-(Matrix other) const {
	Matrix copy{ *this };
	for (int i = 0; i < copy.mNumRows; ++i) {
		for (int j = 0; j < copy.mNumColumns; ++j) {
			copy[{i, j}] -= other[{i, j}];
		}
	}
	return copy;
}

Matrix Matrix::operator-() const {
	Matrix copy{ *this };
	for (int i = 0; i < copy.mNumRows; ++i) {
		for (int j = 0; j < copy.mNumColumns; ++j) {
			copy[{i, j}] = -copy[{i, j}];
		}
	}
	return copy;
}

Matrix Matrix::transpose() {
	Matrix t{ mNumColumns, mNumRows };
	for (int i = 0; i < mNumRows; ++i) {
		for (int j = 0; j < mNumColumns; ++j) {
			t[{j, i}] = this->operator[]({ i, j });
		}
	}
	return t;
}

Matrix Matrix::setOneAtZerosOnDiagonal() {
	Matrix copy{ *this };
	for (int i = 0; i < mNumRows; ++i) {
		for (int j = 0; j < mNumColumns; ++j) {
			if (i == j && mData[i * mNumColumns + j] == 0.0)
				copy[{i, j}] = 1.0;
		}
	}
	return copy;
}

int Matrix::rows() const {
	return mNumRows;
}

int Matrix::columns() const {
	return mNumColumns;
}

void Matrix::forEach(const Matrix::MatrixElem& block) {
	int r = 0;
	int c = 0;
	while (r < mNumRows) {
		while (c < mNumColumns) {
			block(r, c, this->operator[]({ r, c }));
			c++;
		}
		c = 0;
		r++;
	}

}

Matrix Matrix::solveLinearSystemCGS(Matrix& a, Matrix& b, const int maxIt, const double tol) {
	double residue;
	double rho1;
	double rho2;
	double alpha;
	double beta;

	Matrix x = b;

	x[{1, 0}] = 1;

	double normB = b.norm();
	Matrix r = b - a * x;
	Matrix rTilde = r;

	if (normB == 0.0)
		normB = 1.0;

	residue = r.norm() / normB;
	if (residue <= tol) {
		throw std::runtime_error{ "Initial residue is less than tolerance" };
	}

	Matrix u{ r };
	Matrix p{ u };
	Matrix q{ u };
	Matrix pHat{ p };
	Matrix vHat{ p };
	Matrix uHat{ p };
	Matrix qHat{ p };
	for (int i = 1; i <= maxIt; ++i) {
		rho1 = rTilde.dot(r);
		if (rho1 == 0.0) {
			throw std::runtime_error{ "Variable rho1 is equals to zero" };
		}
		if (i == 1) {
			u = r;
			p = u;
		}
		else {
			beta = rho1 / rho2;
			u = r + q * beta;
			p = u + (q + p * beta) * beta;
		}
		pHat = p;
		vHat = a * pHat;

		alpha = rho1 / rTilde.dot(vHat);

		q = u - vHat * alpha;

		uHat = u + q;

		x += uHat * alpha;

		qHat = a * uHat;

		r -= qHat * alpha;

		rho2 = rho1;

		residue = r.norm() / normB;
		if (residue < tol) {
			return x;
		}
	}

	throw std::runtime_error{ "Could not converge to solution" };
}
Matrix Matrix::solveLinearSystemCG(Matrix& a, Matrix& b, const int maxIt, const double tol) {
	double resid;
	Matrix p(2, 1);
	Matrix z = p;
	Matrix q = z;
	double alpha, beta, rho, rho_1;
	//Matrix alpha(2, 1);
	//Matrix	beta(2, 1);
	//Matrix rho(2, 1);
	//Matrix rho_1(2,1);
	double normb = b.norm();

	Matrix x = b;
	Matrix r = b - a * x;

	if (normb == 0.0)
		normb = 1;

	if ((resid = r.norm() / normb) <= tol) {
		throw std::runtime_error{ "Initial residue is less than tolerance" };
	}

	for (int i = 1; i <= maxIt; i++) {
		z = r;
		rho = r.dot(z);
		if (i == 1)
			p = z;
		else {
			beta = rho / rho_1;
			p = z + p * beta;
		}

		q = a * p;
		alpha = rho / q.dot(p);

		x += p * alpha;
		r -= q * alpha;

		resid = r.norm() / normb;

		if ((resid) <= tol) {
			return x;
		}

		rho_1 = rho;
	}

	throw std::runtime_error{ "Could not converge to solution" };

}

double Matrix::norm() {
	double res = 0.0;
	forEach([&](int& r, int& c, double val) {
		res += val * val;
		});
	return sqrt(res);
}

double Matrix::dot(Matrix& matrix) {
	double res = 0.0;
	forEach([&](int& r, int& c, double val) {
		res += val * matrix[{r, c}];
		});
	return res;
}

Matrix Matrix::getSubMatrix(int rowFrom, int rowTo, int colFrom, int colTo) {

	const int m = (rowTo - rowFrom) + 1;
	const int n = (colTo - colFrom) + 1;
	Matrix A(m, n);
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < m; i++)
		{
			const int index = (rowFrom - 1 + i) * mNumColumns + (colFrom - 1 + j);
			A[{i, j}] = mData.at(index);
		}
	}
	return A;
}
void Matrix::setSubMatrix(int rowFrom, int rowTo, int colFrom, int colTo, Matrix & A) {

	const int m = (rowTo - rowFrom) + 1;
	const int n = (colTo - colFrom) + 1;
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < m; i++)
		{
			const int index = (rowFrom - 1 + i) * mNumColumns + (colFrom - 1 + j);
			mData.at(index) = A[{i, j}];
		}
	}
}
void Matrix::show()
{
	std::cout.precision(15);
	for (int i = 0; i < mNumRows; i++)
	{
		std::cout << i << " " << mData.at(i) << "\n";
	}
}
Matrix Matrix::inverse() {
	Matrix invMat = Matrix(rows(), columns());
	for (int i = 0; i < invMat.rows(); i++) {
		Matrix b = Matrix(rows(), 1);
		b[{i, 0}] = 1.0;
		Matrix x = solveLinearSystemCGS(*this, b);
		invMat.forEach([&](int& m, int& n, double val) {
			if (n == i)
				invMat[{m, n}] = x[{m, 0}];
			});
	}
	return invMat;
}
