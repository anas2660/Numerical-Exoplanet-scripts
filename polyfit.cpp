#include "polyfit.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <vector>

// Polyfit implementation gotten from,
// http://www.vilipetek.com/2013/10/17/polynomial-fitting-in-c-not-using-boost/

template <class T> class matrix {
public:
  matrix(unsigned int nRows, unsigned int nCols)
      : m_nRows(nRows), m_nCols(nCols), m_oData(nRows * nCols, 0) {
    if (!nRows || !nCols) {
      throw range_error("invalid matrix size");
    }
  }

  static matrix identity(unsigned int nSize) {
    matrix oResult(nSize, nSize);

    int nCount = 0;
    std::generate(oResult.m_oData.begin(), oResult.m_oData.end(),
                  [&nCount, nSize]() { return !(nCount++ % (nSize + 1)); });

    return oResult;
  }

  inline T &operator()(unsigned int nRow, unsigned int nCol) {
    if (nRow >= m_nRows || nCol >= m_nCols) {
      throw out_of_range("position out of range");
    }

    return m_oData[nCol + m_nCols * nRow];
  }

  inline matrix operator*(matrix &other) {
    if (m_nCols != other.m_nRows) {
      throw domain_error("matrix dimensions are not multiplicable");
    }

    matrix oResult(m_nRows, other.m_nCols);
    for (unsigned int r = 0; r < m_nRows; ++r) {
      for (unsigned int ocol = 0; ocol < other.m_nCols; ++ocol) {
        for (unsigned int c = 0; c < m_nCols; ++c) {
          oResult(r, ocol) += (*this)(r, c) * other(c, ocol);
        }
      }
    }

    return oResult;
  }

  inline matrix transpose() {
    matrix oResult(m_nCols, m_nRows);
    for (unsigned int r = 0; r < m_nRows; ++r) {
      for (unsigned int c = 0; c < m_nCols; ++c) {
        oResult(c, r) += (*this)(r, c);
      }
    }
    return oResult;
  }
  inline unsigned int rows() { return m_nRows; }
  inline unsigned int cols() { return m_nCols; }
  inline vector<T> data() { return m_oData; }
  void print() {
    for (unsigned int r = 0; r < m_nRows; r++) {
      for (unsigned int c = 0; c < m_nCols; c++) {
        std::cout << (*this)(r, c) << "\t";
      }
      std::cout << std::endl;
    }
  }

private:
  std::vector<T> m_oData;

  unsigned int m_nRows;
  unsigned int m_nCols;
};

template <typename T> class Givens {
public:
  Givens() : m_oJ(2, 2), m_oQ(1, 1), m_oR(1, 1) {}
  /*
          Performs QR factorization using Givens rotations.
  */
  void Decompose(matrix<T> &oMatrix) {
    int nRows = oMatrix.rows();
    int nCols = oMatrix.cols();
    if (nRows == nCols) {
      nCols--;
    } else if (nRows < nCols) {
      nCols = nRows - 1;
    }
    m_oQ = matrix<T>::identity(nRows);
    m_oR = oMatrix;
    for (int j = 0; j < nCols; j++) {
      for (int i = j + 1; i < nRows; i++) {
        GivensRotation(m_oR(j, j), m_oR(i, j));
        PreMultiplyGivens(m_oR, j, i);
        PreMultiplyGivens(m_oQ, j, i);
      }
    }

    m_oQ = m_oQ.transpose();
  }

  /*
          Find the solution for a matrix.
          http://en.wikipedia.org/wiki/QR_decomposition#Using_for_solution_to_linear_inverse_problems
  */
  matrix<T> Solve(matrix<T> &oMatrix) {
    matrix<T> oQtM(m_oQ.transpose() * oMatrix);
    int nCols = m_oR.cols();
    matrix<T> oS(1, nCols);
    for (int i = nCols - 1; i >= 0; i--) {
      oS(0, i) = oQtM(i, 0);
      for (int j = i + 1; j < nCols; j++) {
        oS(0, i) -= oS(0, j) * m_oR(i, j);
      }
      oS(0, i) /= m_oR(i, i);
    }

    return oS;
  }
  const matrix<T> &GetQ() { return m_oQ; }
  const matrix<T> &GetR() { return m_oR; }

private:
  /*
          Givens rotation is a rotation in the plane spanned by two coordinates
     axes. http://en.wikipedia.org/wiki/Givens_rotation
  */
  void GivensRotation(T a, T b) {
    T t, s, c;
    if (b == 0) {
      c = (a >= 0) ? 1 : -1;
      s = 0;
    } else if (a == 0) {
      c = 0;
      s = (b >= 0) ? -1 : 1;
    } else if (abs(b) > abs(a)) {
      t = a / b;
      s = -1 / sqrt(1 + t * t);
      c = -s * t;
    } else {
      t = b / a;
      c = 1 / sqrt(1 + t * t);
      s = -c * t;
    }
    m_oJ(0, 0) = c;
    m_oJ(0, 1) = -s;
    m_oJ(1, 0) = s;
    m_oJ(1, 1) = c;
  }
  /*
          Get the premultiplication of a given matrix
          by the Givens rotation.
  */
  void PreMultiplyGivens(matrix<T> &oMatrix, int i, int j) {
    int nRowSize = oMatrix.cols();

    for (int nRow = 0; nRow < nRowSize; nRow++) {
      double nTemp =
          oMatrix(i, nRow) * m_oJ(0, 0) + oMatrix(j, nRow) * m_oJ(0, 1);
      oMatrix(j, nRow) =
          oMatrix(i, nRow) * m_oJ(1, 0) + oMatrix(j, nRow) * m_oJ(1, 1);
      oMatrix(i, nRow) = nTemp;
    }
  }

private:
  matrix<T> m_oQ, m_oR, m_oJ;
};

template <typename T>
vector<T> polyfit(const vector<T> &oX, const vector<T> &oY, int nDegree) {
  if (oX.size() != oY.size())
    throw std::invalid_argument("X and Y vector sizes do not match");

  // more intuative this way
  nDegree++;

  size_t nCount = oX.size();
  matrix<T> oXMatrix(nCount, nDegree);
  matrix<T> oYMatrix(nCount, 1);

  // copy y matrix
  for (size_t i = 0; i < nCount; i++) {
    oYMatrix(i, 0) = oY[i];
  }

  // create the X matrix
  for (size_t nRow = 0; nRow < nCount; nRow++) {
    T nVal = 1.0f;
    for (int nCol = 0; nCol < nDegree; nCol++) {
      oXMatrix(nRow, nCol) = nVal;
      nVal *= oX[nRow];
    }
  }

  // transpose X matrix
  matrix<T> oXtMatrix(oXMatrix.transpose());
  // multiply transposed X matrix with X matrix
  matrix<T> oXtXMatrix(oXtMatrix * oXMatrix);
  // multiply transposed X matrix with Y matrix
  matrix<T> oXtYMatrix(oXtMatrix * oYMatrix);

  Givens<T> oGivens;
  oGivens.Decompose(oXtXMatrix);
  matrix<T> oCoeff = oGivens.Solve(oXtYMatrix);
  // copy the result to coeff
  return oCoeff.data();
}
template vector<double> polyfit<double>(const vector<double> &oX,
                                        const vector<double> &oY, int nDegree);
template vector<float> polyfit<float>(const vector<float> &oX,
                                      const vector<float> &oY, int nDegree);

template <typename T>
vector<T> polyval(const vector<T> &oCoeff, const vector<T> &oX) {
  size_t nCount = oX.size();
  size_t nDegree = oCoeff.size();
  vector<T> oY(nCount);

  for (size_t i = 0; i < nCount; i++) {
    T nY = 0;
    T nXT = 1;
    T nX = oX[i];
    for (size_t j = 0; j < nDegree; j++) {
      // multiply current x by a coefficient
      nY += oCoeff[j] * nXT;
      // power up the X
      nXT *= nX;
    }
    oY[i] = nY;
  }

  return oY;
}

template vector<double> polyval(const vector<double> &oCoeff,
                                const vector<double> &oX);
template vector<float> polyval(const vector<float> &oCoeff,
                               const vector<float> &oX);
