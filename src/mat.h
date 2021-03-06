/* matrix class

Created Jeu 20 jan 2011 17:16:18
Copyright (C) 2011, 2012, 2013      Gilles Villard

This file is part of the hplll Library

The hplll Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The hplll Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the hplll Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */


#ifndef HPLLL_MAT_H
#define HPLLL_MAT_H

#include <sstream>
#include <iostream>


#include "defs.h"
#include "mixed-col.h"


#ifdef NTL_GMP_LIP
#include <NTL/mat_ZZ.h>

using namespace NTL;
#endif

namespace hplll {



#ifdef NTL_GMP_LIP

// Jeu  1 mai 2014 16:06:02 CEST for conversions between mpz_t and NTL
// From : Verifiable Delegation of Computation verifiable-delegation / gmp-utils

void bytes_from_mpz(unsigned char*& outBuf, size_t* l, const mpz_t n)
{
  outBuf = (unsigned char*)mpz_export(outBuf, l, -1, 1, 0, 0, n);
}

void mpz_from_bytes(mpz_t r, const unsigned char* inBuf, size_t l)
{
  mpz_import(r, l, -1, 1, 0, 0, inBuf);
}


void mpztozz(ZZ& r, const mpz_t n)
{
  size_t nBytes;
  /* for this direction, we'll let the bytes_from_mpz function
   * allocate memory for us. */
  unsigned char* bytes = 0;
  bytes_from_mpz(bytes, &nBytes, n);
  ZZFromBytes(r, bytes, nBytes);
  /* the bytes_from_mpz function above allocated memory, so we have
   * to free it here. */
  free(bytes);
}


void zztompz(mpz_t r, const ZZ& n)
{
  size_t nBytes = NumBytes(n);
  unsigned char bytes[nBytes];
  BytesFromZZ(bytes, n, nBytes);
  mpz_from_bytes(r, bytes, nBytes);

}

#endif



// ******************************************************************
// Begin matrix class
// ******************************************************************

// Comme un tableau lin�aire de colonnes

template<class T>
class matrix {

protected:
  int n, d;
  vector<vector<T> >M;


public:
  // An empty matrix
  matrix() : n(0), d(0) {}

  // rows x cols, all elements are initialized with the default constructor of T
  matrix(int rows, int cols) : n(0), d(0) {
    resize(rows, cols);
  }

  /** Sets the dimensions of this matrix, preserving as much as possible of the
      content. The value of new elements is undefined. */

  void clear() {
    n = d = 0;
    M.clear();
  }


  // Met � z�ro
  void resize(int rows, int cols) {

    //int i, j;
    int j;
    n = rows;
    d = cols;
    M.resize(d);

    for (j = 0; j < d; j++) M[j].resize(n);

    //for (j = 0; j < d; j++) // Commented pb no 0 for both Z_NR and FP_NR
    //for (i = 0; i < n; i++) //  Jeu 13 oct 2016 16:36:38 CEST
    //M[j][i] = 0;

  }

  /** Returns the number of rows */
  int getRows() const {
    return n;
  }

  /** Returns the number of columns */
  int getCols() const {
    return d;
  }

  // T& operator()(size_t i, size_t j): returns entry i,j, _rep[i*m+j]
  //- const T& operator()(size_t i, size_t j) const: same thing on a constant matrix

  // a reference to the element (i, j). */
  inline T& operator()(int i, int j) {
    return M[j][i];
  }

  //
  inline  T& get(const int i, const int j)                {return M[j][i];}

  // Dummy function for being generic with MatrixExp
  inline  T& get_non_normalized(const int i, const int j)                {return M[j][i];}

  //
  inline  void set(const int i, const int j, const T a)                {M[j][i] = a;}

  //
  inline  void set(const int i, const int j, const double a)                {M[j][i] = a;}

  // Returns a constant reference to the element (i, j) on constant objects. */
  inline const T& operator()(int i, int j) const {
    return M[j][i];
  }


  //
  inline  T* getcol(const int j, const int i = 0)                {return &M[j][i];}

  const inline T* getcol(const int j, const int i = 0)  const   {return &M[j][i];}

  // Assumed to replace the whole column
  inline  void setcol(const int j, const T* w, const int l)   {

    T* v = getcol(j);

    for (int i = 0; i < l; i++) v[i] = w[i];
  }

// Assumed to replace the whole column
  inline  void setcol(const int j, const T* w, const int beg, const int l)   {

    T* v = getcol(j);

    for (int i = beg; i < beg + l; i++) v[i] = w[i];
  }

  // Assumed to replace the whole column
  inline void setcol(int j, Z_NR<mpz_t>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      set_z(M[j][i], b[i]);

  }

  // Assumed to replace the whole column
  inline void setcol(int j, Z_NR<long>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      set_z(M[j][i], b[i]);

  }

  // Assumed to replace the whole column
  inline void setcol(int j, Z_NR<__int128_t>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      set_z(M[j][i], b[i]);

  }

  // Assumed to replace the whole column
  inline void setcol(int j, Z_NR<double>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      M[j][i].get_data() = b[i].get_data();

  }

  // Assumed to replace the whole column
  inline void setcol(int j, Z_NR<long double>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      M[j][i].get_data() = b[i].get_data();

  }

  // Assumed to replace the whole column
  inline void setcol(int j, FP_NR<mpfr_t>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      set_mpfr(M[j][i], b[i]);

  }

  // Assumed to replace the whole column
  // From mixed matrices
  inline void setcol(int j, const mixed_col<FP_NR<mpfr_t>, Z_NR<mpz_t> > c, int beg, int l) {

    int dR;
    dR = c.dimRT;
    //int mZ;   // Could be used for some check on the value of l
    //mZ=c.dimZT;

    FP_NR<mpfr_t>* vR;
    vR = c.colRT;

    Z_NR<mpz_t>* vZ;
    vZ = c.colZT;

    if (l != 0) {

      // mpfr vector conversion
      for (int i = beg; ((i < beg + l) && (i < dR)); i++)
        set_mpfr(M[j][i], vR[i]);


      // mpz vector conversion
      for (int i = 0; i < l + beg - dR; i++)
        set_z(M[j][i + dR], vZ[i]);

    } // non zero length

  };





  // Assumed to replace the whole column
  inline void setcol(int j, FP_NR<double>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      (M[j][i]).set(b[i]);

  }

  // Assumed to replace the whole column
  inline void setcol(int j, FP_NR<long double>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      (M[j][i]).set(b[i]);

  }

  // Assumed to replace the whole column
  inline void setcol(int j, FP_NR<dpe_t>* b, int beg, int l) {

    for (int i = beg; i < beg + l; i++)
      (M[j][i]).set(b[i]);

  }

  // Swaps col j and col k
  inline void colswap(const int j, const int k)  {     // j < k

    M[j].swap(M[k]);

  }

  // addcol col j := colj + col k

  inline void addcol(const int j, const int k, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = l - 1; i >= 0; i--) {
      v[i].add(v[i], w[i]);
    }

  }

  inline void paddcol(const int j, const int k, int ibeg, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = ibeg; i < l + ibeg; i++) {
      v[i].add(v[i], w[i]);
    }

  }

  // subcol col j := colj - col k

  inline void subcol(const int j, const int k, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = 0; i < l; i++) {
      v[i].sub(v[i], w[i]);
    }


  }

  inline void psubcol(const int j, const int k, const int ibeg, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = ibeg; i < l + ibeg; i++) {
      v[i].sub(v[i], w[i]);
    }


  }

  // subcol col j := colj - a* col k
  inline void submulcol(const int j, const int k,  const T a, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = 0; i < l; i++)
      v[i].submul(a, w[i]);

  }

  inline void psubmulcol(const int j, const int k,  const T a, const int ibeg, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = ibeg; i < l + ibeg; i++)
      v[i].submul(a, w[i]);

  }

// subcol col j := colj - a* col k
  inline void addmulcol_si(const int j, const int k,  const long a, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = 0; i < l; i++)
      v[i].addmul_si(w[i], a);

  }

  inline void paddmulcol_si(const int j, const int k,  const long a, const int ibeg, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = ibeg; i < l + ibeg; i++)
      v[i].addmul_si(w[i], a);

  }


// subcol col j := colj - a* col k
  inline void addmulcol_si_2exp(const int j, const int k,  const long a, const long expo, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    T tmp;


    for (int i = 0; i < l; i++) {

      tmp.mul_si(w[i], a);
      tmp.mul_2si(tmp, expo);

      v[i].add(v[i], tmp);
    }

  }

  inline void paddmulcol_si_2exp(const int j, const int k,  const long a, const long expo, const int ibeg, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    T tmp;


    for (int i = ibeg; i < l + ibeg; i++) {

      tmp.mul_si(w[i], a);
      tmp.mul_2si(tmp, expo);

      v[i].add(v[i], tmp);
    }

  }

// addmulcol col j := colj + a* col k
  inline void addmulcol(const int j, const int k,  const T a, const int l)  {     // j < k

    T* v = getcol(j);
    T* w = getcol(k);

    for (int i = 0; i < l; i++) v[i].addmul(a, w[i]);

  }

  //ICI 29 mars 2013
  //   v := w - z*a
  inline void fmasub(const int j, const int k,  const T* w, const T* z, const T a, const int l)  {     // j < k

    T* v = getcol(j) + k; // Check


    for (int i = 0; i < l; i++)  {

      v[i] = w[i];
      v[i].submul(z[i], a);
    }

  }

  // Vector operation :  div, v := w/a
  inline void div(const int j, const int k, T* w,  const T a, const int l) {

    T* v = getcol(j) + k;

    for (int i = 0; i < l; i++)  v[i].div(w[i], a);

  };


  inline  void rowswap(const int i, const int k)   {

    T tmp;
    for (int j = 0; j < d; j++) {
      tmp = M[j][i];
      M[j][i] = M[j][k];
      M[j][k] = tmp;

    }
  }

  inline  void addrow(const int i, const int k)   {

    for (int j = 0; j < d; j++)
      M[j][i].add(M[j][i], M[j][k]);
  }

  inline  void subrow(const int i, const int k)   {

    for (int j = 0; j < d; j++)
      M[j][i].sub(M[j][i], M[j][k]);
  }

  inline  void addmulrow(const int i, const int k, const T a)   {

    for (int j = 0; j < d; j++)
      M[j][i].addmul(M[j][k], a);
  }


  inline  void shift(const int lalpha, const int nrows) {

    for (int j = 0; j < d; j++)
      for (int i = 0; i < nrows; i++)
        M[j][i].mul_2si(M[j][i], lalpha);

  }

  inline void set(Matrix<T>& A)
  {

    for (int i = 0; i < n; i++)
      for (int j = 0; j < d; j++)
        M[j][i] = A(i, j);

  }


}; // End matrix class


template<class T> void transpose(Matrix<T>& B, Matrix<T> A)
{

  int m, n, i, j;

  m = A.get_rows();
  n = A.get_cols();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      B(j, i) = A(i, j);


};

template<class T> Matrix<T> transpose(Matrix<T> A)
{

  int m, n, i, j;

  m = A.get_rows();
  n = A.get_cols();

  Matrix<T> B;
  B.resize(n, m);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      B(j, i) = A(i, j);

  return B;

};


// **************
// Block accesses
// **************

// SQUARE HERE, bdim divise n
// Cas rectangle ???

template<class T> inline  ZZ_mat<T>  getblock(matrix<Z_NR<T> > A, const int ii, const int jj, const int nbb, const int diagdec)   {

  ZZ_mat<T> B;

  int n = A.getRows();

  long bdim = n / nbb;

  B.resize(bdim, bdim);

  long  si;
  long sj = jj * bdim;

  for (int j = 0; j < bdim; j++) {
    si = ii * bdim;
    for (int i = 0; i < bdim; i++) {
      B(i, j) = A(si + diagdec, sj + diagdec);
      si += 1;
    }
    sj += 1;
  }

  return (B);
};

template<class T> inline  ZZ_mat<T>  getblock(matrix<FP_NR<T> > A, const int ii, const int jj, const int nbb, const int diagdec)   {

  ZZ_mat<T> B;

  int n = A.getRows();

  long bdim = n / nbb;

  B.resize(bdim, bdim);

  long  si;
  long sj = jj * bdim;

  for (int j = 0; j < bdim; j++) {
    si = ii * bdim;
    for (int i = 0; i < bdim; i++) {
      B(i, j).get_data() = A(si + diagdec, sj + diagdec).get_data(); // TO CHECK IN GENERAL
      si += 1;
    }
    sj += 1;
  }

  return (B);
};


// ****************************************************
// Get
// Block i,j for nbr row blocks and nbc column blocks
//   numbers greater than the dimensions
// ****************************************************

template<class T> inline  matrix<Z_NR<T> >  getblock(matrix<Z_NR<T> > A, const int ii, const int jj,
    const int nbr, const int nbc, const int diagdec)   {

  matrix<Z_NR<T> > B;


  int m = A.getRows();
  int n = A.getCols();

  // Block dimensions

  int dr;
  int decr;

  if ((m % nbr) == 0) {
    dr = m / nbr;
    decr = ii * dr;
  }

  else {

    dr = m / nbr;

    int rem;
    rem = m - dr * nbr;

    if (ii < rem) {
      dr += 1;
      decr = ii * dr;
    }
    else
      decr = rem * (dr + 1) + dr * (ii - rem);
  }

  int dc;
  int decc;

  if ((n % nbc) == 0) {
    dc = n / nbc;
    decc = jj * dc;
  }

  else {

    dc = n / nbc;

    int rem;
    rem = n - dc * nbc;

    if (jj < rem) {
      dc += 1;
      decc = jj * dc;
    }
    else
      decc = rem * (dc + 1) + dc * (jj - rem);
  }

  B.resize(dr, dc);

  for (int i = 0; i < dr; i++)
    for (int j = 0; j < dc; j++)
      B(i, j) = A(decr + i, decc + j);

  return (B);
};


template<class T> inline  Matrix<Z_NR<T> >  getblock(Matrix<Z_NR<T> > A, const int ii, const int jj,
    const int nbr, const int nbc, const int diagdec)   {

  Matrix<Z_NR<T> > B;


  int m = A.get_rows();
  int n = A.get_cols();

  // Block dimensions

  int dr;
  int decr;

  if ((m % nbr) == 0) {
    dr = m / nbr;
    decr = ii * dr;
  }

  else {

    dr = m / nbr;

    int rem;
    rem = m - dr * nbr;

    if (ii < rem) {
      dr += 1;
      decr = ii * dr;
    }
    else
      decr = rem * (dr + 1) + dr * (ii - rem);
  }

  int dc;
  int decc;

  if ((n % nbc) == 0) {
    dc = n / nbc;
    decc = jj * dc;
  }

  else {

    dc = n / nbc;

    int rem;
    rem = n - dc * nbc;

    if (jj < rem) {
      dc += 1;
      decc = jj * dc;
    }
    else
      decc = rem * (dc + 1) + dc * (jj - rem);
  }

  B.resize(dr, dc);

  for (int i = 0; i < dr; i++)
    for (int j = 0; j < dc; j++)
      B(i, j) = A(decr + i, decc + j);

  return (B);
};


// ****************************************************
// Put
// Block i,j for nbr row blocks and nbc column blocks
//   numbers greater than the dimensions
// ****************************************************

template<class T> inline  int  putblock(matrix<Z_NR<T> >& B, matrix<Z_NR<T> > A, const int ii, const int jj,
                                        const int nbr, const int nbc, const int diagdec)   {


  int m = B.getRows();
  int n = B.getCols();

  // Block dimensions

  int dr;
  int decr;


  if ((m % nbr) == 0) {
    dr = m / nbr;
    decr = ii * dr;
  }

  else {

    dr = m / nbr;

    int rem;
    rem = m - dr * nbr;

    if (ii < rem) {
      dr += 1;
      decr = ii * dr;
    }
    else
      decr = rem * (dr + 1) + dr * (ii - rem);
  }

  int dc;
  int decc;

  if ((n % nbc) == 0) {
    dc = n / nbc;
    decc = jj * dc;
  }

  else {

    dc = n / nbc;

    int rem;
    rem = n - dc * nbc;

    if (jj < rem) {
      dc += 1;
      decc = jj * dc;
    }
    else
      decc = rem * (dc + 1) + dc * (jj - rem);
  }


  for (int i = 0; i < dr; i++)
    for (int j = 0; j < dc; j++)
      B(decr + i, decc + j) = A(i, j);


  return (0);
};




// SQUARE HERE, bdim divise n
// Cas rectangle ???

template<class T> inline  int putblock(Matrix<T>& A, Matrix<T> B, const int ii, const int jj, const int nbb, const int diagdec)   {


  int n = A.get_rows();

  long bdim = n / nbb;

  long  si;
  long sj = jj * bdim;

  for (int j = 0; j < bdim; j++) {
    si = ii * bdim;
    for (int i = 0; i < bdim; i++) {
      A(si + diagdec, sj + diagdec) = B(i, j);
      si += 1;
    }
    sj += 1;
  }

  return (0);
};

template<class T> inline  int putblock(matrix<T>& A, Matrix<T> B, const int ii, const int jj, const int nbb, const int diagdec)   {


  int n = A.getRows();

  long bdim = n / nbb;

  long  si;
  long sj = jj * bdim;

  for (int j = 0; j < bdim; j++) {
    si = ii * bdim;
    for (int i = 0; i < bdim; i++) {
      A(si + diagdec, sj + diagdec) = B(i, j);
      si += 1;
    }
    sj += 1;
  }

  return (0);
};

template<class T> inline  int putblock(matrix<T>& A, matrix<T> B, const int ii, const int jj, const int nbb, const int diagdec)   {


  int n = A.getRows();

  long bdim = n / nbb;

  long  si;
  long sj = jj * bdim;

  for (int j = 0; j < bdim; j++) {
    si = ii * bdim;
    for (int i = 0; i < bdim; i++) {
      A(si + diagdec, sj + diagdec) = B(i, j);
      si += 1;
    }
    sj += 1;
  }

  return (0);
};


// ******************************************************************
// matrix operations
// ******************************************************************

// Norme_2, flottant cf srqt
template<class T> inline void fp_norm(T&nn, const T* v, const int n)
{

  nn.mul(v[0], v[0]);

  for (int i = 1; i < n; i++)
    nn.addmul(v[i], v[i]);
  nn.sqrt(nn);

};

// Squared Norm_2, flottant cf srqt
// avec calcul de norme partielle pour le test de Lovasz
template<class T> inline void fp_norm_sq(T&nn, const T* v, const int n)
{

  nn.mul(v[0], v[0]);
  for (int i = 1; i < n; i++)  nn.addmul(v[i], v[i]);

};

// Scalar product
template<class T> inline void scalarprod(T&nn, const T* v, const T* w, const int n)
{

  nn.mul(v[0], w[0]);

  for (int i = 1; i < n; i++)
    nn.addmul(v[i], w[i]);
};

// Vector operation :  addmul, v := v + a w
template<class T> inline void vector_addmul(T* v, const T* w, const T a, const int n)
{
  for (int i = 0; i < n; i++)  v[i].addmul(w[i], a);
};

template<class T> inline void vector_addmul_ui(T* v, const T* w, const unsigned long int a, const int n)
{
  for (int i = 0; i < n; i++)  v[i].addmul_ui(w[i], a);
};

// Vector operation :  submul, v := v - a w
template<class T> inline void vector_submul(T* v, const T* w, const T a, const int n)
{
  for (int i = 0; i < n; i++)  v[i].submul(w[i], a);

};

template<class T> inline void vector_submul_ui(T* v, const T* w, const unsigned long int a, const int n)
{
  for (int i = 0; i < n; i++)  v[i].submul_ui(w[i], a);

};

// Vector operation :  div, v := w/a
template<class T> inline void vector_div(T* v,  const T* w, const T a, const int n)
{
  for (int i = 0; i < n; i++)  v[i].div(w[i], a);
};

// Vector operation :  mul, v := w*a
template<class T> inline void vector_mul(T* v,  const T* w, const T a, const int n)
{
  for (int i = 0; i < n; i++)  v[i].mul(w[i], a);
};

// Vector operation :  sub, v := w -z
template<class T> inline void vector_sub(T* v,  const T* w, const T* z, const int n)
{
  for (int i = 0; i < n; i++)  v[i].sub(w[i], z[i]);
};

// Vector operation :  fmasub, v := w - z*a
template<class T> inline void vector_fmasub(T* v,  const T* w, const T* z, const T a, const int n)
{

  for (int i = 0; i < n; i++)  {
    v[i] = w[i];
    v[i].submul(z[i], a);
  }


};


// ******************************************************************
// ******************************************************************

template<class T> inline  void shift_in(Matrix<T>& B, const int lalpha, const int nrows) {

  for (int j = 0; j < B.get_cols(); j++)
    for (int i = 0; i < nrows; i++)
      (B(i, j)).mul_2si(B(i, j), lalpha);

}

template<class T> void print2maple(matrix<T> B, int n, int d)
{
  // Stockage en lignes
  cout << "Matrix([";
  for (int i = 0; i < n; i++) {
    cout << "[";
    for (int j = 0; j < d - 1; j++) {
      cout << (B.get(i, j));
      cout << ", ";
    }
    cout << B.get(i, d - 1);
    if (i < n - 1) cout << "  ],\n"; else  cout << "  ]\n";
  }
  cout << "]);" << endl;

};

// Matrix assignment
// -----------------

// Identity partially on the rows when nbcol >= nbrows

template<class T> void setId(Matrix<T>& A)
{

  int m, n, i, j;

  m = A.get_rows();
  n = A.get_cols();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A(i, j) = 0;

  T one;
  one = 1;

  for (i = 0; i < m; i++)
    A(i, i) = one;

};

template<class T> void setId(matrix<T>& A)
{

  int m, n, i, j;

  m = A.getRows();
  n = A.getCols();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A(i, j) = 0;

  T one;
  one = 1;

  for (i = 0; i < m; i++)
    A(i, i) = one;

};


template<class T> bool isId(Matrix<T> A)
{

  bool answer = 1;

  int m, n, i, j;

  m = A.get_rows();
  n = A.get_cols();

  T one, zero;
  one = 1;
  zero = 0;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {

      if (i == j) {
        if (A(i, j).cmp(one) != 0) answer = 0;
      }
      else {
        if (A(i, j).cmp(zero) != 0) answer = 0;
      }
    }

  return answer;

};


template<class T> void set(Matrix<T>& B, matrix<T> A)
{

  int m, n, i, j;

  m = A.getRows();
  n = A.getCols();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      B(i, j) = A.get(i, j);

};

template<class T> void set(matrix<T>& B, matrix<T> A)
{

  int m, n, i, j;

  m = A.getRows();
  n = A.getCols();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      B.set(i, j, A.get(i, j));

};

template<class T> void set(matrix<T>& B, Matrix<T> A)
{

  int m, n, i, j;

  m = A.getRows();
  n = A.getCols();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      B(i, j) = A(i, j);

};

template<class T> void set(Matrix<T>& B, Matrix<T> A)
{

  int m, n, i, j;

  m = B.get_rows();
  n = B.get_cols();

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      B(i, j) = A(i, j);

};

// On n rows, for d columns
template<class T> void set(Matrix<T>& B, Matrix<T> A, int n, int d)
{

  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < d; j++)
      B(i, j) = A(i, j);

};

// On n rows, after d1 columns, for d2 columns
template<class T> void set(Matrix<T>& B, Matrix<T> A, int n, int d1, int d2)
{

  int i, j;

  for (i = 0; i < n; i++)
    for (j = d1 + 1; j < d2; j++)
      B(i, j) = A(i, j);

};

// ******************************************************************
// Pas en place
// ******************************************************************

template<class T> void matprod(Matrix<T>& C,  Matrix<T> B, Matrix<T> U)
{


  int n, d, dres, i, j, k;

  n = B.get_rows();
  d = B.get_cols();
  dres = U.get_cols();

  for (i = 0; i < n; i++)  {

    for (j = 0; j < dres; j++) {

      C(i, j).mul(B(i, 0), U(0, j));

      for (k = 1; k < d; k++) {
        C(i, j).addmul(B(i, k), U(k, j));
      }
    }
  }


};


template<class T> void matprod_in(Matrix<T>& C, Matrix<T> U)
{

  int m, n, i, j, k;

  m = C.get_rows();
  n = C.get_cols();

  Matrix<T> tmat;
  tmat.resize(m, n);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      tmat(i, j).mul(C(i, 0), U(0, j));
      for (k = 1; k < n; k++) {
        tmat(i, j).addmul(C(i, k), U(k, j));
      }
    }

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C(i, j) = tmat(i, j);
};

inline void matprod_in_int(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> U)
{

  int m, n, i, j, k;

  m = C.get_rows();
  n = C.get_cols();

  ZZ_mat<mpz_t> tmat;
  tmat.resize(m, n);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      tmat(i, j).mul(C(i, 0), U(0, j));
      for (k = 1; k < n; k++) {
        tmat(i, j).addmul(C(i, k), U(k, j));
      }
    }

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C(i, j) = tmat(i, j);
};



inline void pmatprod_in_int(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> U, int S)
{

  int m, n;

  m = C.get_rows();
  n = C.get_cols();


#ifdef _OPENMP
  #pragma omp parallel for shared(C)
#endif

  for (int l = 0; l < S; l++) {

    int mloc;
    mloc = m / S;
    if ((l + 1) <= m % S)
      mloc += 1;

    int ibeg;
    ibeg = (m / S) * l;
    if ((l + 1) <= m % S)
      ibeg += l;
    else
      ibeg += m % S;

    //cout << "l: " << l << "   mloc: " << mloc << "    ibeg: " << ibeg << endl;

    Matrix<Z_NR<mpz_t> > tmat;
    tmat.resize(mloc, n);

    int i, j, k;

    for (i = 0; i < mloc; i++)
      for (j = 0; j < n; j++) {
        tmat(i, j).mul(C(i + ibeg, 0), U(0, j));
        for (k = 1; k < n; k++) {
          tmat(i, j).addmul(C(i + ibeg, k), U(k, j));
        }
      }

    for (i = 0; i < mloc; i++)
      for (j = 0; j < n; j++)
        C(i + ibeg, j) = tmat(i, j);


  } // Parallel loop


};


inline void matprod_in_int(ZZ_mat<mpz_t>& C, ZZ_mat<long int> U)
{

  int m, n, i, j, k;

  m = C.get_rows();
  n = C.get_cols();

  Matrix<Z_NR<mpz_t> > tmat;
  tmat.resize(m, n);


  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      tmat(i, j).mul_si(C(i, 0), U(0, j).get_data());
      for (k = 1; k < n; k++) {
        tmat(i, j).addmul_si(C(i, k), U(k, j).get_data());
      }
    }


  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C(i, j) = tmat(i, j);
};


inline void pmatprod_in_int(ZZ_mat<mpz_t>& C, ZZ_mat<long int> U, int S)
{

  int m, n;

  m = C.get_rows();
  n = C.get_cols();

  int mloc, ibeg;

#ifdef _OPENMP
  #pragma omp parallel for  private(mloc,ibeg) shared(C,U,m,n,S)
#endif

  for (int l = 0; l < S; l++) {


    mloc = m / S;
    if ((l + 1) <= m % S)
      mloc += 1;


    ibeg = (m / S) * l;
    if ((l + 1) <= m % S)
      ibeg += l;
    else
      ibeg += m % S;

    //cout << "l: " << l << "   mloc: " << mloc << "    ibeg: " << ibeg << endl;

    Matrix<Z_NR<mpz_t> > tmat;
    tmat.resize(mloc, n);

    int i, j, k;

    for (i = 0; i < mloc; i++)
      for (j = 0; j < n; j++) {
        tmat(i, j).mul_si(C(i + ibeg, 0), U(0, j).get_data());
        for (k = 1; k < n; k++) {
          tmat(i, j).addmul_si(C(i + ibeg, k), U(k, j).get_data());
        }
      }


    for (i = 0; i < mloc; i++)
      for (j = 0; j < n; j++)
        C(i + ibeg, j) = tmat(i, j);


  } // Parallel loop



};

inline void tpmatprod_in_int(ZZ_mat<mpz_t>& C_in, ZZ_mat<long int> U_in, int S)
{

  int m, n;

  m = C_in.get_rows();
  n = C_in.get_cols();

  ZZ_mat<mpz_t> C;

  C.resize(n, m);
  transpose(C, C_in);

  ZZ_mat<long int> U;

  U.resize(n, n);
  transpose(U, U_in);


  int nloc, ibeg;

#ifdef _OPENMP
  #pragma omp parallel for  private(nloc,ibeg) shared(C,U,m,n,S)
#endif

  for (int l = 0; l < S; l++) {


    nloc = n / S;
    if ((l + 1) <= n % S)
      nloc += 1;


    ibeg = (n / S) * l;
    if ((l + 1) <= n % S)
      ibeg += l;
    else
      ibeg += n % S;

    //cout << "l: " << l << "   mloc: " << mloc << "    ibeg: " << ibeg << endl;

    Matrix<Z_NR<mpz_t> > tmat;
    tmat.resize(nloc, m);

    int i, j, k;

    for (i = 0; i < nloc; i++)
      for (j = 0; j < m; j++) {
        //tmat(i, j).mul_si(C(i + ibeg, 0), U(0, j).get_data());
        tmat(i, j).mul_si(C(0, j), U(i + ibeg, 0).get_data());
        for (k = 1; k < n; k++) {
          //tmat(i, j).addmul_si(C(i + ibeg, k), U(k, j).get_data());
          tmat(i, j).addmul_si(C(k, j), U(i + ibeg, k).get_data());
        }
      }


    for (i = 0; i < nloc; i++)
      for (j = 0; j < m; j++)
        C_in(j, i + ibeg) = tmat(i, j);


  } // Parallel loop



};


// TO OPTIMIZE
inline void matprod_in_int(ZZ_mat<mpz_t>& C, ZZ_mat<__int128_t> U)
{

  int n, i, j;

  n = U.get_rows();

  ZZ_mat<mpz_t> V;
  V.resize(n, n);

  Z_NR<mpz_t> z;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {

      mpz_get_128int(z, U(i, j));
      V(i, j) = z;

    }

  matprod_in(C, V);

};



inline void pmatprod_in_int(ZZ_mat<mpz_t>& C, ZZ_mat<__int128_t> U, int S)
{

  int m, n;

  m = C.get_rows();
  n = C.get_cols();


  ZZ_mat<mpz_t> V;
  V.resize(n, n);

  Z_NR<mpz_t> z;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {

      mpz_get_128int(z, U(i, j));
      V(i, j) = z;

    }


#ifdef _OPENMP
  #pragma omp parallel for  shared(C,V,m,n,S)
#endif

  for (int l = 0; l < S; l++) {

    int mloc;
    mloc = m / S;
    if ((l + 1) <= m % S)
      mloc += 1;

    int ibeg;
    ibeg = (m / S) * l;
    if ((l + 1) <= m % S)
      ibeg += l;
    else
      ibeg += m % S;

    //cout << "l: " << l << "   mloc: " << mloc << "    ibeg: " << ibeg << endl;

    Matrix<Z_NR<mpz_t> > tmat;
    tmat.resize(mloc, n);

    int i, j, k;

    for (i = 0; i < mloc; i++)
      for (j = 0; j < n; j++) {
        tmat(i, j).mul(C(i + ibeg, 0), V(0, j));
        for (k = 1; k < n; k++) {
          tmat(i, j).addmul(C(i + ibeg, k), V(k, j));
        }
      }

    for (i = 0; i < mloc; i++)
      for (j = 0; j < n; j++)
        C(i + ibeg, j) = tmat(i, j);


  } // Parallel loop


};

inline void matprod_in_si(ZZ_mat<mpz_t>& C, ZZ_mat<long int> U)
{

  int m, n, i, j, k;

  m = C.get_rows();
  n = C.get_cols();

  Matrix<Z_NR<mpz_t> > tmat;
  tmat.resize(m, n);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      tmat(i, j).mul_si(C(i, 0), U(0, j).get_data());
      for (k = 1; k < n; k++) {
        tmat(i, j).addmul_si(C(i, k), U(k, j).get_data());
      }
    }

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C(i, j) = tmat(i, j);
};

// TO OPTIMIZE
inline void matprod_in_si(ZZ_mat<mpz_t>& C, ZZ_mat<__int128_t> U)
{

  int n, i, j;

  n = U.get_rows();

  ZZ_mat<mpz_t> V;
  V.resize(n, n);

  Z_NR<mpz_t> z;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {

      mpz_get_128int(z, U(i, j));
      V(i, j) = z;

    }

  matprod_in(C, V);

};


inline void matprod_in_si(matrix<Z_NR<mpz_t> >& C, ZZ_mat<long int> U)
{

  int m, n, i, j, k;

  m = C.getRows();
  n = C.getCols();

  Matrix<Z_NR<mpz_t> > tmat;
  tmat.resize(m, n);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      tmat(i, j).mul_si(C(i, 0), U(0, j).get_data());
      for (k = 1; k < n; k++) {
        tmat(i, j).addmul_si(C(i, k), U(k, j).get_data());
      }
    }

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C(i, j) = tmat(i, j);
};




template<class T> void matprod(matrix<T>& C,  matrix<T> B, matrix<T> U)
{

  int n, d, dres, i, j, k;

  n = B.getRows();
  d = B.getCols();
  dres = U.getCols();

  for (i = 0; i < n; i++)  {

    for (j = 0; j < dres; j++) {

      C(i, j).mul(B(i, 0), U(0, j));


      for (k = 1; k < d; k++) {
        C(i, j).addmul(B(i, k), U(k, j));
      }
    }
  }

};



template<class T> void matprod_in(matrix<T>& C, matrix<T> U)
{

  int m, n, i, j, k;

  m = C.getRows();
  n = C.getCols();

  Matrix<T> tmat;
  tmat.resize(m, n);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      tmat(i, j).mul(C(i, 0), U(0, j));
      for (k = 1; k < n; k++) {
        tmat(i, j).addmul(C(i, k), U(k, j));
      }
    }

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C(i, j) = tmat(i, j);
};


template<class T> void pmatprod_in(matrix<T>& C, Matrix<T> U, int S)  // S blocks
{

  int m, n;

  m = C.getRows();
  n = C.getCols();

  int Sm;

  Sm = m / S;

  int l;

#ifdef _OPENMP
  #pragma omp parallel for shared(C)
#endif

  for (l = 0; l < S; l++) {

    matrix<T> tmat;
    tmat.resize(Sm, n);

    int i, j;

    for (i = 0; i < Sm; i++)
      for (j = 0; j < n; j++)
        tmat.set(i, j, C(l * Sm + i, j));

    matprod_in(tmat, U);

    for (i = 0; i < Sm; i++)
      for (j = 0; j < n; j++)
        C.set(l * Sm + i, j, tmat(i, j));


  } // On the parallel blocks


};


// Product by a block diagonal matrix with S blocks on a
// multiple of S number t of threads
// The column dimension is divisible by 2*S, may not be divisible by t
//   (row dim may not be divisible)

template<class T> void pmaprod_diag_even(matrix<T>& C, Matrix<T> U, int S, int t)  // S blocks
{

  int k;

#ifdef _OPENMP
  #pragma omp parallel for
#endif

  for (k = 0; k < t; k++)  {

    matrix<T> Cb;
    Cb = getblock(C, k / S, k % S, t / S, S, 0);

    Matrix<T> Ub;
    Ub = getblock(U, k % S, k % S, S, S, 0);

    matprod_in(Cb, Ub);

    putblock(C, Cb, k / S, k % S, t / S, S, 0);

  } // On the parallel blocks, k


};


// Product by a block matrix with SxS blocks (square) on a
// multiple of S number t of threads
// The column dimension is divisible by S (not really required), may not be divisible by t
//   (row dim may not be divisible)

template<class T> void pmatprod(matrix<T>& C, Matrix<T> U, int S, int t)  // S blocks
{

  int m, n, k;

  m = C.getRows();
  n = C.getCols();

  matrix<T> L;
  L.resize(m, n * S);



#ifdef _OPENMP
  #pragma omp parallel for shared(C,L,t,S)
#endif

  for (k = 0; k < t; k++)  {

    // Blocks in C, t/S by rows and S by columns
    int I, J;

    I = k / S;
    J = k % S;

    // U in S row slices

    // Times a slice of row J in U

    matrix<T> Cb;
    Cb = getblock(C, I, J, t / S, S, 0);


    Matrix<T> Ub;
    Ub = getblock(U, J, 0, S, 1, 0);


    matprod_in(Cb, Ub);

    // Stored in the J copy in L

    putblock(L, Cb, I, J, t / S, S, 0);


  } // On the parallel blocks, k

#ifdef _OPENMP
  #pragma omp  barrier
#endif

  // Addition of the S copies mxn in L into C

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      C(i, j) = L(i, j);

  for (k = 1; k < S; k++)
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        C(i, j).add(C(i, j), L(i, j + k * n));

};

// Dimensions may change in C

template<class T> void matprod_in(matrix<T>& C, Matrix<T> U)
{

  int m, n, i, j, k;

  m = C.getRows();
  n = C.getCols();

  int nU = U.get_cols();

  Matrix<T> tmat;
  tmat.resize(m, nU);

  for (i = 0; i < m; i++)
    for (j = 0; j < nU; j++) {
      tmat(i, j).mul(C(i, 0), U(0, j));
      for (k = 1; k < n; k++) {
        tmat(i, j).addmul(C(i, k), U(k, j));
      }
    }

  C.resize(m, nU);
  for (i = 0; i < m; i++)
    for (j = 0; j < nU; j++)
      C(i, j) = tmat(i, j);
};


template<class T> void matprod(Matrix<T>& B, Matrix<T> U)
{

  int n, d, i, j, k;

  n = B.get_rows();
  d = B.get_cols();

  Matrix<T> C;
  C.resize(n, d);

  for (i = 0; i < n; i++)  {

    for (j = 0; j < d; j++) {

      C(i, j).mul(B(i, 0), U(0, j));


      for (k = 1; k < d; k++) {
        C(i, j).addmul(B(i, k), U(k, j));
      }
    }
    for (j = 0; j < d; j++)
      B(i, j) = C(i, j);
  }

};




template<class T> void print2maple(Matrix<T> B, int n, int d)
{
  // Stockage en lignes
  cout << "Matrix([";
  for (int i = 0; i < n; i++) {
    cout << "[";
    for (int j = 0; j < d - 1; j++) {
      //B[i][j].print();
      cout << B(i, j);
      cout << ", ";
    }
    //B[i][d-1].print();
    cout << B(i, d - 1);
    if (i < n - 1) cout << "  ],\n"; else  cout << "  ]\n";
  }
  cout << "]);" << endl;

};


// From fplll llldiff
// Up to the sign eg for relation tests

template <class T>
int matcmp (Matrix<T> B1, Matrix<T> B2, int c, int r)
{
  int test = 1, i, j, sg;
  T tmp1;
  T tmp2;

  for (i = 0; i < r; i++) {
    sg = 1;
    tmp1.abs(B1(i, 0));
    tmp2.abs(B2(i, 0));
    if (tmp1.cmp(tmp2) != 0) {
      //cerr << r << ", 0     " << tmp1 << "  " << tmp2 << "\n";
      test = 0;
    }
    if (tmp1.cmp(B1(i, 0)) != 0) sg *= -1;
    if (tmp1.cmp(B2(i, 0)) != 0) sg *= -1;

    if (sg == 1) {
      for (j = 1; j < c; j++) {
        if (B1(i, j).cmp(B2(i, j)) != 0) {
          //cerr << i << " " << j << "     " << B1.Get(i,j) << "  " << B2.Get(i,j) << "\n";
          test = 0;
        }
      }
    }
    else {
      for (j = 1; j < c; j++) {
        tmp1.mul_si(B1(i, j), -1);
        if (tmp1.cmp(B2(i, j)) != 0) {
          //cerr << i << " " << j << "     " << B1.Get(i,j) << "  " << B2.Get(i,j) << "\n";
          test = 0;
        }
      }
    }
  }

  return (test);
}




template<class T> void print2maple(vector<vector<T> > B, int n, int d)
{
  // Stockage en lignes
  cout << "Matrix([";
  for (int i = 0; i < n; i++) {
    cout << "[";
    for (int j = 0; j < d - 1; j++) {
      B[i][j].print();
      cout << ", ";
    }
    B[i][d - 1].print();
    if (i < n - 1) cout << "  ],\n"; else  cout << "  ]\n";
  }
  cout << "]);" << endl;

};

template<class T> void printcol(matrix<T> B, int n, int beg, int k = 0)
{

  cout << "Matrix([";
  for (int i = 0; i < n; i++) {
    cout << "[";
    for (int j = beg; j < beg + k; j++) {
      B(i, j).print();
      cout << ", ";
    }
    B(i, beg + k).print();
    if (i < n - 1) cout << "  ],\n"; else  cout << "  ]\n";
  }
  cout << "]);" << endl;

};

template<class T> void printcol(vector<vector<T> > B, int n, int beg, int k = 0)
{

  cout << "Matrix([";
  for (int i = 0; i < n; i++) {
    cout << "[";
    for (int j = beg; j < beg + k; j++) {
      B[i][j].print();
      cout << ", ";
    }
    B[i][beg + k].print();
    if (i < n - 1) cout << "  ],\n"; else  cout << "  ]\n";
  }
  cout << "]);" << endl;

};



// ********************************************************************
//
//       TRUNCATION OF A BASIS for lift, lehmer, nullspace, L1
//           n = m + d
//           Upper part specific for knapsack or nullspace
//           Lowerpart for the identity or the transformation matrix
//
//  !!! Version  matrix<Z_NR<T> >
//      ==> Also modify ZZ_mat<T>
//          Exactly identical function body
//
// ********************************************************************
// Check the roundings
//
// Division by 2^sigma = 10 ^k the first m=n-d rows
// And truncation, keeping t bits, 2^t = 10 ^ k
// If t=0, keeps everything after the division

/*
template<class T> void trunc_sigma(matrix<Z_NR<T> >& B, ZZ_mat<T> A, long n, long d, long t, long sigma);

template<> void trunc_sigma(matrix<Z_NR<mpz_t> >& B, ZZ_mat<mpz_t> A, long n, long d, long t, long sigma)

{

  int m=n-d;

  // Division of the first m rows
  // and on the fly computation of the log[2] length
  // ------------------------------------------------

  // Computation via pow10 = 10^k = 2^sigma


  mpz_t pow10,ten;
  mpz_init(pow10);
  mpz_init(ten);
  mpz_set_ui(ten, 10);

  long exp10;
  double tmp;

  tmp=((double) sigma)*0.301029995;
  exp10=(long) tmp;
  mpz_pow_ui(pow10, ten, exp10);

  int i,j;

  // Et calcul de la longueur max en base 2

  long max2=0;

  for (i=0; i<m; i++)
    for (j=0; j<d; j++)  {

      mpz_tdiv_q(B(i,j).get_data(), A(i,j).get_data(), pow10);

      max2=max(max2,(long) mpz_sizeinbase(B(i,j).get_data(),2));
    }

  for (i=m; i<n; i++)
    for (j=0; j<d; j++)
      max2=max(max2,(long) mpz_sizeinbase(A(i,j).get_data(),2));


  // On tronque en en gardant msb t
  // ------------------------------

  // To tune with the roundings, here keeps a bit more than t
  // If t=0, keeps everything, TODO: make it clean (do nothing in this case)

  // *****  Negative t
  // An identity part, corresponding transformation from scratch
  if (t< 0) {

    for (i=m; i<n; i++)
      for (j=0; j<d; j++)
  B(i,j)=0;

    for (i=m; i<n; i++)
      B(i,i-m)=1;
  }
  // *****  Zero t
  // One keeps all the bits
  else if (t==0) {

    for (i=m; i<n; i++)
      for (j=0; j<d; j++)
  B(i,j)=A(i,j);

  }
  // **** Positive t
  // One keep t bits
  // Unchanged w.r.t Lehmer
  else {

    tmp=((double) max2-t)*0.301029995;

    exp10= max((long) 0, (long) tmp);
    mpz_pow_ui(pow10, ten, exp10);

    for (i=0; i<m; i++)
      for (j=0; j<d; j++)
  mpz_tdiv_q(B(i,j).get_data(), B(i,j).get_data(), pow10);

    for (i=m; i<n; i++)
      for (j=0; j<d; j++)
  mpz_tdiv_q(B(i,j).get_data(), A(i,j).get_data(), pow10);
  }

  // Check to comment
  //max2=0;
  //for (i=0; i<n; i++)
    //for (j=0; j<d; j++)
      //max2=max(max2,(long) mpz_sizeinbase(B(i,j).get_data(),2));

  //cout << "***** max bit after shift and truncation " << max2 << endl << endl;
  // Compare max2 and t

};
*/



// ********************************************************************
//
//       TRUNCATION OF A BASIS
//
// ********************************************************************
//
// Todo :
//

// Au moins bits en sortie (pour la plus petite colonne)
// si non homog�ne, certaine colonne peuvent rester grandes

// **** If bits = 0 : no truncation


void lift_truncate(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A, long def, long bits) {

  int i, j;

  int n = A.get_rows();
  int d = A.get_cols();


  for (j = 0; j < d; j++)
    C(0, j).mul_2si(A(0, j), def);



  for (i = 1; i < n; i++)
    for (j = 0; j < d; j++)
      C(i, j) = A(i, j);

  // bits > 0 then truncation
  // ------------------------


  if (bits > 0) {

    // Min max norms of the columns
    int mmax, mmin;

    Z_NR<mpz_t> t;

    mmax = 0;
    for (i = 0; i < n; i++) {
      t.abs(C(i, 0));
      mmax = max(mmax, size_in_bits(t));
    }
    mmin = mmax;

    for (j = 1; j < d ; j++) {
      mmax = 0;
      for (i = 0; i < n; i++) {
        t.abs(C(i, j));
        mmax = max(mmax, size_in_bits(t));
      }
      mmin = min(mmin, mmax);
    }

    // Truncation

    if ((mmin - 2) >= bits) {   // We keep a heuristic margin of 2


      long s = bits - mmin + 2;

      for (i = 0; i < n; i++)
        for (j = 0; j < d; j++)
          C(i, j).mul_2si(C(i, j), s);
    }

  } // end bits > 0 / truncation

};



// LONG
void lift_truncate(ZZ_mat<long>& C_out, ZZ_mat<mpz_t> A, long def, long bits) {


  int i, j;

  int n = A.get_rows();
  int d = A.get_cols();

  ZZ_mat<mpz_t> C;
  C.resize(n, d);


  for (j = 0; j < d; j++)
    C(0, j).mul_2si(A(0, j), def);

  for (i = 1; i < n; i++)
    for (j = 0; j < d; j++)
      C(i, j) = A(i, j);



  if (bits > 0) {


    // Min max norms of the columns
    int mmax, mmin, maxloc;

    Z_NR<mpz_t> t;

    mmax = 0;
    for (i = 0; i < n; i++) {
      t.abs(C(i, 0));
      mmax = max(mmax, size_in_bits(t));
    }
    mmin = mmax;


    for (j = 1; j < d ; j++) {
      maxloc = 0;
      for (i = 0; i < n; i++) {
        t.abs(C(i, j));
        maxloc = max(maxloc, size_in_bits(t));
      }
      mmin = min(mmin, maxloc);
      mmax = max(mmax, maxloc);
    }


    // Truncation

    long bbc = mmax; // Bits before conversion

    //cout << "s = " << bits -mmin +2 << endl;
    //cout << "bits = " << bits << endl;
    //cout << "mmin = " << mmin << endl;
    //cout << "mmax = " << mmax << endl;

    if ((mmin - 2) > bits) {

      long s = bits - mmin + 2;
      bbc += s;

//Added by J. Chen on Dec. 21, 2018
      //if (bbc <= 63){


      for (i = 0; i < n; i++)
        for (j = 0; j < d; j++)
          C(i, j).mul_2si(C(i, j), s);


      maxloc = 0;
      for (j = 1; j < d ; j++) {
        for (i = 1; i < n; i++) {
          t.abs(C(i, j));
          maxloc = max(maxloc, size_in_bits(t));
        }
      }

      //cout << "Remains down " << maxloc << endl;
//}// end added by J. Chen
      if (bbc > 63) {

        cout << endl << "** Warning in lift truncate: too many integer digits required for truncation >= 64 ?" << endl;

        //exit(EXIT_FAILURE);
      }

    } // end bits > 0 / truncation

  }

// Conversion to long

  for (i = 0; i < n; i++)
    for (j = 0; j < d; j++)
      C_out(i, j) = mpz_get_si(C(i, j).get_data());

};

// 128 bits
void lift_truncate(ZZ_mat<__int128_t>& C_out, ZZ_mat<mpz_t> A, long def, long bits) {


  int i, j;

  int n = A.get_rows();
  int d = A.get_cols();

  ZZ_mat<mpz_t> C;
  C.resize(n, d);


  for (j = 0; j < d; j++)
    C(0, j).mul_2si(A(0, j), def);

  for (i = 1; i < n; i++)
    for (j = 0; j < d; j++)
      C(i, j) = A(i, j);


  if (bits > 0) {


    // Min max norms of the columns
    int mmax, mmin, maxloc;

    Z_NR<mpz_t> t;

    mmax = 0;
    for (i = 0; i < n; i++) {
      t.abs(C(i, 0));
      mmax = max(mmax, size_in_bits(t));
    }
    mmin = mmax;

    for (j = 1; j < d ; j++) {
      maxloc = 0;
      for (i = 0; i < n; i++) {
        t.abs(C(i, j));
        maxloc = max(maxloc, size_in_bits(t));
      }
      mmin = min(mmin, maxloc);
      mmax = max(mmax, maxloc);
    }

    // Truncation

    long bbc = mmax; // Bits before conversion

    if ((mmin - 2) > bits) {

      long s = bits - mmin + 2;

      bbc += s;

      for (i = 0; i < n; i++)
        for (j = 0; j < d; j++)
          C(i, j).mul_2si(C(i, j), s);


      maxloc = 0;
      for (j = 1; j < d ; j++) {
        for (i = 1; i < n; i++) {
          t.abs(C(i, j));
          maxloc = max(maxloc, size_in_bits(t));
        }

      }

      cout << "Remains down " << maxloc << endl;
    }



    if (bbc > 127) {

      cout << endl << "** Error in lift truncate: too many integer digits required for truncation >= 128 ?" << endl;

      exit(EXIT_FAILURE);

    }


  } // end bits > 0 / truncation


  // Conversion to 128 int

  Z_NR<__int128_t> z;

  for (i = 0; i < n; i++)
    for (j = 0; j < d; j++) {
      mpz_set_128int(z, C(i, j));
      C_out(i, j) = z;
    }
};










// ********************************************************************
//
//       TRUNCATION OF A BASIS for lift, lehmer, nullspace, L1
//           Upper part d rows  n columns
//           Total number of rows m
//             hence lower part m-d rows
//
//       Throw tau bits in the upper part, by division bu 2^tau = 10^k
//         nothing thrown for tau = 0
//       Then global truncation keeping t bits, 2^t = 10^k
//         if t=0 keeps everything
//
//      ==> Also modify matrix<Z_NR<T> >
//          Exactly identical function body
//
// ********************************************************************
//
// Todo : check the roundings
//

template<class ZT, class MatrixZT> void trunc(MatrixZT & B, ZZ_mat<ZT> A, long d, long n, long m, long t, long tau) {

  // Division of the first m rows
  // and on the fly computation of the log[2] length
  // ------------------------------------------------

  // Computation via pow10 = 10^k = 2^tau

  mpz_t pow10, ten;
  mpz_init(pow10);
  mpz_init(ten);
  mpz_set_ui(ten, 10);

  long exp10;
  double tmp;

  int i, j;

  // Division of the upper part
  // and computation of the global max bit length

  long max2 = 0;


  if (tau > 0) {

    tmp = ((double) tau) * 0.301029995;
    exp10 = (long) tmp;
    mpz_pow_ui(pow10, ten, exp10);

    for (i = 0; i < d; i++)
      for (j = 0; j < n; j++)  {
        mpz_tdiv_q(B(i, j).get_data(), A(i, j).get_data(), pow10);
        max2 = max(max2, (long) mpz_sizeinbase(B(i, j).get_data(), 2));
      }
  }
  else {
    for (i = 0; i < d; i++)
      for (j = 0; j < n; j++) {

        B(i, j) = A(i, j);
        //B.set(i,j,A(i,j));
        max2 = max(max2, (long) mpz_sizeinbase(B(i, j).get_data(), 2));
      }
  }

  if (t > 0) {

    for (i = d; i < m; i++)
      for (j = 0; j < n; j++)
        max2 = max(max2, (long) mpz_sizeinbase(A(i, j).get_data(), 2));
  }

  // Truncation that keeps t bits
  // ----------------------------

  // To tune with the roundings, here keeps a bit more than t
  // If t=0, keeps everything

  // *****  Negative t
  // An identity part, corresponding transformation from scratch
  if (t < 0) {

    for (i = d; i < m; i++)
      for (j = 0; j < n; j++)
        B(i, j) = 0;

    for (i = d; i < m; i++)
      B(i, i - d) = 1;
  }
  // *****  Zero t, one keeps evrything
  else if ((t == 0) || (t >= max2)) {

    for (i = d; i < m; i++)
      for (j = 0; j < n; j++)
        B(i, j) = A(i, j);

  }
  // **** Positive t, one keep t bits
  else {

    tmp = ((double) max2 - t) * 0.301029995;
    exp10 = max((long) 0, (long) tmp);
    mpz_pow_ui(pow10, ten, exp10);

    for (i = 0; i < d; i++)
      for (j = 0; j < n; j++)
        mpz_tdiv_q(B(i, j).get_data(), B(i, j).get_data(), pow10);

    for (i = d; i < m; i++)
      for (j = 0; j < n; j++)
        mpz_tdiv_q(B(i, j).get_data(), A(i, j).get_data(), pow10);
  }

};


//template<> void set_f(matrix<Z_NR<mpz_t> >& B, matrix<FP_NR<mpfr_t> > R, long condbits)
void set_f(matrix<Z_NR<mpz_t> >& B, matrix<FP_NR<mpfr_t> > R, long condbits)
{

  int n, d;

  n = B.getRows();
  d = B.getCols();

  FP_NR<mpfr_t> norm, minval;

  int i, j;

  fp_norm(minval, R.getcol(0), n);


  // Avant Mar 29 avr 2014 10:42:12 CEST
  for (j = 1; j < d; j++) {
    fp_norm(norm, R.getcol(j), n);
    if (minval.cmp(norm) > 0) minval = norm;

  }


  Z_NR<mpz_t> z;

  z.set_f(minval);

  //cout << "size : " << size_in_bits(z) << endl;
  //if (size_in_bits(z) > condbits) {
  {
    long s = - (size_in_bits(z) - condbits);

    FP_NR<mpfr_t> bf;

    for (j = 0; j < d; j++)
      for (i = 0; i < n; i++) {

        bf.mul_2si(R(i, j), s);
        B(i, j).set_f(bf);
      }

  }
  /* else { */

  /*   for (j=0; j<d; j++)  */
  /*     for (i=0; i<n; i++)  */


  /*      B(i,j).set_f(R(i,j)); */


  /* }  */
}

//template<> void set_f(matrix<Z_NR<mpz_t> >& B, matrix<FP_NR<mpfr_t> > R, long condbits)
void set_f(ZZ_mat<mpz_t>& B, matrix<FP_NR<mpfr_t> > R, long condbits)
{

  int n, d;

  n = B.get_rows();
  d = B.get_cols();

  FP_NR<mpfr_t> norm, minval;

  int i, j;

  fp_norm(minval, R.getcol(0), n);


  // Avant Mar 29 avr 2014 10:42:12 CEST
  for (j = 1; j < d; j++) {
    fp_norm(norm, R.getcol(j), n);
    if (minval.cmp(norm) > 0) minval = norm;

  }


  Z_NR<mpz_t> z;

  z.set_f(minval);

  //cout << "size : " << size_in_bits(z) << endl;
  //if (size_in_bits(z) > condbits) {
  {
    long s = - (size_in_bits(z) - condbits);

    FP_NR<mpfr_t> bf;

    for (j = 0; j < d; j++)
      for (i = 0; i < n; i++) {

        bf.mul_2si(R(i, j), s);
        B(i, j).set_f(bf);
      }

  }
  /* else { */

  /*   for (j=0; j<d; j++)  */
  /*     for (i=0; i<n; i++)  */


  /*      B(i,j).set_f(R(i,j)); */


  /* }  */
}

/********************************************************/
/* ******         MAXBITSIZE        ******************* */
/********************************************************/

/*
template<class T> int maxbitsize(const ZZ_mat<T>& B);
*/

// Whole matrix
inline  unsigned int maxbitsize(const ZZ_mat<mpz_t>& B) {

  int l = 0;

  int n = B.get_rows();
  int d = B.get_cols();

  for (int i = 0; i < n ; i++)
    for (int j = 0; j < d; j++)
      l = max(l, (int) mpz_sizeinbase(B(i, j).get_data(), 2));

  return l;

}

// Some rows
template<class T>
int maxbitsize(const ZZ_mat<T>& B, int d0, int d, int n) {

  int l = 0;

  for (int i = d0; i < d ; i++)
    for (int j = 0; j < n; j++)
      l = max(l, size_in_bits(B(i, j)));

  return l;

}

// ICI
template<class T>
unsigned int maxbitsize(const matrix<T>& B) {

  int l = 0;

  int n = B.getRows();
  int d = B.getCols();

  for (int i = 0; i < n ; i++)
    for (int j = 0; j < d; j++)
      l = max(l, size_in_bits(B(i, j)));

  return l;

}



// ********************************************************************
//
//       MATRIX CAST
//
// ********************************************************************

// Dimensions may change in C

template<class ZT> void matrix_cast(ZZ_mat<ZT>& B, ZZ_mat<mpz_t> A);

template<class ZT> void matrix_cast(matrix<Z_NR<ZT> > & B,  matrix<Z_NR<mpz_t> >  A);

template<> void matrix_cast(ZZ_mat<mpz_t>& B, ZZ_mat<mpz_t> A) {

  int n = A.get_rows();
  int d = A.get_cols();

  B.resize(n, d);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++)
      B(i, j) = A(i, j);
}

template<> void matrix_cast(matrix<Z_NR<mpz_t> > & B,  matrix<Z_NR<mpz_t> >  A) {

  int n = A.getRows();
  int d = A.getCols();

  B.resize(n, d);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++)
      B(i, j) = A(i, j);
}


template<> void matrix_cast(ZZ_mat<long>& B, ZZ_mat<mpz_t> A) {

  int n = A.get_rows();
  int d = A.get_cols();


  B.resize(n, d);

  if (maxbitsize(A) > ((8 * sizeof(B(0, 0).get_data())) - 1))
    cerr << endl << "** Error in matrix cast **" << endl;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++)
      B(i, j) = A(i, j).get_si();


}

void matrix_cast(ZZ_mat<mpz_t>& B, ZZ_mat<long> A) {

  int n = A.get_rows();
  int d = A.get_cols();


  B.resize(n, d);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++)
      B(i, j) = A(i, j).get_data();

}

template<> void matrix_cast(ZZ_mat<__int128_t>& B, ZZ_mat<mpz_t> A) {

  int n = A.get_rows();
  int d = A.get_cols();


  B.resize(n, d);

  if (maxbitsize(A) > ((8 * sizeof(B(0, 0).get_data())) - 1))
    cerr << endl << "** Error in matrix cast **" << endl;

  Z_NR<__int128_t> tz;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++) {
      mpz_set_128int(tz, A(i, j));
      B(i, j) = tz;
    }
}

template<> void matrix_cast(matrix<Z_NR<__int128_t> >& B, matrix<Z_NR<mpz_t> > A) {

  int n = A.getRows();
  int d = A.getCols();


  B.resize(n, d);

  if (maxbitsize(A) > ((8 * sizeof(B(0, 0).get_data())) - 1))
    cerr << endl << "** Error in matrix cast **" << endl;

  Z_NR<__int128_t> tz;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++) {
      mpz_set_128int(tz, A(i, j));
      B(i, j) = tz;
    }
}



void matrix_cast(ZZ_mat<mpz_t>& B, ZZ_mat<__int128_t> A) {

  int n = A.get_rows();
  int d = A.get_cols();


  B.resize(n, d);

  Z_NR<mpz_t> tz;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++) {
      mpz_get_128int(tz, A(i, j));
      B(i, j) = tz;
    }

}


template<> void matrix_cast(matrix<Z_NR<long> >& B, matrix<Z_NR<mpz_t> >  A) {

  int n = A.getRows();
  int d = A.getCols();


  B.resize(n, d);

  if (maxbitsize(A) > ((8 * sizeof(B(0, 0).get_data())) - 1))
    cerr << endl << "** Error in matrix cast **" << endl;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++)
      B(i, j) = A(i, j).get_si();


}


inline int
utime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec * 1000000 + rus.ru_utime.tv_usec ;
}

inline int
utimesec ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec ;
}

// Discover a triangular structure, integer matrix
/******
    By JC on Dec. 7, 2018
    The input matrix B is with n rows and d columns.
    ------------------------------------------------------------------------------
    vargrind reports "Conditional jump or move depends on uninitialised value(s)"
    at: void hplll::matrix_structure<fplll::Z_NR<long> >
    (std::vector<int, std::allocator<int> >&, hplll::matrix<fplll::Z_NR<long> >,
    int, int)
    ------------------------------------------------------------------------------

*/
template<class T> inline void matrix_structure(vector<int>& structure, matrix<T> B, int n, int d)
{
  int i, k;
  structure.resize(d);



  // Bottom zeros
  for (k = 0; k < d; k++) {
    for (i = n - 1; (i >= 0) && (B(i, k) == 0); i--) { }
    structure[k] = i;
  }

  // Make it triangular for limiting changes during the computation
  for (k = 1; k < d; k++)   structure[k] = max(structure[k - 1], structure[k]);
}

// p two times greater than 2^bits

inline void next2prime(Z_NR<mpz_t>& p, long bits) {

  mpz_t mpzp;
  mpz_init(mpzp);
  mpz_set_si(mpzp, 1);
  mpz_mul_2exp(mpzp, mpzp, bits + 1);
  mpz_nextprime(p.get_data(), mpzp);

}

#ifdef NTL_GMP_LIP

// ZZ_mat to NTL Mat<zz_p>
// -----------------------


inline  void  zzmat_to_ntlp(mat_ZZ_p & Ap, const ZZ_mat<mpz_t> A, Z_NR<mpz_t> p) {

  mpz_t zp;
  mpz_init(zp);
  mpz_set(zp, p.get_data());

  // Mod p reduction and to NTL
  // --------------------------

  int n, m;
  n = A.get_rows();
  m = A.get_cols();

  int i, j;

  ZZ x;
  ZZ_p y;

  ZZ_p zero;
  zero = 0;

  mpz_t zx;
  mpz_init(zx);

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      if (A(i, j).sgn() == 0)
        Ap(i + 1, j + 1) = zero;
      else {
        mpz_mod(zx, A(i, j).get_data(), zp);
        if (mpz_sgn(zx) == 0)
          Ap(i + 1, j + 1) = zero;
        else {
          mpztozz(x, zx);
          conv(y, x);
          Ap(i + 1, j + 1) = y;
        }
      }
    }
}


// NTL Mat<zz_p> to ZZ_mat
// -----------------------


inline  void  ntlp_to_zzmat(ZZ_mat<mpz_t>& A, mat_ZZ_p Ap, Z_NR<mpz_t> p, long bits) {

  int n, m;
  n = A.getRows();
  m = A.getCols();


  int i, j;

  ZZ x;
  ZZ_p y;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      x = rep(Ap(i + 1, j + 1));
      if (IsZero(x))
        A(i, j) = 0;
      else
        zztompz(A(i, j).get_data(), x);
    }



  Z_NR<mpz_t> max;
  max = 1;
  max.mul_2si(max, bits);

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      if (A(i, j).cmp(max) >= 0) A(i, j).sub(A(i, j), p);

    }

}



// Matrix inversion through NTL
// ----------------------------

inline  void  NTL_inv(ZZ_mat<mpz_t>& V, const ZZ_mat<mpz_t> U) {

  int n = U.getRows();

  int i, j;

  mat_ZZ Untl;

  Untl.SetDims(n, n);

  int size = 0;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      size = max(size, size_in_bits(U(i, j)));

  char str[size + 20];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {

      mpz_get_str (str, 10, U(i, j).get_data());
      Untl(i + 1, j + 1) = to_ZZ(str);

    }

  mat_ZZ Vntl;
  Vntl.SetDims(n, n);

  ZZ det;

  int start = utime();
  inv(det, Vntl, Untl);
  start = utime() - start;
  cout << endl << "  Inv NTL time " << start / 1000 << " ms" << endl;

  size = 0;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      size = max(size, Vntl(i + 1, j + 1).size());

  ostringstream osmat;
  std::string k;
  char y[size * 32 + 20];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {

      osmat.str("");
      osmat << Vntl(i + 1, j + 1);
      k = osmat.str();
      strcpy(y, k.c_str());

      mpz_set_str((V(i, j).get_data()), y, 10);
    }


}

#endif // NTL

//

} // end namespace hplll

#endif
