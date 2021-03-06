/* integer relations test file

Copyright (C) 2013      Gilles Villard

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

#include "matgen.h"
#include "relations.h"

using namespace hplll;

/* ***********************************************

          MAIN

   ********************************************** */


int main(int argc, char *argv[])  {

  //int status = 0;

  int n;
  int r, s;

  int found;

  int difference;

  int succeed = 0;
  int nbtest = 0;

  filebuf fb;
  iostream os(&fb);

  int nbrel;

  long setprec;

  //  *****************************************************
  cout <<  "Testing relation finder" << endl;
  //  *****************************************************



  typedef mpz_t integer_t;

  ZZ_mat<integer_t> C;  // Output relations
  ZZ_mat<integer_t> Ccheck;  // Output relations

  //  -------------------- TEST i --------------------------------
  nbtest += 1;

  { // for changing the mpfr prec in A

    vector<FP_NR<mpfr_t> > fpv;   // Input fp vector

    r = 4;
    s = 4;
    n = r * s + 1;


    setprec = 280;
    mpfr_set_default_prec(setprec);

    gen3r2s(fpv, n, r, s);


    nbrel = 1;

    cout << "     Relation test, dim = " << n << ", " << setprec << " bits " << endl;

    FPTuple<mpz_t, dpe_t, MatrixPE<double, dpe_t> > L(fpv);

    found = L.relation(C, setprec, 20, 20, -1, FPLLL);

    //found = relation<mpz_t, dpe_t, MatrixPE<double, dpe_t> >(C, A, setprec, 20, 20, -1, FPLLL);
    //found=relation_f<long, double>(C, A,setprec,20,20,10,FPLLL);

    cout << C << endl;

    Ccheck.resize(n, 1);
    fb.open ("C1_out", ios::in);
    os >> Ccheck ;
    fb.close();

    if (found != 1)
      cerr << "*** Problem in relation test, no relation found" << endl;

    if (nbrel == 1) {
      difference = !matcmp(C, Ccheck, 1, n);
      if (difference) {
        cerr << "*** Invalid matrix comparison in relation test" << endl;
      }
      else
        succeed += 1;
    }

  }

  //  -------------------- TEST i --------------------------------
  nbtest += 1;

  { // for changing the mpfr prec in A


    vector<FP_NR<mpfr_t> > fpv;   // Input fp vector

    r = 7;
    s = 7;
    n = r * s + 1;

    setprec = 1800;
    mpfr_set_default_prec(setprec);

    gen3r2s(fpv, n, r, s);


    nbrel = 1;

    cout << "     Relation test, dim = " << n << ", " << setprec << " bits " << endl;

    FPTuple<mpz_t, dpe_t, MatrixPE<double, dpe_t> > L(fpv);

    found = L.relation(C, setprec, 60, 200, -1, FPLLL);

    //found = relation<mpz_t, dpe_t, MatrixPE<double, dpe_t> >(C, A, setprec, 60, 200, -1, FPLLL);
    //found=relation_f<long, double>(C, A,setprec,20,20,10,FPLLL);

    cout << C << endl;

    Ccheck.resize(n, 1);
    fb.open ("C2_out", ios::in);
    os >> Ccheck ;
    fb.close();

    if (found != 1)
      cerr << "*** Problem in relation test, no relation found" << endl;

    if (nbrel == 1) {
      difference = !matcmp(C, Ccheck, 1, n);
      if (difference) {
        cerr << "*** Invalid matrix comparison in relation test" << endl;
      }
      else
        succeed += 1;
    }

  }

//  -------------------- TEST i --------------------------------
  nbtest += 1;

  { // for changing the mpfr prec in A


    vector<FP_NR<mpfr_t> > fpv;   // Input fp vector

    r = 7;
    s = 7;
    n = r * s + 1;

    setprec = 1800;
    mpfr_set_default_prec(setprec);

    gen3r2s(fpv, n, r, s);


    nbrel = 1;

    cout << "     Relation test, dim = " << n << ", " << setprec << " bits " << endl;

    FPTuple<long, double, matrix<FP_NR<double> > >  L(fpv);

    found = L.relation(C, setprec, 20, 20, 40, FPLLL);


    cout << C << endl;

    Ccheck.resize(n, 1);
    fb.open ("C2_out", ios::in);
    os >> Ccheck ;
    fb.close();

    if (found != 1)
      cerr << "*** Problem in relation test, no relation found" << endl;

    if (nbrel == 1) {
      difference = !matcmp(C, Ccheck, 1, n);
      if (difference) {
        cerr << "*** Invalid matrix comparison in relation test" << endl;
      }
      else
        succeed += 1;
    }

  }

  //  -------------------- TEST i --------------------------------
  nbtest += 1;

  { // for changing the mpfr prec in A


    vector<FP_NR<mpfr_t> > fpv;   // Input fp vector

    static string s;

    fb.open ("C3_in", ios::in);

    os >> setprec ;
    os >> n;

    mpfr_set_default_prec(setprec);

    fpv.resize(n);
    for (int i = 0; i < n; i++) {
      os >> s;
      mpfr_set_str (fpv[i].get_data(), s.c_str(), 10, GMP_RNDN);
    }

    fb.close();



    nbrel = 1;

    cout << "     Relation test, dim = " << n << ", " << setprec << " bits " << endl;

    FPTuple_f<long, double> L(fpv);

    found = L.relation(C, setprec, 80, 20, 10, HLLL);

    cout << C << endl;

    Ccheck.resize(n, 1);
    fb.open ("C3_out", ios::in);
    os >> Ccheck ;
    fb.close();

    if (found != 1)
      cerr << "*** Problem in relation test, no relation found" << endl;

    if (nbrel == 1) {
      difference = !matcmp(C, Ccheck, 1, n);
      if (difference) {
        cerr << "*** Invalid matrix comparison in relation test" << endl;
      }
      else
        succeed += 1;
    }

  }

//  -------------------- TEST i --------------------------------
  nbtest += 1;

  { // for changing the mpfr prec in A


    vector<FP_NR<mpfr_t> > fpv;   // Input fp vector

    static string s;

    fb.open ("C3_in", ios::in);

    os >> setprec ;
    os >> n;

    mpfr_set_default_prec(setprec);

    fpv.resize(n);
    for (int i = 0; i < n; i++) {
      os >> s;
      mpfr_set_str (fpv[i].get_data(), s.c_str(), 10, GMP_RNDN);
    }

    fb.close();



    nbrel = 1;

    cout << "     Relation test, dim = " << n << ", " << setprec << " bits " << endl;

    FPTuple<long, double, matrix<FP_NR<double> > > L(fpv);

    found = L.relation(C, setprec, 20, 20, 40, HLLL);

    cout << C << endl;

    Ccheck.resize(n, 1);
    fb.open ("C3_out", ios::in);
    os >> Ccheck ;
    fb.close();

    if (found != 1)
      cerr << "*** Problem in relation test, no relation found" << endl;

    if (nbrel == 1) {
      difference = !matcmp(C, Ccheck, 1, n);
      if (difference) {
        cerr << "*** Invalid matrix comparison in relation test" << endl;
      }
      else
        succeed += 1;
    }

  }

  //  -------------------- TEST i --------------------------------
  nbtest += 1;

  { // for changing the mpfr prec in A


    vector<FP_NR<mpfr_t> > fpv;   // Input fp vector

    r = 8;
    s = 8;
    n = r * s + 1;

    setprec = 2800;
    mpfr_set_default_prec(setprec);

    gen3r2s(fpv, n, r, s);


    nbrel = 1;

    cout << "     Relation test, dim = " << n << ", " << setprec << " bits " << endl;

    FPTuple_f<long, double> L(fpv);

    found = L.relation(C, setprec, 60, 200, 20, FPLLL);

    cout << C << endl;

    Ccheck.resize(n, 1);
    fb.open ("C4_out", ios::in);
    os >> Ccheck ;
    fb.close();

    if (found != 1)
      cerr << "*** Problem in relation test, no relation found" << endl;

    if (nbrel == 1) {
      difference = !matcmp(C, Ccheck, 1, n);
      if (difference) {
        cerr << "*** Invalid matrix comparison in relation test" << endl;
      }
      else
        succeed += 1;
    }

  }

  


  //  *****************************************************
  cout << endl << "     " << succeed << " relations tests ok over " << nbtest << endl;
  //  *****************************************************

  if (succeed == nbtest)
    return 0;
  else
    return -1;

  //return status;
}





