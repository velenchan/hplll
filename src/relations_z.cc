/* LLL based relation algorithms

Created Jeu  7 mar 2013 15:01:41 CET
        Ven 27 mar 2015 14:18:38 CET

Copyright (C) 2013      Gilles Villard
Modified Mar 12 sep 2017 15:34:12 CEST

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



#ifndef HPLLL_RELATIONS_Z_CC
#define HPLLL_RELATIONS_Z_CC



namespace hplll {


/************************************************************************************

   Temporary wrapping the call to fplll for no compilation error

   e.g. until  Z_NR<__int128_t> available

**************************************************************************************/



template<> int
FPTuple<long, double, matrix<FP_NR<double> > >::call_fplll(int &n_swaps, ZZ_mat<long> &b, ZZ_mat<long> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;
  ZZ_mat<long> u_invZ;


  status = lll_reduction(n_swaps, b, u, delta, eta, method, FT_DOUBLE, precision, 0);


  cout << "status 1 = "<< status << endl;

  return status;
}


template<> int
FPTuple<long, long double, matrix<FP_NR<long double> > >::call_fplll(int &n_swaps, ZZ_mat<long> &b, ZZ_mat<long> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(n_swaps, b, u, delta, eta, method, FT_LONG_DOUBLE, precision, 0);

    cout << "status 2 = "<< status << endl;

  return status;


}


template<> int
FPTuple<long, dd_real, matrix<FP_NR<dd_real> > >::call_fplll(int &n_swaps, ZZ_mat<long> &b, ZZ_mat<long> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(n_swaps, b, u, delta, eta, method, FT_DD, precision, 0);

  cout << "status 3 = "<< status << endl;

  return status;


}

template<> int
FPTuple<__int128_t, double, matrix<FP_NR<double> > >::call_fplll(int &n_swaps, ZZ_mat<__int128_t> &b, ZZ_mat<__int128_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  cerr << endl << "** Error in relations: no fplll with Z_NR<__int128_t> **" << endl;

  exit(EXIT_FAILURE);

  return 0;
}

template<> int
FPTuple<__int128_t, long double, matrix<FP_NR<long double> > >::call_fplll(int &n_swaps, ZZ_mat<__int128_t> &b, ZZ_mat<__int128_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  cerr << endl << "** Error in relations: no fplll with Z_NR<__int128_t> **" << endl;

  exit(EXIT_FAILURE);

  return 0;
}

template<> int
FPTuple<__int128_t, dd_real, matrix<FP_NR<dd_real> > >::call_fplll(int &n_swaps, ZZ_mat<__int128_t> &b, ZZ_mat<__int128_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  cerr << endl << "** Error in relations: no fplll with Z_NR<__int128_t> **" << endl;

  exit(EXIT_FAILURE);

  return 0;
}

template<> int
FPTuple<mpz_t, double, matrix<FP_NR<double> > >::call_fplll(int &n_swaps, ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(n_swaps, b, u, delta, eta, method, FT_DOUBLE, precision, 0);

  cout << "status 4= "<< status << endl;

  return status;

}

template<> int
FPTuple<mpz_t, long double, matrix<FP_NR<long double> > >::call_fplll(int &n_swaps, ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(n_swaps, b, u, delta, eta, method, FT_LONG_DOUBLE, precision, 0);

  cout << "status 5 = "<< status << endl;

  return status;

}



template<> int
FPTuple<mpz_t, dpe_t, MatrixPE<double, dpe_t> >::call_fplll(int &n_swaps, ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;


  status = lll_reduction(n_swaps, b, u, delta, eta, method, FT_DOUBLE, precision, 0);

  cout << "status 6 = "<< status << endl;

  return status;

}


template<> int
FPTuple<mpz_t, ldpe_t, MatrixPE<long double, ldpe_t> >::call_fplll(int &n_swaps, ZZ_mat<mpz_t> &b, ZZ_mat<mpz_t> &u, double delta, double eta,  \
    LLLMethod method, FloatType floatType,               \
    int precision, int flags) {

  int status;

  status = lll_reduction(n_swaps,b, u, delta, eta, method, FT_LONG_DOUBLE, precision, 0);

  cout << "status 7 = "<< status << endl;

  return status;

}




/***********************************************************************************

   Construction

   Fix the precision inside or let outside ?

**************************************************************************************/

// template<class ZT, class FT>
// FPTuple<ZT, FT>::FPTuple(vector<FP_NR<mpfr_t> > fpvin) {

//  d = fpvin.size();

//  fpv.resize(d);

//  for (int i = 0; i < d; i++)
//    fpv[i] = fpvin[i];


// }



template<class ZT, class FT, class MatrixFT>
FPTuple<ZT, FT, MatrixFT>::FPTuple(vector<FP_NR<mpfr_t> > fpvin) {

  d = fpvin.size();

  fpv.resize(d);

  for (int i = 0; i < d; i++)
    fpv[i] = fpvin[i];

  long emax, emin;

  emin = fpv[0].exponent();
  emax = fpv[0].exponent();


  for (int i = 0; i < d; i++) {

    emin = min(emin, fpv[i].exponent());
    emax = max(emax, fpv[i].exponent());

  }

  inputgap = emax - emin;

#ifdef _OPENMP
  omp_set_num_threads(S);
#endif

}


#ifdef _OPENMP
template<class ZT, class FT, class MatrixFT> int
FPTuple<ZT, FT, MatrixFT>::set_num_threads(int nbt) {

  S = nbt;

  omp_set_num_threads(S);

  return S;

}
#endif




/***********************************************************************************

    Relation_lll

    Calls LLL with successive lifts on mpz_t and FT bases

    alpha correct bits

 TODO: tune the parameters according to input data types

  truncate = -1 : no truncation
           = 0 : automatic choice
           > 0 : this number of bits

**************************************************************************************/


template<class ZT, class FT, class MatrixFT>  int
FPTuple<ZT, FT, MatrixFT>::relation(ZZ_mat<mpz_t>& C,  long alpha,
                                    long confidence_gap,  long shift, int truncate, int lllmethod,  int sizemethod, double delta) {


  ZZ_mat<mpz_t> L;
  L.resize(1, d);

  FP_NR<mpfr_t> t;

  for (int j = 0; j < d; j++) {
    t.mul_2si( fpv[j], alpha);
    L(0, j).set_f(t);
  }

  int found;



  // Shift should be maximized for efficiency, then the overall precision is at least truncate - shift  (truncate > shift)
  //
  // Hence enough precision should be ensured by truncate - shift
  //
  // However, truncate cannot be too large : the ratio between column norms after the shift will
  //    give extra size  : truncate + extra  <= integer size


  found = relation_lll(C, L, alpha, confidence_gap, shift, truncate, lllmethod, sizemethod, delta);

  return found;

}




/***********************************************************************************

    Companion to relation_lll
    Relation from an integer matrix

    Calls LLL with successive lifts on mpz_t and FT bases

    alpha correct bits

**************************************************************************************/

template<class ZT, class FT, class MatrixFT>  int
FPTuple<ZT, FT, MatrixFT>::relation_lll(ZZ_mat<mpz_t>& C, ZZ_mat<mpz_t> A, long alpha,
                                        long confidence_gap, long shift, int truncate, int lllmethod,  int sizemethod, double delta) {


    //cout << "\n begin relation_lll().\n"<<endl;

  Timer time;

  Timer tlll;
  Timer tprod;
  Timer ttrunc;


    //unsigned int swapnum = 0;
    //unsigned int testnum = 0;


#ifdef _OPENMP
  double en, st;
  double lllt = 0.0;
  double prodt = 0.0;
  double trunct = 0.0;
#endif

  int m, d;
  int i, j;

  m = A.get_rows();
  d = A.get_cols();

  ZZ_mat<mpz_t> A_in;
  A_in.resize(m + d, d);


  // **** m=1 for the moment

  for (j = 0; j < d; j++)
    A_in(0, j) = A(0, j);

  for (i = 0; i < d; i++)
    A_in(m + i, i) = 1;


  int bitsize = maxbitsize(A, 0, m, d);

  // For assigning the truncated basis at each step



  int def = -bitsize;

  int target_def = -bitsize + alpha;

  int found = 0;

  FP_NR<mpfr_t> quot, new_quot, tf;
  new_quot = 1.0;  // new_quot should be bigger after the first iteration of the loop
  new_quot.mul_2si(new_quot, alpha);

  FP_NR<mpfr_t> gap;
  gap = 1.0;

  FP_NR<mpfr_t> confidence;
  // For testing 1/gap < confidence
  confidence = 1.0;
  // relié, plus petit,  au shift sur S (ex 80)
  confidence.mul_2si(confidence, -confidence_gap - shift); // La baisse absolue est plus petite que le shift

  //cout << "confidence = " << confidence << endl;
  //cout << "confidence.exponent = " << confidence.exponent() << endl;

  FP_NR<mpfr_t> epsilon;
  epsilon = 10.0;


  Z_NR<mpz_t> tz, t, maxcol;

  long foundcol;

  ZZ_mat<ZT> T, TT;

  T.resize(m + d, d);
  TT.resize(d, m + d);


  ZZ_mat<ZT> U, UT;

  U.resize(d, d);
  UT.resize(d, d);




  Lattice<ZT, FT, matrix<Z_NR<ZT> >,  MatrixFT > Bp(T, TRANSFORM, sizemethod);

    int n_swaps = 0;
    int nb_swaps = 0;

    //cout << " position 1. \n" << endl;



/*
    ZZ_mat<ZT> VT;
    VT.resize(d, d);
    int gso_flags = 2; // This option means LM_FAST is used.
    int flags = 0;

    cout << " position 2. \n" << endl;



    MatGSO<Z_NR<ZT>, FP_NR<FT>> m_gso(T, UT, VT, gso_flags);
    LLLReduction<Z_NR<ZT>, FP_NR<FT>> lll_obj(m_gso, LLL_DEF_DELTA, LLL_DEF_ETA, flags);

    cout << " position 3. \n" << endl;
*/


  // Main loop on the shifts
  // -----------------------

  while (def < target_def) {



    //HPLLL_INFO("Current default: ", def);

    if ((target_def - def) <= shift)
      def = target_def;
    else def += shift;


    //HPLLL_INFO("now default: ", def);
#ifdef _OPENMP
    st = omp_get_wtime();
#else
    time.start();
#endif


    if (truncate == -1)
      lift_truncate(T, A_in, def, 0);
    else if (truncate == 0)
      lift_truncate(T, A_in, def, shift + 2 * d);
    else
      lift_truncate(T, A_in, def, truncate);

#ifdef _OPENMP
    en = omp_get_wtime();
    trunct += (en - st);
#else
    time.stop();
    ttrunc += time;
#endif


    // HLLL and DOUBLES

    if (lllmethod == HLLL) {


#ifdef _OPENMP
      st = omp_get_wtime();
#else
      time.start();
#endif

      Bp.assign(T);

      Bp.hlll(delta);




#ifdef _OPENMP
      en = omp_get_wtime();
      lllt += (en - st);
#else
      time.stop();
      tlll += time;
#endif



#ifdef _OPENMP
      st = omp_get_wtime();

      if (S > 1)
        pmatprod_in_int(A_in, Bp.getU(), S);
      else
        matprod_in_int(A_in, Bp.getU());

      en = omp_get_wtime();
      prodt += (en - st);
#else
      time.start();

      matprod_in_int(A_in, Bp.getU());

      time.stop();
      tprod += time;
#endif


      //cout << "sizeof U: " << maxbitsize(Bp.getU(), 0, d, d) << endl;
      //cout << "sizeof basis: " << maxbitsize(A_in, 1, d + 1, d) << endl << endl;



    } // end HLLL


    // FPLLL

    else if (lllmethod == FPLLL) {
    cout << " FPLLL starts \n" << endl;


#ifdef _OPENMP
      st = omp_get_wtime();
#else
      time.start();
#endif

      transpose(TT, T);

      setId(UT);

      /* The following is only for rs_5_5 to rs_10_10. The main goal
        is to count the number of LLL-swaps after using lift-truncate.
       */
        call_fplll(nb_swaps, TT, UT, delta, 0.51, LM_FAST, FT_DOUBLE, 0);
        //lll_reduction(nb_swaps, TT, UT, delta, 0.51, LM_FAST, FT_DOUBLE, 0, 0);
        n_swaps += nb_swaps;



#ifdef _OPENMP
      en = omp_get_wtime();
      lllt += (en - st);
#else
      time.stop();
      tlll += time;
#endif




#ifdef _OPENMP
      st = omp_get_wtime();

      transpose(U, UT);

      if (S > 1)
        pmatprod_in_int(A_in, U, S);
      else
        matprod_in_int(A_in, U);

      en = omp_get_wtime();
      prodt += (en - st);
#else
      time.start();

      transpose(U, UT);



      matprod_in_int(A_in, U);

      time.stop();
      tprod += time;
#endif



      //cout << "sizeof U: " << maxbitsize(Bp.getU(), 0, d, d) << endl;
      //cout << "sizeof A_in: " << maxbitsize(A_in, 0, d + 1, d) << endl << endl;
      //cout << "sizeof the unimodular basis: " << maxbitsize(A_in, 1, d + 1, d) << endl << endl;


    } // end FPLLL

    if (lllmethod == HLLL) {
        n_swaps =  Bp.nbswaps;
    }

    // Test
    // ----

    /* GV Mer 20 sep 2017 12:57:27 CEST
        Mofified, the small entry may appear in any column, if not in the first one and
         if integer in fixed precision (ex long) then the lift truncate may lead to an error
         for the detection (bad truncation and crash)
    */


/*
 * For tracing rs_12_13, with 20, 20, 40 parameter.
 *
 * ZZ_mat<mpz_t> tmp;//16889777313727009140128160;

	tmp.resize(1,3);
	tmp(0, 0) = 1;
	tmp(0, 1) = 104099653179802559;
	tmp(0, 2) = 162246240;
	tmp(0, 0).mul(tmp(0, 1), tmp(0, 2));

	cout << "the product is " << tmp(0, 0) << endl;

	Z_NR<mpz_t> ttt;

    for (int i = m; i < d+m; i++){
    	for (int kk = 0; kk < d; kk++){
		ttt.abs(A_in(i, kk));
		if (tmp(0,0).cmp(ttt) == 0){
			cout << "YES, FOUND HERE!" << endl;
		}
	}
    }
    cout << "NOT FOUND YET. " << endl;
*/


    quot = new_quot;
    //cout << " quot = " << quot << endl;

    foundcol = 0;

    if (A_in(0, 0).sgn() == 0) { // For making the gap pertinent even if 0
      tz = 1;
    }
    else
      tz.abs(A_in(0, 0));


    for (int i = 1; i < d; i++) {

      t.abs(A_in(0, i));

      if (t.cmp(tz) == -1) {
        tz = t;
        foundcol = i;
      }

    }// now, tz = min |A_in(0, i)| = |A_in(0, foundcol)|.

    new_quot.set_z(tz);

    //cout << "new_quot = tz = " << new_quot << endl;


    maxcol.abs(A_in(0, foundcol));
    maxcol.mul_2si(maxcol, def);


    //cout << " initial maxcol = " << maxcol << endl;

    for (i = 0; i < d; i++) {
      tz.abs(A_in(m + i, foundcol));
      if (tz.cmp(maxcol) == 1) maxcol = tz;
    }

    tf.set_z(maxcol);

    //cout << "tf = maxcol = " << tf << endl;

    new_quot.div(new_quot, tf);

    //cout << "new_quot = tz/tf = " << new_quot << endl;

    gap.div(new_quot, quot);
    gap.abs(gap);

    //cout << "gap = |new_quot/quot| = " << gap;


    /* GV Mer 20 sep 2017 12:57:27 CEST
        added the test using the input gap otherwise the computation may terminate
        without having shifted for discovering all non zero entries
    */

    //cout << " expondent of gap = " << gap.exponent() << endl;


    if (def > inputgap - bitsize) {

      //if ((gap.cmp(confidence) == -1) && (new_quot.cmp(epsilon) == -1)) {
      // Si epsilon mettre à la valeur max quotient des nombres au départ

      if (gap.cmp(confidence) == -1) {
        C.resize(d, 1);
        for (j = 0; j < d; j++)
          C(j, 0) = A_in(m + j, foundcol);

        gap.mul_2si(gap, shift);

        cout << endl;
        HPLLL_INFO(def + bitsize, " bits used");
        HPLLL_INFO("Candidate relation found with bit confidence: ", -gap.exponent());

        //cout << endl << "Time lll: " << tlll << endl;
        //cout << "Time products: " << tprod << endl << endl;

#ifdef _OPENMP
        cout << endl << "Time trunc: " << trunct << endl;
        cout << endl << "Time lll: " << lllt << endl;
        cout << "Time products: " << prodt << endl << endl;
#else
        cout << endl << "Time trunc: " << ttrunc << endl;
        cout << endl << "Time lll: " << tlll << endl;
        cout << "Time products: " << tprod << endl << endl;
#endif



        cout << endl << "The total number of swaps is: " << n_swaps <<endl;

        return 1;

      }

    }




  } // End while on the shift

  HPLLL_INFO(alpha, " digits used");

// Relation bound
// --------------

  unsigned oldprec;
  oldprec = mpfr_get_default_prec();

  mpfr_set_default_prec(2 * d);

  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > B(A_in, NO_TRANSFORM, DEF_REDUCTION);

  B.householder();

  matrix<FP_NR<mpfr_t> > R;

  R = B.getR();

  FP_NR<mpfr_t> minr, tr;

  minr.abs(R(0, 0));
  for (i = 1; i < d; i++) {
    tr.abs(R(i, i));
    if (minr.cmp(tr) > 0) minr = tr;
  }

  cerr << endl << "** No relation found with min Rii = " << minr << endl;

  mpfr_set_default_prec(oldprec);


  return found; // 0 here


};



/***********************************************************************************

    Direct LLL reduction

    Note from J. Chen: ZT would be better with Z_NR<mpz_t>, or else the input data will be wrong,
          due to the large scale parameter (2^alpha) of the integer relation lattice.

**************************************************************************************/


template<class ZT, class FT, class MatrixFT>
int FPTuple<ZT, FT, MatrixFT>::lll(ZZ_mat<mpz_t>& C,  long alpha,
                                    int lllmethod,  int sizemethod, double delta) {



  int i, j;

  int m = 1;

  ZZ_mat<ZT> L;
  L.resize(1, d);

  FP_NR<mpfr_t> t;

  for (int j = 0; j < d; j++) {
    t.mul_2si( fpv[j], alpha);
    L(0, j).set_f(t);
  }

    ZZ_mat<ZT> U;
    U.resize(d, d);



  int found = 0;

  ZZ_mat<ZT> A;
  A.resize(m + d, d);

  // **** m=1 for the moment
  for (j = 0; j < d; j++)
    A(0, j) = L(0, j);

  for (i = 0; i < d; i++){
    A(m + i, i) = 1;
  }

    //cout<< " before LLL reduction, the input basis is " << endl << A <<endl;


  if (lllmethod == HLLL) {



    Lattice<ZT, FT, matrix<Z_NR<ZT> >,  MatrixFT > B(A);

    cout << " begin direct HLLL" << endl;

    B.hlll(delta);

    cout << endl << "The total number of swaps is: " << B.nbswaps <<endl;

    A = B.getbase();





    /*
        Added by J. Chen on Nov. 29, 2018

        If we use the following specification:
            FPTuple<mpz_t, double, matrix<FP_NR<double> > > L(fpv);
            L.lll(C, alpha, HLLL, DEF_REDUCTION, 0.99);
        then when alpha > 511, say alpha = 512,
            ./relation < ../data/rs_5_5
        will run with #swaps up to 4294967295, which is
            nblov_max,
        set in the constructor of the class Lattice in hlll.cc.



    */
  }


  else if (lllmethod == FPLLL) {

    ZZ_mat<ZT> T;
    int n_swaps;
    T.resize(d, d + m);

    transpose(T, A);

    cout << " begin direct FPLLL" << endl;

    // Put the wrapper also if large examples
    lll_reduction(n_swaps, T, delta, 0.51, LM_FAST, FT_DEFAULT, 0, 0);

    transpose(A, T);

    cout << "the number of swaps is " << n_swaps << endl <<  endl;


  } // end FPLLL


  int foundcol = 0;

    for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
            U(i, j) = A(i + 1, j);



  matprod_in(L, U);


  Z_NR<ZT> tz, tt, maxcol;

  tz.abs(L(0, 0));

  for (i = 1; i < d; i++) {

    tt.abs(L(0, i));

    if (tt.cmp(tz) == -1) {
      tz = tt;
      foundcol = i;
    }

  }




  C.resize(d, 1);
  for (j = 0; j < d; j++)
    C(j, 0) = A(m + j, foundcol);



  // cout << endl;
  // HPLLL_INFO(def + bitsize, " bits used");
  // HPLLL_INFO("Candidate relation found with bit confidence: ", -gap.exponent());

  // cout << endl << "Time lll: " << tlll << endl;
  // cout << "Time products: " << tprod << endl << endl;
  // return 1;


  // TODO: get status of LLL and assign found accordingly
  return found; // 0 here


};




} // end namespace hplll


#endif











