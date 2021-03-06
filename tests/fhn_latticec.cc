/* 
Ven  3 jui 2016 15:04:49 CEST
Copyright (C) 2016      Gilles Villard 

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


#include "hplll.h"


#include <NTL/LLL.h>

using namespace NTL;

/* ***********************************************

          MAIN   

   ********************************************** */

using namespace hplll; 

int main(int argc, char *argv[])  {

  char results[]="benchmarks_results/fhn_latticec.results";    // ******** SPECIALIZE
  
  filebuf fb;
  iostream os(&fb);
  fb.open (results,ios::out);           


  
  ZZ_mat<mpz_t> A; // For hpLLL 
  ZZ_mat<mpz_t> AT,tmpmat;  // fpLLL  

  // ---------------------------------------------------------------------

  int k,K;

  vector<char*> s(100);
 
  k=0;

  //------------

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-200",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-450",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-650",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-800",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-900",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1000",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1100",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1200",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1300",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1400",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1500",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1600",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1700",255);
  k+=1;

  s[k]=(char *) malloc(256);
  strncpy(s[k], "collection/latticec/challenge-1800",255);
  k+=1;

    
  cout << s[0] << endl; 
    
  //-------------

  K=k;

  double delta=0.99;

  int run=0;
  
  Timer time;

  int status=0;

    os << endl << "FPLLL, HPLLL, NTL running times / lattice challenge bases / long" << endl;   // ******** SPECIALIZE
    os << endl << "FPLLL stopped at dim 800, otherwise wrapper with mpfr too costly" << endl;
    os << endl << "To see: to tune with hlll wrapper if cond too big" << endl; 
                                                                                    
    os <<         "----------------------------------------------------------------" << endl << endl;
 
    for (int k=0; k<K; k++) { 


      /*****************************************************************************/
      /*   i-th run  */
      /*****************************************************************************/
      
      run+=1;

      // Input basis 
      fb.close();
      fb.open (s[k],ios::in);
      os >>  AT ;
      fb.close();
      fb.open (results,ios::app);    

      int d=AT.get_rows();
      int n=AT.get_cols();

      A.resize(n,d);
      
      transpose(A,AT);

      ZZ_mat<long> Along;
      matrix_cast(Along,A);

      cout << n << "   " << d << endl; 
      
      cout << "--------------  HLLL" << endl << endl; 

      {

      	os << endl << "------------------------------------------------ " << endl ;

	Lattice<long, double, matrix<Z_NR<long> >,  matrix<FP_NR<double> > > B(Along,NO_TRANSFORM,DEF_REDUCTION);  //* name 

	verboseDepth=1;
      	time.start();
      	status=B.hlll(delta); //* name
      	time.stop();
	verboseDepth=0;
 
      	os << "Run " << run << "  with n,d = " << n << "  " << d << ",    delta = " << delta <<  endl << endl;
      	os << "    hlll: " << time << endl ;
      	time.print(os);
      	os << endl;

	matrix_cast(A,B.getbase());
	
      	if (status ==0) {
      	  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A,NO_TRANSFORM,DEF_REDUCTION); //* names

      	  T.isreduced(delta-0.1); //* name

	  double t,u,v,w;

	  hplll::ratio<mpz_t>(A,t,u,v,w);

	  cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
	  cout << ".. Average diagonal ratio: " << u << endl;
	  cout << ".. Max diagonal ratio: " << v << endl;
	  cout << ".. First vector quality: " << w << endl;
      	} 
      	cout << endl; 

      	cout << "--------------  FPLLL WRAPPER VERBOSE " << endl << endl; 

	if (n <=800) {
	  time.start();
	  lll_reduction(AT, delta, 0.501, LM_WRAPPER,FT_DEFAULT,0,LLL_VERBOSE);
	  time.stop();
	  

	  os << "   fplll: " << time << endl << endl ;
	  time.print(os);
	  os << endl;
	  
	  transpose(A,AT);
	  Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(A,NO_TRANSFORM,DEF_REDUCTION); //* name
	  T2.isreduced(delta-0.1); //* name
	}
	cout << "--------------  NTL  " << endl << endl; 

	Mat<ZZ> BN; 

	// Input basis 
	fb.close();
	fb.open (s[k],ios::in);
	os >>  BN ;
	fb.close();
	fb.open (results,ios::app);

	time.start();	
       
	LLL_FP(BN,0.99,0,0,1); 
	
	time.stop();

	fb.close();
	fb.open ("tmp.txt",ios::out);
	os <<  BN ;
	fb.close();
	fb.open ("tmp.txt",ios::in);
	os >> AT;
	fb.close();
	system("rm tmp.txt");
	fb.open (results,ios::app);

	transpose(A,AT);

	os << "   ntl: " << time << endl << endl ;
      	time.print(os);
      	os << endl;
	
	Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(A,NO_TRANSFORM,DEF_REDUCTION);
	verboseDepth=1;
	T.isreduced(delta-0.1);
	verboseDepth=0;
	
	os << "   ntl after check: " << time << endl << endl ;
      	time.print(os);
      	os << endl;
	
      } 
   
    }// End on runs, k loop


    // END 
    fb.close();


 
  return 0;
}
