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

    char results[] = "benchmarks_results/fhn_table_knapsack.results";  // ******** SPECIALIZE

    filebuf fb;
    iostream os(&fb);
    fb.open (results, ios::out);



    ZZ_mat<mpz_t> A; // For hpLLL


    ZZ_mat<mpz_t> AT, tmpmat, TAT; // fpLLL

// ---------------------------------------------------------------------

    int k, K;


    k = 0;

    vector<int> d(100);

    k = 0;

//------------

    d[k] = 20;
    k += 1;

    d[k] = 200;
    k += 1;

    d[k] = 240;
    k += 1;

    d[k] = 280;
    k += 1;

    d[k] = 300;
    k += 1;

    d[k] = 320;
    k += 1;

    d[k] = 340;
    k += 1;

//-------------

    K = k;

    double delta = 0.99;

    int run = 0;

    Timer time;

    int status;

    os << endl << "FPLLL, HPLLL, NTL running times / Knapsack bases " << endl;   // ******** SPECIALIZE

    os <<         "-------------------------------------------------" << endl << endl;

    for (int k = 0; k < K; k++) {



        /*****************************************************************************/
        /*   i-th run  */
        /*****************************************************************************/

        run += 1;



        A.resize(d[k] + 1, d[k]);
        AT.resize(d[k], d[k] + 1);

        tmpmat.resize(d[k] + 1, d[k]);

        AT.gen_intrel(d[k] * 100);
        transpose(A, AT);


        cout << d[k] <<  endl;

        cout << "--------------  HLLL" << endl << endl;

        os << endl << "------------------------------------------------ " << endl ;


        // Temporairement hlll comme dans un wrapper avec fplll et pas isolément


        //Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  B(A, NO_TRANSFORM, DEF_REDUCTION); //* name

        //time.start();
        //status = B.hlll(delta); //* name
        //time.stop();

        TAT.resize(d[k], d[k] + 1);
        transpose(TAT, A);

        time.start();

        status = lll_reduction(TAT, delta, 0.501, LM_FAST, FT_DOUBLE, 0, LLL_VERBOSE);

        transpose(A, TAT);

        Lattice<mpz_t, dpe_t, matrix<Z_NR<mpz_t> >, MatrixPE<double, dpe_t> >  B(A, NO_TRANSFORM, DEF_REDUCTION); //* name

        if (status != 0) {

            status = B.hlll(delta); //* name

        }


        time.stop();


        os << "Run " << run << "  with dim = " << d[k] << ",  bits = " << d[k] * 100 << ",   delta = " << delta <<  endl << endl;
        os << "    hlll: " << time << endl ;
        time.print(os);
        os << endl;


        if (status == 0) {

            Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(B.getbase(), NO_TRANSFORM, DEF_REDUCTION); //* names

            T.isreduced(delta - 0.1); //* name

            double t, u, v, w;

            hplll::ratio<mpz_t>(B.getbase(), t, u, v, w);

            cout << endl << ".. log 2 Frobenius norm cond: " << t << endl;
            cout << ".. Average diagonal ratio: " << u << endl;
            cout << ".. Max diagonal ratio: " << v << endl;
            cout << ".. First vector quality: " << w << endl;

        }

        cout << endl;

        cout << "--------------  FPLLL WRAPPER VERBOSE " << endl << endl;

        time.start();
        lll_reduction(AT, delta, 0.501, LM_WRAPPER, FT_DEFAULT, 0, LLL_VERBOSE);
        time.stop();


        os << "   fplll: " << time << endl << endl ;
        time.print(os);
        os << endl;


        transpose(tmpmat, AT);
        Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T2(tmpmat, NO_TRANSFORM, DEF_REDUCTION); //* name
        T2.isreduced(delta - 0.1); //* name

        cout << "--------------  NTL  " << endl << endl;

        // Mat<ZZ> BN;

        // // Input basis

        // fb.close();
        // fb.open ("tmp.txt", ios::out);
        // os <<  transpose(A) ;
        // fb.close();
        // fb.open ("tmp.txt", ios::in);
        // os >> BN;
        // fb.close();
        // system("rm tmp.txt");
        // fb.open (results, ios::app);


        // time.start();

        // LLL_XD(BN, 0.99, 0, 0, 1);

        // time.stop();

        // fb.close();
        // fb.open ("tmp.txt", ios::out);
        // os <<  BN ;
        // fb.close();
        // fb.open ("tmp.txt", ios::in);
        // os >> AT;
        // fb.close();
        // system("rm tmp.txt");
        // fb.open (results, ios::app);



        // transpose(tmpmat, AT);

        // Lattice<mpz_t, mpfr_t, matrix<Z_NR<mpz_t> >, matrix<FP_NR<mpfr_t> > > T(tmpmat, NO_TRANSFORM, DEF_REDUCTION);
        // verboseDepth = 0;
        // T.isreduced(delta - 0.1);


        // os << "   ntl: " << time << endl << endl ;
        // time.print(os);
        // os << endl;


    }// End on runs, k loop


// END
    fb.close();



    return 0;
}
