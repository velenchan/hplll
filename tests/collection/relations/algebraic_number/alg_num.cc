/*
 * This code 
 * (1) generates a random univariate polynomial of degreen n,
 * (2) computes an approximation of a real root of the polynomial, and
 * (3) constructs the input form for hplll integer relation finding.
 *
 * Can be compiled with:
 *	g++ -Wall -O2 -o alg_num_gen alg_num.cc -lm -lgmp -lmps -lpthread
 *
 * Author: Jingwei Chen <jingwei.chen@outlook.com>
 *
 */

#include <iostream>
#include <mps/mps.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <string.h>

#define LOG_2_10 3.3219281

using namespace std;

int
main (int argc, char **argv)
{
  /* n is the degree of the polynomial,
   * i is used as counter */
  long int n = 0, i;
  long int prec, binary_prec, old_prec;
  mps_monomial_poly *p;
  mps_context *s;

  time_t t1;
  time(&t1);

  //unsigned long seed;


  // for gmp random integer
  gmp_randstate_t rstate;
  gmp_randinit_mt(rstate);



  mpq_t one, m_one, zero, t;
  mpz_t r, x, base, seed;

  mpz_init(x);
  mpz_init(r);
  mpz_init(base);
  mpz_init(seed);

  mpq_init (t);
  mpq_init (m_one);
  mpq_init(one);
  mpq_init (zero);

  mpq_set_si (m_one, -1, 1);
  mpq_set_si (zero, 0, 1);
  mpq_set_si (one, 1, 1);


  mpz_set_si(base, 2);

  /* Get n from command line */
  if (argc > 1)
    {
      n = atoi (argv[1]);
    }

  /* If parsing failed set n = 5 */
  if (n == 0)
    {
      n = 5;
    }

  /* 
  * 'prec' is to control the number of decimal digits of the root.
  * The 'binary_prec' is is the target precision of alpha^i.
  * For example, if 'binary_prec' is fixed, then 
  *	1. solve a root to decimal precision binary_prec/LOG_2_10 + 2*n;
  *	2. set 'mpf_prec' as binary_prec + 2*n;
  *	3. compute alpha^i for i= 2..n;
  *	4. output alpha^i with binary_prec/LOG_2_10 decimal precision. 
  */

  prec = (n+1)*(n+1);
  /* 
   * for the case of log H = n, binary_prec = 6.6*(n+1)^2
   * for the case of log H = n^1.5, binary_prec = 6.6*(n+1)^2.5
  */

  //prec *= (n+1); // remove this comment for the case of log H = n^1.5.
  prec *= 6.6;
  binary_prec = prec;
  prec /= LOG_2_10;

  //cout << "\n digital prec is "<<prec<<endl;

  s = mps_context_new ();
  p = mps_monomial_poly_new (s, n);

  mps_monomial_poly_set_coefficient_q (s, p, 0, m_one, zero);
    //f = mpq_out_str(stream, 10, m_one);
    //cout <<" "<<endl;

  //size_t f;
  //FILE *stream;


  unsigned long tmp;
  tmp = (unsigned long)t1;
  mpz_set_ui(seed, tmp);
  gmp_randseed(rstate, seed);

  for (i=1;i<=n;i++){

    mpz_rrandomb(x, rstate, 1+i);
    mpz_urandomm(r, rstate, base);
    if (mpz_cmp_si(r, 1) < 0) {
      mpz_mul_si(x, x, -1);
    }

    mpq_set_z(t, x);

    //f = mpq_out_str(stream, 10, t);
    //cout <<" "<<endl;
    mps_monomial_poly_set_coefficient_q (s, p, i, t, zero);
  }

  mpz_clear(r);
  //mpz_clear(x);
  //mpz_clear(base);
  mpz_clear(seed);
  gmp_randclear(rstate);
  // End of random number generation
  //cout<<"random coefficient finished\n"<<endl;




  mpz_pow_ui(x, base, n+1);
  mpq_set_z(t, x);
  //f = mpq_out_str(stream, 10, t);
  //cout<<"\n "<<endl;//cout<<n<<endl;
  //cout<<"coefficient generation finished"<<endl;
  mps_monomial_poly_set_coefficient_q (s, p, n, t, zero);
  // polynomial generation complete

  //cout<<"coefficient setting finished\n"<<endl;

  //mps_context_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA);
  mps_context_select_algorithm(s, MPS_ALGORITHM_STANDARD_MPSOLVE);


  /* Set the input polynomial */
  mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));

  /* Set the output precision*/
  mps_context_set_output_prec (s, prec+ 2*n + 4);

  /* Set the output goal*/
  mps_context_set_output_goal (s, MPS_OUTPUT_GOAL_APPROXIMATE);




  /* Set the outpout format*/
  mps_context_set_output_format (s, MPS_OUTPUT_FORMAT_FULL);


  old_prec = mpf_get_default_prec();
  mpf_set_default_prec(binary_prec+2*n + 4);
  //cout<<"ready to solve"<<endl;


  /* Allocate space to hold the results. We check only floating point results
   * in here */
  //mpc_t *results = mpc_valloc (n);
  mpc_t *results;



  /* Actually solve the polynomial */
  mps_mpsolve (s);
  //mps_msort (s);
  //

  //int degree;
  //degree = mps_context_get_degree (s);

  //cout<<"root's been found, degree is "<<degree<<endl;

  /* Save roots computed in the vector results */
  mps_context_get_roots_m (s, &results, NULL);


  //cout<<"result has been set"<<endl;

  mpf_t alpha, beta;
  mpf_init(alpha);
  mpf_init(beta);
  mpf_set_d(beta, 1e-52);

  /* Print out roots*/
  for (i = 0; i < n; i++)
    {
      //n_digits = mpc_out_str_2 (stdout, 10, n_digits, n_digits, results[i]);
      mpf_set (alpha, mpc_Im(results[i]));
      mpf_abs(alpha, alpha);
      if (mpf_cmp(beta, alpha)>0) {
        mpf_set (alpha, mpc_Re(results[i]));
        break;
      }
    }


  mpf_set_default_prec(old_prec);
  // print the sequence for  algebraic number reconstruction


  mpf_set_si (beta, 1);
  cout<<binary_prec<<endl;
  cout<< n + 1<<endl;
  cout<< 1 <<endl;

  for (i=1; i< n+1; i++){
    mpf_mul(beta, beta, alpha);
    mpf_out_str(stdout, 10, prec, beta);
    cout<<" "<<endl;
  }

  mpf_clear(alpha);
  mpf_clear(beta);

  //cout<<"YES"<<endl;

  //mpz_clear(r);
  mpz_clear(x);
  mpz_clear(base);
  //mpz_clear(seed);
  //gmp_randclear(state);

  free (results);
  //free (stream);
  free(s);
  free (p);

  return EXIT_SUCCESS;
}
