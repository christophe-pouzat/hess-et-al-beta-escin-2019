#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

/** @brief Returns the Anderson-Darling statistics assuming a
 *  standard normal distribution for the data
 *
 *  The data are contained in the `gsl_vector` pointed to
 *  by `data`. If the content is not sorted (`sorted==false`)
 *  the data are first copied before being sorted.
 *
 *  @param[in] data pointer to a `gsl_vector` containing the data
 *  @param[in] sorted a boolean indicated if the `data` content is
 *             already sorted (`true`) or not (`false`)
 *  @returns a double with the Anderson-Darling statistics
*/
double AndersonDarling_W2(gsl_vector * data, bool sorted)
{
  gsl_vector * data_s;
  if (sorted == false)
  {
    data_s = gsl_vector_alloc(data->size);
    gsl_vector_memcpy(data_s,data);
    gsl_sort_vector(data_s);
  }
  else
  {
    data_s = data;
  }
  size_t n = data->size;
  double n_d = (double) n;
  double A=0.;
  for (size_t i=1; i<=n; i++)
  {
    double y = gsl_vector_get(data_s,i-1);
    double Phi_at_y = gsl_cdf_gaussian_P(y,1.0);
    A += (2*i-1.)*log(Phi_at_y)+(2*(n_d-i)+1)*log(1-Phi_at_y);
  }
  A /= -n_d;
  if (sorted == false)
    gsl_vector_free(data_s);
  return A-n_d;
}

/** @brief Returns the asymptotic cdf of the Anderson-Darling
    statistics.
 *
 *  Adaptation of function `adinf` of Marsaglia & Marsaglia (2004)
 *  [J. Stat. Software 9(2): 1-5](https://www.jstatsoft.org/article/view/v009i02).
 *
 *  @param[in] z a double the observed statistics value
 *  @returns a double Prob{W2 <= z}
*/
double adinf(double z)
{
  if(z<2.)
    return exp(-1.2337141/z)/sqrt(z)*(2.00012+(.247105-(.0649821-(.0347962-(.011672-.00168691*z)*z)*z)*z)*z);
  /* max |error| < .000002 for z<2, (p=.90816...) */
 return exp(-exp(1.0776-(2.30695-(.43424-(.082433-(.008056 -.0003146*z)*z)*z)*z)*z));
 /* max |error|<.0000008 for 4<z<infinity */
}

/** @brief Returns the finite sample cdf of the Anderson-Darling
    statistics.
 
   Adaptation of function `AD` of Marsaglia & Marsaglia (2004)
   [J. Stat. Software 9(2): 1-5](https://www.jstatsoft.org/article/view/v009i02).
 
   @param[in] n an integer with the sample size
   @param[in] z a double the observed statistics value
   @returns a double Prob{W2 <= z}
*/
double AD_cdf_P(int n,double z)
{
  double v;
  double x=adinf(z);
  if(x>.8)
  {
    v=(-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)*x)*x)*x)*x)/n;
    return x+v;
  }
  double c=.01265+.1757/n;
  if(x<c)
  {
    v=x/c;
    v=sqrt(v)*(1.-v)*(49*v-102);
    return x+v*(.0037/(n*n)+.00078/n+.00006)/n;
  }
  v=(x-c)/(.8-c);
  v=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*v)*v)*v)*v)*v;
  return x+v*(.04213+.01365/n)/n;
}

int main()
{
  gsl_rng * r;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  double w90 = 1.93295783274159;
  double w95 = 2.4923671600494096;
  double w99 = 3.8781250216053948;

  printf("********************************************************\n");
  printf("Test asymptotic cdf using Marsaglia & Marsaglia values:\n");
  printf("Pr(W2 < 1.93295783274159) = %g, the theoretical value is 0.90\n", adinf(w90));
  printf("Pr(W2 < 2.4923671600494096) = %g, the theoretical value is 0.95\n", adinf(w95));
  printf("Pr(W2 < 3.8781250216053948) = %g, the theoretical value is 0.99\n", adinf(w99));
  printf("\n");
  printf("Generating a sample of size 100 with Box-Müller (mean 0 and var 1)\n");
  gsl_vector * sample = gsl_vector_alloc(100);
  for (size_t i=0; i < sample->size; i++)
    gsl_vector_set(sample,i,gsl_ran_gaussian(r,1.0));
  double W2 = AndersonDarling_W2(sample,false);
  printf("The Anderson-Darling statistics, W2, value is: %g\n", W2);
  printf("Prob(W2 <= %g) = %g\n",W2, AD_cdf_P(sample->size,W2));
  gsl_vector_free(sample);
  printf("Generating a sample of size 750 with the ratio of uniforms (mean 0 and var 1)\n");
  sample = gsl_vector_alloc(750);
  for (size_t i=0; i < sample->size; i++)
    gsl_vector_set(sample,i,gsl_ran_ugaussian_ratio_method(r));
  W2 = AndersonDarling_W2(sample,false);
  printf("The Anderson-Darling statistics, W2, value is: %g\n", W2);
  printf("Prob(W2 <= %g) = %g\n",W2, AD_cdf_P(sample->size,W2));
  gsl_vector_free(sample);
  printf("Generating a sample of size 500 with the ziggurat (mean 0 and var 1.2)\n");
  sample = gsl_vector_alloc(500);
  for (size_t i=0; i < sample->size; i++)
    gsl_vector_set(sample,i,gsl_ran_gaussian_ziggurat(r, 1.2));
  W2 = AndersonDarling_W2(sample,false);
  printf("The Anderson-Darling statistics, W2, value is: %g\n", W2);
  printf("Prob(W2 <= %g) = %g\n",W2, AD_cdf_P(sample->size,W2));
  printf("\n");
  printf("Generating now 10,000 samples of size 100 form a standard normal with the ziggurat\n");
  printf("method and computing the cdf of the W2 statistics...\n");
  double W2_cdf_sample[10000];
  sample = gsl_vector_alloc(100);
  for (size_t rep=0; rep < 10000; rep++) {
    for (size_t i=0; i < sample->size; i++)
      gsl_vector_set(sample,i,gsl_ran_gaussian_ziggurat(r, 1.0));
    W2_cdf_sample[rep] = AD_cdf_P(sample->size,AndersonDarling_W2(sample,false));
  }
  gsl_vector_free(sample);
  gsl_sort(W2_cdf_sample,1,10000);
  double decile[9];
  printf("The theoretical deciles of the cdf of the W2 statistics are:\n");
  printf("0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9.\n");
  printf("The empirical deciles of the cdf of the W2 statistics are:\n");
  for (size_t i=0; i < 9; i++) 
    decile[i] = W2_cdf_sample[(i+1)*1000-1];
  printf("%g, %g, %g, %g, %g, %g, %g, %g, %g.\n", decile[0],
    decile[1],decile[2],decile[3],decile[4],decile[5],decile[6],
    decile[7],decile[8]);
  printf("\n");
  printf("Generating now 10,000 samples of size 100 form a normal(0,1.05) with the ziggurat\n");
  printf("method and computing the cdf of the W2 statistics...\n");
  sample = gsl_vector_alloc(100);
  for (size_t rep=0; rep < 10000; rep++) {
    for (size_t i=0; i < sample->size; i++)
      gsl_vector_set(sample,i,gsl_ran_gaussian_ziggurat(r, 1.05));
    W2_cdf_sample[rep] = AD_cdf_P(sample->size,AndersonDarling_W2(sample,false));
  }
  gsl_vector_free(sample);
  gsl_sort(W2_cdf_sample,1,10000);
  printf("The deciles of the cdf of the W2 statistics are:\n");
  for (size_t i=0; i < 9; i++) 
    decile[i] = W2_cdf_sample[(i+1)*1000-1];
  printf("%g, %g, %g, %g, %g, %g, %g, %g, %g\n", decile[0],
    decile[1],decile[2],decile[3],decile[4],decile[5],decile[6],
    decile[7],decile[8]);
  printf("********************************************************\n");
  printf("\n");
  
  gsl_rng_free(r);
  
}
