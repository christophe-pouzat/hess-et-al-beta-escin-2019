/** \file gof.c
    \brief Functions defining goodness of fit tests.
 */

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
