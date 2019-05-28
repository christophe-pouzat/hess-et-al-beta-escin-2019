/** \file ratio.c
    \brief Functions to work with the ratiometric estimator.
 */
#include "abaa.h"
/** \fn ratio * ratio_alloc(size_t n_obs)
    \brief Allocates a [ratio](@ref ratio) structure
 
    The function allocates memory for an [ratio](@ref ratio) structure
 
    \param[in] n_obs the number of measurements / obserations 
    \returns a pointer to an allocated ratio
*/
ratio * ratio_alloc(size_t n_obs) {
  ratio * res = malloc(sizeof(ratio));
  res->RATIO = gsl_vector_alloc(n_obs);
  res->RATIO_SE = gsl_vector_alloc(n_obs);
  res->TIME = gsl_vector_alloc(n_obs);
  return res;
}

/** \fn int ratio_free(ratio * ratio_ptr)
    \brief Frees a [ratio](@ref ratio) structure
 
    \param ratio_ptr a pointer to an allocated [ratio](@ref ratio) structure
    \returns 0 if everything goes fine 
*/
int ratio_free(ratio * ratio_ptr) {
  gsl_vector_free(ratio_ptr->RATIO);
  gsl_vector_free(ratio_ptr->RATIO_SE);
  gsl_vector_free(ratio_ptr->TIME);
  free(ratio_ptr);
  return 0;
}

/** \fn int ratio_fprintf(FILE* fp, ratio * pratio)
    \brief Prints [ratio](@ref ratio) content to fp
 
    \param[in] fp file pointer
    \param[in] pratio a pointer to a [ratio](@ref ratio) structure
    \return 0 if everything goes fine
*/
int ratio_fprintf(FILE* fp, ratio * pratio) {
  size_t nobs=(pratio)->TIME->size;
  fprintf(fp,"# Ratiometric estimator:\n");
  fprintf(fp,"#    Time    Ca (est)    Ca (se)\n");
  for (size_t i=0; i<nobs; i++) {
    double time = ratio_get((pratio),TIME,i);
    double ratio = ratio_get((pratio),RATIO,i);
    double ratio_se = ratio_get((pratio),RATIO_SE,i);
    fprintf(fp,"%9.9g %9.9g %9.9g\n",
	   time, ratio, ratio_se);
  }
  fprintf(fp,"\n\n");
  return 0;
}

/** \fn ratio * ratio_est(adu * padu, dye * pdye, illumination * plight, ccd * pccd, size_t nrep)
    \brief Computes ratiometric estimator and estimates its standard error

    \param[in] padu a pointer to an [adu](@ref adu) structure
    \param[in] pdye a pointer to a [dye](@ref dye) structure
    \param[in] plight a pointer to an [illumination](@ref illumination) structure
    \param[in] nrep the number of replicates used in the MC estimation of the SE
    \return a pointer to an initialized [ratio](@ref ratio) structure
 */
ratio * ratio_est(adu * padu,
		  dye * pdye,
		  illumination * plight,
		  ccd * pccd,
		  size_t nrep)
{
  // get the dye parameters
  double Rmin = pdye->R_min_hat;
  double Rmax = pdye->R_max_hat;
  double Keff = pdye->K_eff_hat;
  // get the illumination parameters
  double T340 = plight->T_340;
  double T380 = plight->T_380;
  // get the CCD chip parameters
  double nP = (double) pccd->P;
  double nPB = (double) pccd->P_B;

  // macro that returns the "ratio"
#define R(ADU340,ADU340B,ADU380,ADU380B)			\
  ((ADU340/nP)-(ADU340B/nPB))*T380/\
  ((ADU380/nP)-(ADU380B/nPB))/T340
  // macro the returns the estimated free [Ca2+]
  // for a given "ratio"
#define Ca(r) Keff*(r-Rmin)/(Rmax-r)

  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng_env_setup();
  gsl_rng * rng = gsl_rng_alloc (T);
  size_t nobs = (padu)->TIME->size;
  ratio * res = ratio_alloc(nobs);
  // do the MC standard error estimation
  for (size_t i=0; i<nobs; i++) {
    double adu340 = adu_get(padu,ADU340,i);
    double adu340B = adu_get(padu,ADU340B,i);
    double adu380 = adu_get(padu,ADU380,i);
    double adu380B = adu_get(padu,ADU380B,i);
    double r = R(adu340,adu340B,adu380,adu380B);
    ratio_set(res,TIME,i,adu_get(padu,TIME,i));
    ratio_set(res,RATIO,i,Ca(r));
    double g = pccd->gain;
    double s2 = g*g*nP*pccd->s2;
    double s2B = g*g*nPB*pccd->s2;
    double Ca_rep[nrep];
    for (size_t j=0; j<nrep; j++) {
      double adu340r = adu340+gsl_ran_gaussian_ziggurat(rng,sqrt(g*adu340+s2));
      double adu340Br = adu340B+gsl_ran_gaussian_ziggurat(rng,sqrt(g*adu340B+s2B));
      double adu380r = adu380+gsl_ran_gaussian_ziggurat(rng,sqrt(g*adu380+s2));
      double adu380Br = adu380B+gsl_ran_gaussian_ziggurat(rng,sqrt(g*adu380B+s2B));
      r = R(adu340r,adu340Br,adu380r,adu380Br);
      Ca_rep[j] = Ca(r);
    }
    ratio_set(res,RATIO_SE,i,gsl_stats_sd(Ca_rep,1,nrep));
  }
  gsl_rng_free (rng);
  return res;
}

/** \fn size_t ratio_find_fit_start(ratio * pratio, double remaining_fraction, size_t baseline)
    \brief Locates the first point of a transient whose amplitude is smaller than the baseline 
    plus a fraction of the stimulation induced jump.
    
    \param[in] pratio a pointer to an initialized [ratio](@ref ratio) structure
    \param[in] remaining_fraction the remaining fraction of the jump amplitude
    from which to start
    \param[in] baseline the length (in sample points) of the baseline
    \return a `size_t` with the searched location 
 */
size_t ratio_find_fit_start(ratio * pratio, double remaining_fraction, size_t baseline)
{
  double baseline_mean = 0.0;
  for (size_t i=0; i<baseline; i++)
    baseline_mean += gsl_vector_get(pratio->RATIO,i);
  baseline_mean /= baseline;
  size_t max_idx = gsl_vector_max_index(pratio->RATIO);
  double DeltaCa = gsl_vector_get(pratio->RATIO,max_idx)-baseline_mean;
  double threshold = DeltaCa*remaining_fraction+baseline_mean;
  size_t i=max_idx+1;
  while ((gsl_vector_get(pratio->RATIO,i) > threshold) & (i < pratio->RATIO->size))
    i++;
  return i;
}

/** \fn int ratio_residuals (const gsl_vector * x, void *data, gsl_vector * f)
    \brief Returns the (weighted) residuals of a mono-exponential fit suitable 
    for use with the `gsl` [nonlinear least-squares solver](https://www.gnu.org/software/gsl/doc/html/nls.html#providing-the-function-to-be-minimized)

    \param[in] x pointer to a [`gsl_vector`](https://www.gnu.org/software/gsl/doc/html/vectors.html#vectors) 
    holding the present model parameter values
    \param[in] data pointer to a structure holding the "data" to fit
    \param[out] f pointer to a [`gsl_vector`](https://www.gnu.org/software/gsl/doc/html/vectors.html#vectors) 
    where the residuals get stored.
    \returns `GSL_SUCCESS` if everything goes fine
 */
int ratio_residuals (const gsl_vector * x,
		     void *data, 
		     gsl_vector * f)
{
  size_t i0 = ((ratio_for_fit *)data)->i0;
  size_t baseline_length = ((ratio_for_fit *)data)->baseline_length;
  ratio *pratio = ((ratio_for_fit *)data)->pratio;
  gsl_vector * time = pratio->TIME;
  gsl_vector * Ca = pratio->RATIO;
  gsl_vector * SE = pratio->RATIO_SE;
  size_t n = time->size;

  double baseline = gsl_vector_get (x, 0);
  double delta = gsl_vector_get (x, 1);
  double tau = gsl_vector_get (x, 2);

  size_t i;

  for (i=0; i<baseline_length; i++)
    gsl_vector_set (f, i, (baseline - gsl_vector_get(Ca,i))/gsl_vector_get(SE,i));

  double t0 = gsl_vector_get(time,i0);
  for (size_t j=i0; j<n; j++){
    /* Model Yj = delta * exp(-(tj-t0)/tau) + baseline */
    double dt = gsl_vector_get(time,j) - t0;
    double Yj = delta * exp (-dt/tau) + baseline;
    gsl_vector_set (f, i, (Yj - gsl_vector_get(Ca,j))/gsl_vector_get(SE,j));
    i++;
  }

  return GSL_SUCCESS;
}

/** \fn void ratio_fit_callback(const size_t iter, void *params,
    const gsl_multifit_nlinear_workspace *w)
    \brief A callback function printing progress during 
    [nonlinear least-squares fitting](https://www.gnu.org/software/gsl/doc/html/nls.html#high-level-driver)

    \param[in] iter the current iteration
    \param[in] params set to `NULL` when used but required by the [solver](https://www.gnu.org/software/gsl/doc/html/nls.html#high-level-driver)
    \param[in] w a pointer to the [`gsl_multifit_nlinear_workspace`](https://www.gnu.org/software/gsl/doc/html/nls.html#initializing-the-solver)
*/
void ratio_fit_callback(const size_t iter,
			void *params,
			const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  
  fprintf(stderr, "iter %2zu: baseline = %.4f, delta = %.4f, tau = %.4f, RSS = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
          gsl_pow_2(gsl_blas_dnrm2(f)));
}

/** \fn int mono_exp_fit_res_fprintf(FILE * fp, mono_exp_fit_res * str, ratio * pratio)
    \brief Prints content of a [mono_exp_fit_res](@ref mono_exp_fit_res) structure 
    together with the fitted data

    \param[in] fp pointer to the file where printing is performed
    \param[in] str pointer to the [mono_exp_fit_res](@ref mono_exp_fit_res) 
    structure to print
    \param[in] pratio pointer to a [ratio](@ref ratio) structure containing
    the fitted data
    \returns 0 if everything goes fine
 */
int mono_exp_fit_res_fprintf(FILE * fp, mono_exp_fit_res * str, ratio * pratio)
{
  fprintf(fp,"# Fitted model Ca = baseline+delta*exp(-(t-t0)/tau)\n");
  fprintf(fp,"# nobs = %d\n", (int) str->nobs);
  fprintf(fp,"# number of degrees of freedom = %d\n", ((int) str->nobs)-3);
  fprintf(fp,"# baseline length = %d\n", (int) str->baseline_length);
  fprintf(fp,"# fit started from point %d\n", (int) str->fit_start);
  fprintf(fp,"# estimated baseline %g and standard error %g\n", str->baseline,str->baseline_se);
  fprintf(fp,"# estimated delta %g and standard error %g\n", str->delta,str->delta_se);
  fprintf(fp,"# estimated tau %g and standard error %g\n", str->tau,str->tau_se);
  fprintf(fp,"# residual sum of squares: %g\n", str->rss);
  double dof = (((double) str->nobs)-3.0);
  fprintf(fp,"# RSS per degree of freedom: %g\n", str->rss/dof);
  double p_larger = gsl_cdf_chisq_Q(str->rss,dof);
  fprintf(fp,"# Probability of observing a larger of equal RSS per DOF under the null hypothesis: %g\n",p_larger);
  if (p_larger < 0.01)
    fprintf(fp,"# WARNING: THE FIT IS NOT GOOD!\n");
  fprintf(fp,"\n");
  fprintf(fp,"# rss per degree of freedom: %g\n\n", str->rss/(((double) str->nobs)-3.0));
  fprintf(fp,"# Time    Ca  Prediction  Residual\n");
  for (size_t i=0; i<str->baseline_length; i++) {
    double time = gsl_vector_get(pratio->TIME,i);
    double y = gsl_vector_get(pratio->RATIO,i);
    double pred = str->baseline;
    double resid = (y-pred)/gsl_vector_get(pratio->RATIO_SE,i);
    fprintf(fp,"%g %g %g %g\n",time,y,pred,resid);
  }
  fprintf(fp,"\n");
  double t0 = gsl_vector_get(pratio->TIME,str->fit_start);
  for (size_t i=str->fit_start; i<pratio->TIME->size; i++) {
    double time = gsl_vector_get(pratio->TIME,i);
    double dt = time-t0;
    double y = gsl_vector_get(pratio->RATIO,i);
    double pred = str->baseline;
    pred += str->delta*exp(-dt/str->tau);
    double resid = (y-pred)/gsl_vector_get(pratio->RATIO_SE,i);
    fprintf(fp,"%g %g %g %g\n",time,y,pred,resid);
  }
  return 0;
}

/** \fn mono_exp_fit_res ratio_fit(ratio * pratio, size_t baseline,
    size_t start_fit, size_t maxit)
    \brief Fits a mono-exponential to the decaying part of a stimulation

    \param[in] pratio a pointer to a [ratio](@ref ratio) structure
    holding the "data"
    \param[in] baseline baseline length
    \param[in] start_fit the sampling point from which we start
    fitting the decay
    \param[in] the maximal number of iterations
    \returns a [mono_exp_fit_res](@ref mono_exp_fit_res) structure 
    holding the results
 */
mono_exp_fit_res ratio_fit(ratio * pratio,
			   size_t baseline,
			   size_t start_fit,
			   size_t maxit)
{
  // nobs: number of observations to fit
  const size_t nobs = baseline+pratio->TIME->size-start_fit;
  // p: number of model parameters
  const size_t p = 3;
  // data: a ratio_for_fit structure with the data and parameters
  //       required for computing the residuals
  ratio_for_fit data = {.i0=start_fit,
			.baseline_length=baseline,
			.pratio=pratio};

  // par: contains the estimated model parameters baseline,
  //      delta and tau
  gsl_vector * par = gsl_vector_alloc(p);
  // par initialization
  gsl_vector_set(par,0, gsl_vector_get(pratio->RATIO,0));
  gsl_vector_set(par,1, gsl_vector_get(pratio->RATIO,start_fit)-\
                 gsl_vector_get(pratio->RATIO,0));
  gsl_vector_set(par,2, gsl_vector_get(pratio->TIME,pratio->TIME->size-1)-\
  	       gsl_vector_get(pratio->TIME,start_fit));

  // define the function to be minimized
  gsl_multifit_nlinear_fdf fdf = {.f = ratio_residuals,
  				.df = NULL,  // set to NULL for finite-difference Jacobian
  				.fvv = NULL, // not using geodesic acceleration 
  				.n = nobs,
  				.p = p,
  				.params = &data};

  // allocate solver workspace with default parameters
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
  gsl_multifit_nlinear_workspace *w = gsl_multifit_nlinear_alloc (T, &fdf_params, nobs, p);
  
  // initialize solver with starting point and weights
  gsl_multifit_nlinear_init (par, &fdf, w);
  
  // compute initial rss value
  gsl_vector *residuals = gsl_multifit_nlinear_residual(w);
  double rss0;
  gsl_blas_ddot(residuals, residuals, &rss0);

  // solve the system with a maximum of maxit iterations
  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;
  int info;
  int status = gsl_multifit_nlinear_driver(maxit, xtol, gtol, ftol,
  					 ratio_fit_callback, NULL,
  					 &info, w);  
  
  gsl_matrix *J;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  // compute covariance of best fit parameters
  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar (J, 0.0, covar);
  
  // compute final cost
  double rss;
  gsl_blas_ddot(residuals, residuals, &rss);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
  
  // print end info
  fprintf(stderr, "Fitted model Ca = baseline+delta*exp(-(t-t0)/tau)\n");
  fprintf(stderr, "Summary from method '%s/%s'\n",
  	gsl_multifit_nlinear_name(w),
  	gsl_multifit_nlinear_trs_name(w));
  fprintf(stderr, "number of iterations: %zu\n",
  	gsl_multifit_nlinear_niter(w));
  fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",
  	(info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial RSS = %f\n", rss0);
  fprintf(stderr, "final   RSS = %f\n", rss);
  fprintf(stderr, "\n");
  fprintf(stderr,"Number of observation: %d\n", (int) nobs);
  fprintf(stderr,"Number of degrees of freedom: %d\n", ((int) nobs)-3);
  fprintf(stderr,"Baseline length: %d\n", (int) baseline);
  fprintf(stderr,"Fit started from point %d\n", (int) start_fit);
  fprintf(stderr,"Estimated baseline %g and standard error %g\n", FIT(0),ERR(0));
  fprintf(stderr,"Estimated delta %g and standard error %g\n", FIT(1),ERR(1));
  fprintf(stderr,"Estimated tau %g and standard error %g\n", FIT(2),ERR(2));
  double dof = (((double) nobs)-3.0);
  fprintf(stderr,"RSS per degree of freedom: %g\n", rss/dof);
  double p_larger = gsl_cdf_chisq_Q(rss,dof);
  fprintf(stderr,"Probability of observing a larger of equal RSS per DOF under the null hypothesis: %g\n",p_larger);
  if (p_larger < 0.01)
    fprintf(stderr,"WARNING: THE FIT IS NOT GOOD!\n");
  fprintf(stderr,"\n");  

  // Prepare the output  
  mono_exp_fit_res res = {.nobs=nobs,
                          .baseline_length=baseline, 
			  .fit_start=start_fit,
			  .baseline=FIT(0),
			  .baseline_se=ERR(0),
			  .delta=FIT(1),
			  .delta_se=ERR(1),
			  .tau=FIT(2),
			  .tau_se=ERR(2),
			  .rss=rss,
			  .status=status};
  // Free allocated memory
  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);
  gsl_vector_free (par);

  return res;
}
