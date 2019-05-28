/** \file abaa.h
 *  \brief Header file for "added buffer approach analysis" 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

/** @brief Structure holding arrays of gsl_vectors each vector contains
 *         ADU measurements from at a specific wavelength from a 
 *         specific loaction (ROI or ROB).
*/
typedef struct
{
  gsl_vector * ADU340; //!< measurements at 340 nm from ROI
  gsl_vector * ADU340B; //!< measurements at 340 nm from ROB
  gsl_vector * ADU360; //!< measurements at 360 nm from ROI
  gsl_vector * ADU360B; //!< measurements at 360 nm from ROB
  gsl_vector * ADU380; //!< measurements at 380 nm from ROI
  gsl_vector * ADU380B; //!< measurements at 380 nm from ROB
  gsl_vector * TIME; //!< time (in s) of measurements
} adu;

adu * adu_alloc(size_t n_obs);
int adu_free(adu * adu_ptr);
/** @def adu_get(adu,member,i)
 *  @brief A macro that returns value at index \a i of
 *         member \a member from \a adu structure 
 */
#define adu_get(adu,member,i) gsl_vector_get(adu->member,i)
/** @def adu_set(adu,member,i,x)
 *  @brief A macro that sets value at index \a i of
 *         member \a member from \a adu structure at \a x 
 */
#define adu_set(adu,member,i,x) gsl_vector_set(adu->member,i,x)
int adu_printf(adu * padu);

/** @brief Structure holding arrays of `adu` structures .
*/
typedef struct
{
  size_t nelt; //!< number of elements in the vector
  adu ** adu_v; //!< array of pointers to adu structures
} adu_vector;
adu_vector * adu_vector_alloc(size_t nelt);
int adu_vector_free(adu_vector * adu_vector_ptr);
int adu_vector_printf(adu_vector * padu_vector);
size_t data_get_nelt(hid_t file_id);
adu_vector * adu_vector_read_from_file(hid_t file_id);

/** @brief Structure holding dye parameters.
*/
typedef struct
{
  double R_min_hat; //!< estimated `R_min` parameter
  double R_min_se; //!< estimated `R_min` standard error
  double R_max_hat; //!< estimated `R_max` parameter
  double R_max_se; //!< estimated `R_max` standard error
  double K_eff_hat; //!< estimated `K_eff` parameter in \f$\mu{}M\f$
  double K_eff_se; //!< estimated `K_eff` standard error in \f$\mu{}M\f$
  double K_d_hat; //!< estimated `K_d` parameter in \f$\mu{}M\f$
  double K_d_se; //!< estimated `K_d` standard error in \f$\mu{}M\f$
  double pipette_concentration;//!< dye concentration in the pipette in \f$\mu{}M\f$ 
} dye;
dye dye_read_from_file(hid_t file_id);
int dye_printf(dye * pdye);

/** @brief Structure holding illumination parameters.
*/
typedef struct
{
  double T_340; //!< illumination duration at 340 nm (s)
  double T_360; //!< illumination duration at 360 nm (s)
  double T_380; //!< illumination duration at 380 nm (s)
} illumination;
illumination illumination_read_from_file(hid_t file_id);
int illumination_printf(illumination * pillumination);

/** @brief Structure holding ccd parameters.
*/
typedef struct
{
  double gain; //!< CCD chip gain
  double s2; //!< CCD chip read-out variance
  size_t P; //!< number of pixels in ROI
  size_t P_B; //!< number of pixels in ROB
} ccd;
ccd ccd_read_from_file(hid_t file_id);
int ccd_printf(ccd * pccd);

/** \brief Structure holding all the data
 */
typedef struct 
{
  adu_vector * data; //!< a pointer to an adu_vector
  dye dye; //!< dye parameters
  illumination light; //!< illumination parameters
  ccd ccd; //!< ccd chip parameters
} aba;
aba * aba_alloc();
int aba_free(aba * paba);
aba * aba_read_from_file(hid_t file_id);
int aba_printf(aba * paba);

/** \struct ratio
 *  \brief Structure holding arrays of gsl_vectors each vector contains
 *         ratiometric estimates.
*/
typedef struct
{
  gsl_vector * RATIO; //!< ratiometric estimator
  gsl_vector * RATIO_SE; //!< standard error of RATIO
  gsl_vector * TIME; //!< time (in s) of measurements
} ratio;
ratio * ratio_alloc(size_t n_obs);
int ratio_free(ratio * ratio_ptr);
/** \def ratio_get(ratio,member,i)
    \brief A macro that returns value at index \a i of
           member \a member from \a [ratio](@ref ratio) structure 
*/
#define ratio_get(ratio,member,i) gsl_vector_get(ratio->member,i)
/** \def ratio_set(ratio,member,i,x)
    \brief A macro that sets value at index \a i of
           member \a member from \a [ratio](@ref ratio) structure at \a x 
*/
#define ratio_set(ratio,member,i,x) gsl_vector_set(ratio->member,i,x)
int ratio_fprintf(FILE * fp,ratio * pratio);
ratio * ratio_est(adu * padu, dye * pdye,illumination * plight,ccd * pccd,size_t nrep);
size_t ratio_find_fit_start(ratio * pratio, double remaining_fraction, size_t baseline);
/** \struct mono_exp_fit_res
    \brief Structure holding result of a mono-exponential
    fit Y = baseline + delta * exp(-(t-t0)/tau)
 
    Depending on the input data the nonlinear least-squares
    can be weighted. That's the case when working with the
    ratiometric estimator.
 */
typedef struct
{
  size_t nobs; //!< number of observation used for the fit
  size_t baseline_length; //!< length of baseline region
  size_t fit_start; //!< first point used for the fit
  double baseline; //!< fitted baseline value
  double baseline_se; //!< standard error of the baseline
  double delta; //!< fitted jump amplitude
  double delta_se; //!< standard error of delta
  double tau; //!< fitted time constant
  double tau_se; //!< standard error of tau
  double rss; //!< RSS value
  int status; //!< solver's status when fit stops
} mono_exp_fit_res;
int mono_exp_fit_res_fprintf(FILE * fp, mono_exp_fit_res * str, ratio * pratio);
/** \struct ratio_for_fit
    \brief Structure suitable for fitting a mono-exponential decay to 
    the ratio metric estimates stored in a [ratio](@ref ratio) structure 
 */
typedef struct
{
  size_t i0; //!< index at which the decay fit starts
  size_t baseline_length; //!< baseline length
  ratio * pratio; //!< pointer to a ratio structure
} ratio_for_fit;
int ratio_residuals (const gsl_vector * x, void *data, gsl_vector * f);
void ratio_fit_callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
mono_exp_fit_res ratio_fit(ratio * pratio, size_t baseline, size_t start_fit, size_t maxit);

/** \struct ts
    \brief A two members structure (two pointeres to gsl_vectors
           holding time and amplitude from a time series.
 */
typedef struct
{
  gsl_vector * TIME; //!< a pointer to a gsl_vector holding times
  gsl_vector * AMPLITUDE; //!< a pointer to a gsl_vector holding amplitudes
} ts;

/** \struct ts_vector
    \brief Structure holding arrays of [ts](@ref ts) structures .
*/
typedef struct
{
  size_t nelt; //!< number of elements in the vector
  ts ** ts_v; //!< array of pointers to [ts](@ref ts) structures
} ts_vector;

ts * ts_alloc(size_t n_obs);

int ts_free(ts * pts);

int ts_fprintf(FILE * stream, ts * pts);

ts_vector * ts_vector_alloc(size_t nelt);

int ts_vector_free(ts_vector * pts_vector);

int ts_vector_fprintf(FILE * stream, ts_vector * pts_vector);

ts_vector * fura_est(aba * paba);
