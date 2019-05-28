/** \file ratio_test.c
    \brief Program testing the ratiometric functions
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
/** @brief Allocates an [adu](@ref adu)
 *
 *  The function allocates memory for an [adu](@ref adu) structure
 *
 *  @param[in] n_obs the number of measurements / obserations 
 *  @returns a pointer to an allocated [adu](@ref adu)
*/
adu * adu_alloc(size_t n_obs) {
  adu * res = malloc(sizeof(adu));
  res->ADU340 = gsl_vector_alloc(n_obs);
  res->ADU340B = gsl_vector_alloc(n_obs);
  res->ADU360 = gsl_vector_alloc(n_obs);
  res->ADU360B = gsl_vector_alloc(n_obs);
  res->ADU380 = gsl_vector_alloc(n_obs);
  res->ADU380B = gsl_vector_alloc(n_obs);
  res->TIME = gsl_vector_alloc(n_obs);
  return res;
}
/** @brief Frees an [adu](@ref adu)
 
    @param[in,out] adu_ptr a pointer to an allocated [adu](@ref adu) structure
    @returns 0 if everything goes fine 
*/
int adu_free(adu * adu_ptr) {
  gsl_vector_free(adu_ptr->ADU340);
  gsl_vector_free(adu_ptr->ADU340B);
  gsl_vector_free(adu_ptr->ADU360);
  gsl_vector_free(adu_ptr->ADU360B);
  gsl_vector_free(adu_ptr->ADU380);
  gsl_vector_free(adu_ptr->ADU380B);
  gsl_vector_free(adu_ptr->TIME);
  free(adu_ptr);
  return 0;
}
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
/** @brief Structure holding arrays of `adu` structures .
*/
typedef struct
{
  size_t nelt; //!< number of elements in the vector
  adu ** adu_v; //!< array of pointers to adu structures
} adu_vector;
/** @brief Allocates an adu_vector
 *
 *  The function allocates memory for an adu_vector structure
 *
 *  @param[in] nelt the number of stimulation 
 *  @returns a pointer to an allocated adu_vector
*/
adu_vector * adu_vector_alloc(size_t nelt) {
  adu_vector * res = malloc(sizeof(adu_vector));
  res->nelt = nelt;
  res->adu_v = malloc(nelt*sizeof(adu));
  return res;
}
/** @brief Frees an adu_vector
 *
 *  @param[in,out] adu_vector_ptr a pointer to an allocated adu_vector structure
 *  @returns 0 if everything goes fine 
*/
int adu_vector_free(adu_vector * adu_vector_ptr) {
  for (size_t d_idx=0; d_idx<adu_vector_ptr->nelt; d_idx++) 
    adu_free(adu_vector_ptr->adu_v[d_idx]);
  free(adu_vector_ptr->adu_v);
  free(adu_vector_ptr);
  return 0;
}
/** @brief Returns number of DataSets in Group DATA
 *
 *  The number returned equal 1 plus the number of stim
 *
 *  @param[in] file_id an HFD5 file identifier
 *  @return a size_t, the number of DataSets 
 */
size_t data_get_nelt(hid_t file_id) {
  char dset[] = "/DATA";
  hid_t gid = H5Gopen(file_id,dset,H5P_DEFAULT);
  // Get info on group DATA
  H5G_info_t group_info;
  H5Gget_info(gid, &group_info);
  size_t n_elt = (size_t) group_info.nlinks;
  // Close Group
  H5Gclose(gid);
  return n_elt;
}
/** @brief Allocates and initializes an `adu_vector` structure
 *         read from Group DATA in a file
 *
 *  @param[in] file_id HDF5 file identifier
 *  @return an allocated and initialized pointer to an adu_vector structure
 */
adu_vector * adu_vector_read_from_file(hid_t file_id) {
  char *dsets[] = {"/DATA/load","/DATA/stim1","/DATA/stim2",
		   "/DATA/stim3","/DATA/stim4","/DATA/stim5",
		   "/DATA/stim6","/DATA/stim7","/DATA/stim8"};
  char STIM[256],DELTA[256],OFFSET[256];
  size_t n_elt = data_get_nelt(file_id);
  if (n_elt > 9) {
    fprintf(stderr,"Too many data sets (>9).\n");
    return NULL;
  }
  adu_vector * data = adu_vector_alloc(n_elt);
  for (size_t d_idx=0; d_idx<n_elt; d_idx++) { 
    STIM[0] = '\0';
    strcat(STIM,dsets[d_idx]);
    strcat(STIM,"/ADU");
    // load DataSet
    hsize_t dims[2];
    H5LTget_dataset_info(file_id,STIM,dims,NULL,NULL);
    size_t nobs = (size_t) dims[0];
    size_t ncol = (size_t) dims[1];
    int *ADU = malloc(nobs*ncol*sizeof(int));
    H5LTread_dataset_int(file_id,STIM,ADU);
    DELTA[0] = '\0';
    strcat(DELTA,dsets[d_idx]);
    strcat(DELTA,"/TIME_DELTA");
    double delta;
    H5LTread_dataset_double(file_id,DELTA,&delta);
    OFFSET[0] = '\0';
    strcat(OFFSET,dsets[d_idx]);
    strcat(OFFSET,"/TIME_OFFSET");
    double offset;
    H5LTread_dataset_double(file_id,OFFSET,&offset);
    data->adu_v[d_idx] = adu_alloc(nobs);
    for (size_t i=0; i<nobs; i++) {
      adu_set((data->adu_v[d_idx]),TIME,i,offset+delta*((double) ADU[i*ncol]));
      adu_set((data->adu_v[d_idx]),ADU340,i,(double) ADU[i*ncol + 1]);
      adu_set((data->adu_v[d_idx]),ADU340B,i,(double) ADU[i*ncol + 2]);
      adu_set((data->adu_v[d_idx]),ADU360,i,(double) ADU[i*ncol + 3]);
      adu_set((data->adu_v[d_idx]),ADU360B,i,(double) ADU[i*ncol + 4]);
      adu_set((data->adu_v[d_idx]),ADU380,i,(double) ADU[i*ncol + 5]);
      adu_set((data->adu_v[d_idx]),ADU380B,i,(double) ADU[i*ncol + 6]);
    }
    free(ADU);
  }
  return data;
}
/** @brief Prints [adu](@ref adu) content to stdout
 
    @param[in] padu a pointer to an [adu](@ref adu) structure
    @return 0 if everything goes fine
 */
int adu_printf(adu * padu) {
  size_t nobs=(padu)->TIME->size;
  printf("#    Time   ADU340  ADU340B   ADU360  ADU360B   ADU380  ADU380B\n");
  for (size_t i=0; i<nobs; i++) {
    printf("%9.9g %8d %8d %8d %8d %8d %8d\n",
	   adu_get((padu),TIME,i),
	   (int) adu_get(padu,ADU340,i),
	   (int) adu_get(padu,ADU340B,i),
	   (int) adu_get(padu,ADU360,i),
	   (int) adu_get(padu,ADU360B,i),
	   (int) adu_get(padu,ADU380,i),
	   (int) adu_get(padu,ADU380B,i));
  }
  printf("\n\n");
  return 0;
}
/** @brief Prints adu_vector content to stdout
 *
 *  @param[in] padu_vector a pointer to an adu_vector structure
 *  @return 0 if everything goes fine
 */
int adu_vector_printf(adu_vector * padu_vector) {
  for (size_t d_idx=0; d_idx<padu_vector->nelt; d_idx++) {
    size_t nobs=(padu_vector->adu_v[d_idx])->TIME->size;
    if (d_idx == 0) {
      printf("# Loading curve with %d elements\n", (int) nobs);
    } else {
      printf("# Stim %d with %d elements\n", (int) d_idx, (int) nobs);
    }
    adu_printf(padu_vector->adu_v[d_idx]);
  }
  return 0;
}
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
/** @brief Returns a `dye` structure
 *         read from Group DYE in a file
 *
 *  @param[in] file_id HDF5 file identifier
 *  @return a `dye` structure
 */
dye dye_read_from_file(hid_t file_id) {
  dye res;
  H5LTread_dataset_double(file_id,"/DYE/R_min_hat",&res.R_min_hat);
  H5LTread_dataset_double(file_id,"/DYE/R_min_se",&res.R_min_se);
  H5LTread_dataset_double(file_id,"/DYE/R_max_hat",&res.R_max_hat);
  H5LTread_dataset_double(file_id,"/DYE/R_max_se",&res.R_max_se);
  H5LTread_dataset_double(file_id,"/DYE/K_eff_hat",&res.K_eff_hat);
  H5LTread_dataset_double(file_id,"/DYE/K_eff_se",&res.K_eff_se);
  H5LTread_dataset_double(file_id,"/DYE/K_d_hat",&res.K_d_hat);
  H5LTread_dataset_double(file_id,"/DYE/K_d_se",&res.K_d_se);
  H5LTread_dataset_double(file_id,"/DYE/pipette_concentration",&res.pipette_concentration);
  return res;
}
/** @brief Prints dye structure content to stdout
 *
 *  @param[in] pdye a pointer to a `dye` structure
 *  @return 0 if everything goes fine
 */
int dye_printf(dye * pdye) {
  printf("# DYE parameters\n");
  printf("# R_min_hat: %g\n",pdye->R_min_hat);
  printf("# R_min_se: %g\n",pdye->R_min_se);
  printf("# R_max_hat: %g\n",pdye->R_max_hat);
  printf("# R_max_se: %g\n",pdye->R_max_se);
  printf("# K_eff_hat: %g\n",pdye->K_eff_hat);
  printf("# K_eff_se: %g\n",pdye->K_eff_se);
  printf("# K_d_hat: %g\n",pdye->K_d_hat);
  printf("# K_d_se: %g\n",pdye->K_d_se);
  printf("# pipette_concentration: %g\n\n",pdye->pipette_concentration);
  return 0;
}
/** @brief Structure holding illumination parameters.
*/
typedef struct
{
  double T_340; //!< illumination duration at 340 nm (s)
  double T_360; //!< illumination duration at 360 nm (s)
  double T_380; //!< illumination duration at 380 nm (s)
} illumination;
/** @brief Returns a `illumination` structure
 *         read from Group ILLUMINATION in a file
 *
 *  @param[in] file_id HDF5 file identifier
 *  @return a `illumination` structure
 */
illumination illumination_read_from_file(hid_t file_id) {
  illumination res;
  H5LTread_dataset_double(file_id,"/ILLUMINATION/T_340",&res.T_340);
  H5LTread_dataset_double(file_id,"/ILLUMINATION/T_360",&res.T_360);
  H5LTread_dataset_double(file_id,"/ILLUMINATION/T_380",&res.T_380);
  return res;
}
/** @brief Prints illumination structure content to stdout
 *
 *  @param[in] pillumination a pointer to a `illumination` structure
 *  @return 0 if everything goes fine
 */
int illumination_printf(illumination * pillumination) {
  printf("# ILLUMINATION parameters\n");
  printf("# T_340: %g\n",pillumination->T_340);
  printf("# T_360: %g\n",pillumination->T_360);
  printf("# T_380: %g\n\n",pillumination->T_380);
  return 0;
}
/** @brief Structure holding ccd parameters.
*/
typedef struct
{
  double gain; //!< CCD chip gain
  double s2; //!< CCD chip read-out variance
  size_t P; //!< number of pixels in ROI
  size_t P_B; //!< number of pixels in ROB
} ccd;
/** @brief Returns a `ccd` structure
 *         read from Group CCD in a file
 *
 *  @param[in] file_id HDF5 file identifier
 *  @return a `ccd` structure
 */
ccd ccd_read_from_file(hid_t file_id) {
  ccd res;
  H5LTread_dataset_double(file_id,"/CCD/GAIN",&res.gain);
  double value;
  H5LTread_dataset_double(file_id,"/CCD/S_RO",&value);
  res.s2 = value*value;
  int n;
  H5LTread_dataset_int(file_id,"/CCD/P",&n);
  res.P = (size_t) n;
  H5LTread_dataset_int(file_id,"/CCD/P_B",&n);
  res.P_B = (size_t) n;
  return res;
}
/** @brief Prints `ccd` structure content to stdout
 *
 *  @param[in] pccd a pointer to a `ccd` structure
 *  @return 0 if everything goes fine
 */
int ccd_printf(ccd * pccd) {
  printf("# CCD parameters\n");
  printf("# Gain: %g\n",pccd->gain);
  printf("# Read-out variance: %g\n",pccd->s2);
  printf("# P: %d\n", (int) pccd->P);
  printf("# P_B: %d\n\n", (int) pccd->P_B);
  return 0;
}
/** \brief Structure holding all the data
 */
typedef struct 
{
  adu_vector * data; //!< a pointer to an adu_vector
  dye dye; //!< dye parameters
  illumination light; //!< illumination parameters
  ccd ccd; //!< ccd chip parameters
} aba;
/** \brief Allocates an [aba](@ref aba) structure
 
    The function allocates memory for an [aba](@ref aba) structure
  
    @returns a pointer to an allocated [aba](@ref aba)
*/
aba * aba_alloc()
{
  aba * res = malloc(sizeof(aba));
  return res;
}
/** @brief Frees an [aba](@ref aba)
 
    @param[in,out] aba_ptr a pointer to an allocated [aba](@ref aba) structure
    @returns 0 if everything goes fine 
*/
int aba_free(aba * paba)
{
  adu_vector_free(paba->data);
  free(paba);
  return 0;
}
/** \brief Reads the content of an opened HDF5 file
           and stores the result in an [aba](@ref aba)

    \param[in] file_id pointer to an opened HDF5 file
    \return a pointer to an allocated and intialized [aba](@ref aba)
 */
aba * aba_read_from_file(hid_t file_id)
{
  aba * res = aba_alloc();
  res->data = adu_vector_read_from_file(file_id);
  res->dye = dye_read_from_file(file_id);
  res->light = illumination_read_from_file(file_id);
  res->ccd = ccd_read_from_file(file_id);
  return res;
}
/** @brief Prints [aba](@ref aba) content to stdout
 
    @param[in] paba a pointer to an [aba](@ref aba) structure
    @return 0 if everything goes fine
 */
int aba_printf(aba * paba)
{
  dye_printf(&(paba->dye));
  ccd_printf(&(paba->ccd));
  illumination_printf(&(paba->light));
  adu_vector_printf(paba->data);
  return 0;
}
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
#define FNAME "data_paper/data_beta_escin/DA_121219_E1.h5"
int main()
{
  // Open FILE
  char fname[] = FNAME;
  hid_t fid = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  aba * paba = aba_read_from_file(fid);
  // Close file
  H5Fclose (fid);
  ratio * pratio = ratio_est(paba->data->adu_v[1],&paba->dye,&paba->light,&paba->ccd,1000);
  aba_free(paba);
  ratio_fprintf(stdout,pratio);
  printf("\n\n");
  size_t start_pos = ratio_find_fit_start(pratio,0.5,15);
  mono_exp_fit_res res = ratio_fit(pratio,15,start_pos,25);
  mono_exp_fit_res_fprintf(stdout,&res,pratio);
  ratio_free(pratio);
  return 0;
}
