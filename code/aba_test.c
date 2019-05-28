/** @file aba_test.c
 *  @brief test functions reading data from an HDF5 file
 *         and printing them to the stdout 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_vector.h>
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
#define FNAME "data_paper/data_beta_escin/DA_121219_E1.h5"
int main()
{
  // Open FILE
  char fname[] = FNAME;
  hid_t fid = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  aba * paba = aba_read_from_file(fid);
  // Close file
  H5Fclose (fid);
  aba_printf(paba);
  aba_free(paba);
  return 0;
}
