/** \file adu_vector.c
    \brief Function definitions for [adu_vector](@ref adu_vector) structures
 */
#include "abaa.h"
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
