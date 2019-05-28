/** \file fura.c
    \brief Function definitions for [Fura] related calculations
 */
#include "abaa.h"

/** \brief Allocates a [ts](@ref ts)
 
    The function allocates memory for a [ts](@ref ts) structure
    
    \param[in] n_obs the number of obserations 
    \returns a pointer to an allocated [ts](@ref ts)
*/
ts * ts_alloc(size_t n_obs) {
  ts * res = malloc(sizeof(ts));
  res->TIME = gsl_vector_alloc(n_obs);
  res->AMPLITUDE = gsl_vector_alloc(n_obs);
  return res;
}

/** @fn int ts_free(ts * pts)
    @brief Frees memory taken up by a [ts](@ref ts)
    
    @param[in,out] ptd a pointer to a [ts](@ref ts)
    @returns 0 if everything goes fine
*/
int ts_free(ts * pts) {
  gsl_vector_free(pts->TIME);
  gsl_vector_free(pts->AMPLITUDE);
  free(pts);
  return 0;
}

/** @brief Prints  [ts](@ref ts) content to `stream`
 
    @param[in] stream a pointer to an output "file"
    @param[in] pts a pointer to a [ts](@ref ts) structure
    @return 0 if everything goes fine
 */
int ts_fprintf(FILE * stream, ts * pts) {
  size_t nobs=(pts)->TIME->size;
  fprintf(stream,"#    Time   AMPLITUDE\n");
  for (size_t i=0; i<nobs; i++) {
    fprintf(stream,
	    "%9.9g %9.9g\n",
	    adu_get((pts),TIME,i),
	    adu_get((pts),AMPLITUDE,i));
  }
  fprintf(stream,"\n\n");
  return 0;
}

/** @brief Allocates an [ts_vector](@ref ts_vector)
 *
 *  The function allocates memory for a [ts_vector](@ref ts_vector) structure
 *
 *  @param[in] nelt the number of stimulation 
 *  @returns a pointer to an allocated [ts_vector](@ref ts_vector)
*/
ts_vector * ts_vector_alloc(size_t nelt) {
  ts_vector * res = malloc(sizeof(ts_vector));
  res->nelt = nelt;
  res->ts_v = malloc(nelt*sizeof(ts));
  return res;
}

/** @brief Frees a [ts_vector](@ref ts_vector)
 *
 *  @param[in,out] pts_vector a pointer to an allocated [ts_vector](@ref ts_vector) structure
 *  @returns 0 if everything goes fine 
*/
int ts_vector_free(ts_vector * pts_vector) {
  for (size_t d_idx=0; d_idx<pts_vector->nelt; d_idx++) 
    ts_free(pts_vector->ts_v[d_idx]);
  free(pts_vector->ts_v);
  free(pts_vector);
  return 0;
}

/** @brief Prints [ts_vector](@ref ts_vector) content to `stream`
 
    @param[in] stream a pointer to an output "file"
    @param[in] pts_vector a pointer to a [ts_vector](@ref ts_vector) structure
    @return 0 if everything goes fine
 */
int ts_vector_fprintf(FILE * stream, ts_vector * pts_vector) {
  for (size_t d_idx=0; d_idx<pts_vector->nelt; d_idx++) {
    size_t nobs=(pts_vector->ts_v[d_idx])->TIME->size;
    if (d_idx == 0) {
      fprintf(stream,"# Loading curve with %d elements\n", (int) nobs);
    } else {
      fprintf(stream,"# Stim %d with %d elements\n", (int) d_idx, (int) nobs);
    }
    ts_fprintf(stream,pts_vector->ts_v[d_idx]);
  }
  return 0;
}

/** \fn ts_vector * fura_est(aba * paba)
    \brief Get [Fura] time series from [aba](@ref aba) structure

    __It is assumed that when the maximal background subtracted fluorescence at
    360 nm is observed, the fura concentration in the cell and in the pipette
    are identical (no extrapolation based on a fit is performed).__

    \param[in] paba a pointer to an [aba](@ref aba) structure
    \return a pointer to an allocated and initialized [ts_vector](@ref ts_vector)
 */
ts_vector * fura_est(aba * paba) {
  // Get the number of stim + 1 in paba
  size_t nelt=paba->data->nelt;
  // Allocates result
  ts_vector *res=ts_vector_alloc(nelt);
  // Get [Fura] in the pipette
  double F_p = paba->dye.pipette_concentration;
  // Get the inverse of the number of pixels in ROI
  double inv_P = 1.0/(double)paba->ccd.P;
  // Get the inverse of the number of pixels in ROB
  double inv_P_B = 1.0/(double)paba->ccd.P_B;
  double max_adu;
  for (size_t i=0; i<nelt; i++) {
    gsl_vector *adu360 = paba->data->adu_v[i]->ADU360;
    gsl_vector *adu360B = paba->data->adu_v[i]->ADU360B;
    gsl_vector *time = paba->data->adu_v[i]->TIME;
    res->ts_v[i] = ts_alloc(adu360->size);
    for (size_t j=0; j<adu360->size; j++) {
      double fura_c = gsl_vector_get(adu360,j)*inv_P-
	gsl_vector_get(adu360B,j)*inv_P_B;
      gsl_vector_set(res->ts_v[i]->AMPLITUDE,j,fura_c);
      gsl_vector_set(res->ts_v[i]->TIME,j,gsl_vector_get(time,j));
    }
    if (i == 0)
      max_adu = gsl_vector_max(res->ts_v[i]->AMPLITUDE);
    gsl_vector_scale(res->ts_v[i]->AMPLITUDE,F_p/max_adu);
  }
  return res;
}
