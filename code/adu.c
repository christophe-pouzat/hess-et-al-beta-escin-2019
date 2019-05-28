/** \file adu.c
    \brief Definitions of functions related to [adu](@ref adu) structures
 */
#include "abaa.h"
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
