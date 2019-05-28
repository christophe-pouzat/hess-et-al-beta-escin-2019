/** @file adu_test.c
 *  @brief Test program for adu structure and related
 *         functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
int main()
{
  // Allocate an adu with 2 observations
  adu * adu_ptr = adu_alloc(2);
  printf("Allocated adu_ptr with two elements.\n");
  printf("Setting values of ADU340 field.\n");
  adu_set(adu_ptr,ADU340,0,1.0);
  adu_set(adu_ptr,ADU340,1,2.0);
  printf("Setting values of ADU340B field.\n");
  adu_set(adu_ptr,ADU340B,0,3.0);
  adu_set(adu_ptr,ADU340B,1,4.0);
  printf("Setting values of ADU360 field.\n");
  adu_set(adu_ptr,ADU360,0,5.0);
  adu_set(adu_ptr,ADU360,1,6.0);
  printf("Setting values of ADU360B field.\n");
  adu_set(adu_ptr,ADU360B,0,7.0);
  adu_set(adu_ptr,ADU360B,1,8.0);
  printf("Setting values of ADU380 field.\n");
  adu_set(adu_ptr,ADU380,0,9.0);
  adu_set(adu_ptr,ADU380,1,10.0);
  printf("Setting values of ADU380B field.\n");
  adu_set(adu_ptr,ADU380B,0,11.0);
  adu_set(adu_ptr,ADU380B,1,12.0);
  printf("Setting values of TIME field.\n");
  adu_set(adu_ptr,TIME,0,1.5);
  adu_set(adu_ptr,TIME,1,2.5);
  printf("The content of the structure is.\n");
  adu_printf(adu_ptr);
  // free allocated adu
  adu_free(adu_ptr);
  printf("Freed adu_ptr. Don't forget running valgrind!\n");
  return 0;
}
