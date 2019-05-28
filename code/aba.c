/** \file aba.c
    \brief Function definitions for [aba](@ref aba) structures. 
 */
#include "abaa.h"
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
