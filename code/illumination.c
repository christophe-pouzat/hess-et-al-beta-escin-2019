/** \file illumination.c
    \brief Function definitions for [illumination](@ref illumination) structures.
 */
#include "abaa.h"
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
