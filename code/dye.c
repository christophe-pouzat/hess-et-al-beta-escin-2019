/** \file dye.c
    \brief Function definitions for [dye](@ref dye) structures
 */
#include "abaa.h"
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
