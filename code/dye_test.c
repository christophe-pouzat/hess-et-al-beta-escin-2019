/** @file dye_test.c
 *  @brief test functions reading DYE data from an HDF5 file
 *         and printing them to the stdout 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_vector.h>
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
#define FNAME "data_paper/data_beta_escin/DA_121219_E1.h5"
int main()
{
  // Open FILE
  char fname[] = FNAME;
  hid_t fid = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  dye dye_st = dye_read_from_file(fid);
  // Close file
  H5Fclose (fid);
  dye_printf(&dye_st);
  return 0;
}
