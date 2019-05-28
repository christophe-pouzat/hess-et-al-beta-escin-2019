/** @file illumination_test.c
 *  @brief test functions reading ILLUMINATION data from an HDF5 file
 *         and printing them to the stdout 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_vector.h>
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
#define FNAME "data_paper/data_beta_escin/DA_121219_E1.h5"
int main()
{
  // Open FILE
  char fname[] = FNAME;
  hid_t fid = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  illumination illumination_st = illumination_read_from_file(fid);
  // Close file
  H5Fclose (fid);
  illumination_printf(&illumination_st);
  return 0;
}
