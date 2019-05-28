/** @file ccd_test.c
 *  @brief test functions reading CCD data from an HDF5 file
 *         and printing them to the stdout 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_vector.h>
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
#define FNAME "data_paper/data_beta_escin/DA_121219_E1.h5"
int main()
{
  // Open FILE
  char fname[] = FNAME;
  hid_t fid = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  ccd ccd_st = ccd_read_from_file(fid);
  // Close file
  H5Fclose (fid);
  ccd_printf(&ccd_st);
  return 0;
}
