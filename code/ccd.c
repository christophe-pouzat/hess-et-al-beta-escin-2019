/** \file ccd.c
    \brif Function definitions for [ccd](@ref ccd) structures.
 */
#include "abaa.h"
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
