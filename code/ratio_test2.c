/** \file ratio_test2.c
    \brief Test program for the ratiometric estimator functions and
           structures of the abaa library.
 */
#include "abaa.h"
#define FNAME "data_paper/data_beta_escin/DA_121219_E1.h5"
int main()
{
  // Open FILE
  char fname[] = FNAME;
  hid_t fid = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  aba * paba = aba_read_from_file(fid);
  // Close file
  H5Fclose (fid);
  ratio * pratio = ratio_est(paba->data->adu_v[3],&paba->dye,&paba->light,&paba->ccd,1000);
  aba_free(paba);
  ratio_fprintf(stdout,pratio);
  printf("\n\n");
  size_t start_pos = ratio_find_fit_start(pratio,0.5,15);
  mono_exp_fit_res res = ratio_fit(pratio,15,start_pos,25);
  mono_exp_fit_res_fprintf(stdout,&res,pratio);
  ratio_free(pratio);
  return 0;
}
