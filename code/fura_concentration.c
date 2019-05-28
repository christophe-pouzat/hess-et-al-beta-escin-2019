/** \file fura_concentration.c
    \brief Program computing the [Fura] and returning it to the stdout
 */
#include "abaa.h"
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s -i --input=string\n"
    "  -i --input <character string>: data file name (e.g. data_paper/data_beta_escin/DA_121219_E1.h5)\n\n"
    " The program opens 'input' file, estimates the [Fura] within the cell\n"
    " and prints it to the 'stdout'. The loading curve measurements appear first\n"
    " followed by the individual stimalation separated by two blank lines.\n"
    " Results are printed on two columns: Time and [Fura].\n\n"
    " IT IS ASSUMED THAT WHEN THE MAXIMAL BACKGROUND SUBTRACTED FLUORESCENCE AT\n"
    " 360 nm IS OBSERVED, THE FURA CONCENTRATION IN THE CELL AND IN THE PIPETTE\n"
    " ARE IDENTICAL (NO EXTRAPOLATION BASED ON A FIT IS PERFORMED).\n\n";

  char *filename;
  {int opt;
    static struct option long_options[] = {
      {"input",required_argument,NULL,'i'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hi:s:f:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 'i':
      {
        filename = optarg;
      }
      break;
      case 'h': printf(usage,argv[0]);
        return -1;
      default : fprintf(stderr,usage,argv[0]);
        return -1;
      }
    }
  }

  hid_t fid = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  aba * paba = aba_read_from_file(fid);
  // Close file
  H5Fclose (fid);
  ts_vector *pfura = fura_est(paba);
  aba_free(paba);
  ts_vector_fprintf(stdout,pfura);
  ts_vector_free(pfura);
  return 0;
}
