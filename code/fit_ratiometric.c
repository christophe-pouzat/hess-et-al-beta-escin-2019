/** \file fit_ratiometric.c
    \brief Fits mono-exponential model to ratiometric estimator
 */
#include "abaa.h"
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s -i --input=string -o --output=string [-s --stim=integer] ...\n"
    "          ... [-m --maxiter=integer] [-b --baseline_length=integer] ...\n"
    "          ... [-f --start_fit=real] [-g --graphic]\n\n"
    "  -i --input <character string>: data file name (e.g. data_paper/data_beta_escin/DA_121219_E1.h5)\n"
    "  -o --output <character string>: output file name prefix (e.g. 'DA_121219_E1_s1_ratio')\n"
    "  -s --stim <positive integer>: the stimulation to fit (default 1)\n"
    "  -m --max_iter <positive integer>: maximal number of iterations performed by\n"
    "       the nonlinear least-squares solver (default 50)\n"
    "  -b --baseline_length <positive integer>: baseline length in samples (default 15)\n"
    "  -f --start_fit <positive real>: where decay fit starts, if > 1 interpreted as the\n"
    "       number of samples after the peak, if 0<f<1 interpreted as the remaining\n"
    "       fraction of the jump and fitting starts when the transient reaches that level\n"
    "  -g --graphics: if set, a gnuplot script file named output.gp is generated\n\n"
    " The program opens 'input' file, computes the ratiometric estimator from the raw data\n"
    " and fits the following model to it: Ca = baseline + delta * exp (-(t-t0)/tau)\n"
    " where t0 is the time on the decay phase at which the fit starts (set by parameter 'f').\n"
    " A constant is also fitted to the baseline region: Ca = baseline.\n"
    " Standard errors for the ratiometric estimator are obtained by Monte-Carlo simulation.\n"
    " While the program performs the least-squares optimization relevant information get printed\n"
    " to the stderr. The ratiometric estimator and its standard error are printed to a text\n"
    " file named 'output_CaRatio'. The fitted parameters, fitted values and residuals are\n"
    " printed to a text file named 'output_RatioFit'. If 'g' is selected, a gnuplot script\n"
    " file whose name is the 'output' name with the '.gp' suffix is generated. A graphical display\n"
    " of the result can then be obtained with: gnuplot -persist 'output.gp'\n\n";

  char *filename;
  char *output;
  size_t stim=1, maxiter=50, baseline_length=15;
  double start_fit=0.5;
  int g_script=0;
  {int opt;
    static struct option long_options[] = {
      {"input",required_argument,NULL,'i'},
      {"output",required_argument,NULL,'o'},
      {"stim",optional_argument,NULL,'s'},
      {"graphic",optional_argument,NULL,'g'},
      {"maxiter",optional_argument,NULL,'m'},
      {"baseline_length",optional_argument,NULL,'b'},
      {"start_fit",optional_argument,NULL,'f'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hgi:o:s:m:b:f:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 'i':
      {
        filename = optarg;
      }
      break;
      case 'o':
      {
        output = optarg;
      }
      break;
      case 'g':
      {
        g_script=1;
      }
      break;
      case 's':
      {
        stim = (size_t) atoi(optarg);
      }
      break;
      case 'm':
      {
        maxiter = (size_t) atoi(optarg);
      }
      break;
      case 'b':
      {
        baseline_length = (size_t) atoi(optarg);
      }
      break;
      case 'f':
      {
        start_fit = (double) atof(optarg);
        if (start_fit<=0) {
          fprintf(stderr,"start_fit should be > 0.\n");
          return -1;
        }
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
  ratio * pratio = ratio_est(paba->data->adu_v[stim],&paba->dye,&paba->light,&paba->ccd,1000);
  aba_free(paba);
  char out[512];
  strcpy(out,output);
  strcat(out,"_CaRatio");
  FILE *fp = fopen(out,"w");
  fprintf(fp,"# Data from %s stim %d\n",filename,(int) stim);
  fprintf(fp,"# Ratiometric estimator\n");
  ratio_fprintf(fp,pratio);
  fclose(fp);
  size_t start_pos;
  if (start_fit >= 1)
    start_pos = (size_t) start_fit;
  else
    start_pos = ratio_find_fit_start(pratio,start_fit,baseline_length);
  mono_exp_fit_res res = ratio_fit(pratio,baseline_length,start_pos,maxiter);
  strcpy(out,output);
  strcat(out,"_RatioFit");
  fp = fopen(out,"w");
  fprintf(fp,"# Data from %s stim %d\n",filename,(int) stim);
  fprintf(fp,"# Ratiometric estimator mono-exponential fit\n");
  mono_exp_fit_res_fprintf(fp,&res,pratio);
  fclose(fp);
  ratio_free(pratio);
  if (g_script) {
    strcpy(out,output);
    strcat(out,".gp");
    fp = fopen(out,"w");
    char RatioFit[512];
    strcpy(RatioFit,output);
    strcat(RatioFit,"_RatioFit");
    char CaRatio[512];
    strcpy(CaRatio,output);
    strcat(CaRatio,"_CaRatio");
    char format[] = "%g";
    fprintf(fp,
    	"unset key\n"
    	"set grid\n"
    	"set multiplot title 'Data from %s stim %d' noenhanced layout 2,1 margins 0.1,0.9,0.1,0.95 spacing 0,0\n"
    	"set ylabel 'Normalized residuals'\n"
    	"set xtics format ''\n"
    	"plot '%s' using 1:4 with lines linecolor rgb 'red' linewidth 2\n"
    	"set format x '%s'\n"
    	"set xlabel 'Time (s)'\n"
    	"set ylabel 'Estimated [Ca2+]'\n"
    	"plot '%s' using 1:2:($3*1.96) with yerrorlines	\\\n"
    	"     linecolor rgb 'black' linewidth 1,		\\\n"
    	"     '%s' using 1:3 with lines linecolor rgb 'red'\\\n"
    	"     linewidth 2\n"
    	"unset multiplot",
    	filename, (int) stim, RatioFit, format, CaRatio, RatioFit);
    fclose(fp);
  }
  return 0;
}
