/** \file aba_ratio.c
    \brief Added buffer approach using the ratiometric estimator

    This program does "everything at once":
    1. Reads the data.
    2. Get the loading curve.
    3. Fits a mono-exponential decay to the ratiometric estimators of the user specified stimulations.
    4. Fits a weighted linear regression to the decay time constant vs \f$\kappa_{F}\f$ and get the two key parameters: \f$\gamma/v\f$ and \f$\kappa_{S}\f$.
    
    The last fit, time constant vs \f$\kappa_{F}\f$, is performed in 3 different ways using respectively the mean \f$\kappa_{F}\f$
    value during the fitted part of the transient, the minimal or the maximal value. 
    A text file (or several ones, like in step 3) is/are generated at each step.
 */
#include "abaa.h"

/** \def NREP 10000
    \brief Number of replicates in MC/Bootstrap simulations
 */
#define NREP 10000



int wls_fit(const gsl_vector * kappa,
	    mono_exp_fit_res * fit_res,
	    size_t n_stim,
	    char * filename,
	    char * output,
	    gsl_rng * rng,
	    int what,
	    size_t *stim);

int robust_fit(const gsl_vector * kappa,
	       mono_exp_fit_res * fit_res,
	       size_t n_stim,
	       char * filename,
	       char * output,
	       gsl_rng * rng,
	       int what,
	       size_t *stim);

int main(int argc, char *argv[])
{
  /*
    This code block defines a character strings containing the
    help of the program.
   */
  static char usage[] = \
    "usage: %s -i --input=string [-o --output=string] ...\n"
    "          ... [-s --stim=integer,integer,integer] ...\n"
    "          ... [-m --maxiter=integer] [-b --baseline_length=integer] ...\n"
    "          ... [-f --start_fit=real] [-r --robust] [-g --graphic]\n\n"
    "  -i --input <character string>: data file name (e.g. data_paper/data_beta_escin/DA_121219_E1.h5)\n"
    "  -o --output <character string>: output file name prefix (e.g. 'DA_121219_E1_s1_ratio');\n"
    "       if not specificied, the '.h5' is stripped from 'input' and '_aba' is added\n"
    "  -s --stim <positive integer>: coma separated list of stimulations to fit\n"
    "       (default all stimulation considered).\n"
    "  -m --max_iter <positive integer>: maximal number of iterations performed by\n"
    "       the nonlinear least-squares solver (default 50).\n"
    "  -b --baseline_length <positive integer>: baseline length in samples (default 15).\n"
    "  -f --start_fit <positive real>: where decay fit starts, if > 1 interpreted as the\n"
    "       number of samples after the peak, if 0<f<1 interpreted as the remaining\n"
    "       fraction of the jump and fitting starts when the transient reaches that level.\n"
    "       Default value: 0.5.\n"
    "  -r --robust: if set, a robust linear regression of tau vs kappa_Fura is performed\n"
    "       using a bi-square weight function.\n"
    "  -g --graphics: if set, gnuplot script files are generated.\n\n"
    " The program opens 'input' file, get the loading curve and prints the result to a file\n"
    " named output_loading_curve; then for every stimulation specified in the 'stim' list \n"
    " it computes the ratiometric estimator from the raw data\n"
    " and fits the following model to it: Ca = baseline + delta * exp (-(t-t0)/tau)\n"
    " where t0 is the time on the decay phase at which the fit starts (set by parameter 'f').\n"
    " A constant is also fitted to the baseline region: Ca = baseline.\n"
    " Standard errors for the ratiometric estimator are obtained by Monte-Carlo simulation.\n"
    " While the program performs the least-squares optimization relevant information get printed\n"
    " to the stderr. The ratiometric estimator and its standard error are printed to a text\n"
    " file named 'output_CaRatio_sX' (where 'X' stands for the stimulation number. The fitted\n"
    " parameters, fitted values and residuals are printed to a text file named 'output_RatioFit_sX'.\n"
    " If 'g' is selected, a gnuplot script file whose name is the 'output_YY_sX.gp' where 'YY' is\n"
    " either 'loading_curve' or 'RatioFit' 'is generated. A graphical display\n"
    " of the result can then be obtained with: gnuplot -persist 'output_YY_sX.gp'\n\n"
    " Once each specified stimulation has been fitted, a regression of tau vs kappa is performed\n"
    " using three different values for kappa: the mean value during the fitted part of the transient,\n"
    " the minimal and the maximal values. Information and fit results get printed to the 'stderr' as\n"
    " well as to files named 'output_aba_tau_vs_mean/min/max_kappa', if option '-g' is set, gnuplot\n"
    " script files are also generated with names 'output_aba_tau_vs_mean/min/max_kappa.gp'\n\n"
    " When a robust linear regression is requested (by using optional argument '-r' or '--robust'\n"
    " when calling the program), a bi-square weight function is used (see:\n"
    " https://www.gnu.org/software/gsl/doc/html/lls.html#robust-linear-regression). Since\n"
    " the robust methods implemented in the GSL do not allow for the use of the standard errors on\n"
    " the dependent variable, we multiply the observations in order to have approximately the\n"
    " right standard error ratio: that is, if an observation has a standard error of 0.13 and\n"
    " another one has 0.091, we will include two identical copies of the second since (0.13/0.091)^2\n"
    " is approximately 2. This use of multiple copies is reported to the user.\n\n";

  /*
     This code block reads the program parameters. It assigns (after allocating memory if necessary):
     - filename: pointer to a character string with the name of the data file.
     - output: pointer to a character string with the name of the output file name.
     - stim: an array with the numbers of the stimulations to use.
     - maxiter, baseline_length: maximal number of iterations and baseline length.
     - start_fit: the number used to set from where the decay gets fitted.
     - g_script: indicates if a gnuplot script should be generated.
    */  
  char *filename;
  char output[512];
  size_t maxiter=50, baseline_length=15;
  double start_fit=0.5;
  int g_script=0;
  size_t n_stim = 0; // initialize the stimulation counter
  size_t *stim;
  int do_robust = 0; 
  int out_unset=1; // Indicator of output specification  
  {int opt;
    static struct option long_options[] = {
      {"input",required_argument,NULL,'i'},
      {"output",optional_argument,NULL,'o'},
      {"stim",optional_argument,NULL,'s'},
      {"graphic",optional_argument,NULL,'g'},
      {"maxiter",optional_argument,NULL,'m'},
      {"baseline_length",optional_argument,NULL,'b'},
      {"start_fit",optional_argument,NULL,'f'},
      {"robust",optional_argument,NULL,'r'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hgri:o:s:m:b:f:",
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
        out_unset = 0;
        strcpy(output,optarg);
      }
      break;
      case 'g':
      {
        g_script=1;
      }
      break;
      case 'r':
      {
        do_robust=1;
      }
      break;
      case 's':
      {
        char *start = strdup(optarg); // duplicate optarg content
        char *running;
        running = start;
        char *token  = strsep(&running, ",");
        while (token != NULL) {
  	token = strsep (&running, ","); // split optarg at each ","
  	n_stim++;
        }
        free(start);
        // The number of stimulation is now known
        // Allocate memory for the vector of stim indexes
        stim=malloc(n_stim*sizeof(size_t));
        start = strdup(optarg); // duplicate optarg content again
        running = start;
        // Get the index of each stimulation
        for (size_t i=0; i<n_stim; i++) {
  	token = strsep (&running, ",");
  	stim[i] = (size_t) atoi(token);
        }
        free(start);
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
  // Set output prefix name if not given
  if (out_unset) {
    char *start = strdup(filename); // duplicate filename content
    char *p = strtok(start,"."); // stop at ".h5"
    strcpy(output,p);
    strcat(output,"_aba");
    free(start);
  }

  /*
     This code block opens HDF5 file `filename`. 
     `paba` is a pointer to an aba structure that gets
     allocated and initialized.
     Do not forget to free the memory pointed to by `paba` before
     program exit. 
  */
  hid_t fid = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  aba * paba = aba_read_from_file(fid);
  // Close file
  H5Fclose (fid);
  // Take care of the default for stim
  if (n_stim == 0) {
    n_stim = paba->data->nelt-1;
    stim = malloc(n_stim*sizeof(size_t));
    for (size_t s_idx=0; s_idx<n_stim; s_idx++)
      stim[s_idx]=s_idx+1;
  }

  /*
     This code block allocates and initializes a ts_vector
     structure pointed to by `pfura`.
     A char array, `out`, with 512 elements is declared.
     A file pointer, `fp`, is declared.
     The loading curve is printed to a file named `output_loading_curve`.
     Do not forget to free `pfura` before program exit.
  */
  ts_vector *pfura = fura_est(paba);
  char out[512];
  strcpy(out,output);
  strcat(out,"_loading_curve");
  FILE *fp = fopen(out,"w");
  fprintf(fp,"# Cellular [Fura] for dataset %s\n\n",filename);
  ts_vector_fprintf(fp,pfura);
  fclose(fp);

  if (g_script) {
    /*
      This code block generated a gnuplot script
      making the loading curve figure.
      A char array, `file4gp`, with 512 elements is declared.
     */
    char file4gp[512];
    strcpy(file4gp,output);
    strcat(file4gp,"_loading_curve");
    strcpy(out,output);
    strcat(out,"_loading_curve.gp");
    fp = fopen(out,"w");
    fprintf(fp,
    	"unset key\n"
    	"set grid\n"
    	"set ylabel '[Fura] (µM)'\n"
    	"set xlabel 'Time (s)'\n"
    	"plot '%s' index 0 using 1:2 linecolor rgb 'black',\\\n",
    	file4gp);
    for (size_t i=0; i<n_stim; i++)
      fprintf(fp,
    	  "'' index %d using 1:2 with lines linecolor rgb 'red',\\\n",
    	  (int) i+1);
    fclose(fp);
  }

  /*
  The next code blocks get the ratiometric estimator
  and to the mono-exponential fit for each selected
  stimulation. 
  */
  mono_exp_fit_res fit_res[n_stim];
  for (size_t s_idx=0; s_idx<n_stim; s_idx++) {
    fprintf(stderr,
    	"**********************************\n"
    	"* Doing now stimulation %d\n"
    	"**********************************\n",
    	(int) stim[s_idx]);
    ratio * pratio = ratio_est(paba->data->adu_v[stim[s_idx]],
    			   &paba->dye,&paba->light,
    			   &paba->ccd,NREP);
    strcpy(out,output);
    strcat(out,"_CaRatio_s");
    char stim_nb[99];
    sprintf(stim_nb,"%d", (int) stim[s_idx]);
    strcat(out,stim_nb);
    FILE *fp = fopen(out,"w");
    fprintf(fp,"# Data from %s stim %d\n",filename,(int) stim[s_idx]);
    fprintf(fp,"# Ratiometric estimator\n");
    ratio_fprintf(fp,pratio);
    fclose(fp);
    size_t start_pos;
    if (start_fit >= 1)
      start_pos = (size_t) start_fit;
    else
      start_pos = ratio_find_fit_start(pratio,start_fit,baseline_length);
    fit_res[s_idx] = ratio_fit(pratio,baseline_length,start_pos,maxiter);
    strcpy(out,output);
    strcat(out,"_RatioFit_s");
    strcat(out,stim_nb);
    fp = fopen(out,"w");
    fprintf(fp,"# Data from %s stim %d\n",filename,(int) stim[s_idx]);
    fprintf(fp,"# Ratiometric estimator mono-exponential fit\n");
    mono_exp_fit_res_fprintf(fp,&fit_res[s_idx],pratio);
    fclose(fp);
    ratio_free(pratio);
    fprintf(stderr,
    	"**********************************\n"
    	"* Stimulation %d done\n"
    	"**********************************\n\n",
    	(int) stim[s_idx]);
    if (g_script) {
      strcpy(out,output);
      strcat(out,"_RatioFit_s");
      strcat(out,stim_nb);
      strcat(out,".gp");
      fp = fopen(out,"w");
      char RatioFit[512];
      strcpy(RatioFit,output);
      strcat(RatioFit,"_RatioFit_s");
      strcat(RatioFit,stim_nb);
      char CaRatio[512];
      strcpy(CaRatio,output);
      strcat(CaRatio,"_CaRatio_s");
      strcat(CaRatio,stim_nb);
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
              "plot '%s' using 1:2:($3*1.96) with yerrorlines \\\n"
              "     linecolor rgb 'black' linewidth 1,                \\\n"
              "     '%s' using 1:3 with lines linecolor rgb 'red'\\\n"
              "     linewidth 2\n"
              "unset multiplot",
              filename, (int) stim[s_idx], RatioFit, format, CaRatio, RatioFit);
      fclose(fp);
    }
  }

  /*
  The next block gets the mean, min and max kappa_Fura 
  values at each stimulation. The result is stored
  in three gsl_vector called kappa_mu, kappa_inf and kappa_sup.
  These vectors must be freed before exiting the program.
  */
  gsl_vector * kappa_mu = gsl_vector_alloc(n_stim);
  gsl_vector * kappa_inf = gsl_vector_alloc(n_stim);
  gsl_vector * kappa_sup = gsl_vector_alloc(n_stim);
  for (size_t s_idx=0; s_idx<n_stim; s_idx++) {
    double K_d=paba->dye.K_d_hat;
    double Ca = (fit_res[s_idx]).baseline;
    double res=K_d/gsl_pow_2(K_d+Ca);
    double res_inf,res_sup;
    gsl_vector * fura_conc=pfura->ts_v[stim[s_idx]]->AMPLITUDE;
    gsl_vector_const_view Fura=gsl_vector_const_subvector(fura_conc,
  							(fit_res[s_idx]).fit_start,
  							fura_conc->size-(fit_res[s_idx]).fit_start);
    // We take as [Fura] the mean [Fura] over the fitted decay period
    res *= gsl_stats_mean(Fura.vector.data,1,Fura.vector.size);
    gsl_stats_minmax(&res_inf,&res_sup,Fura.vector.data,1,Fura.vector.size);
    gsl_vector_set(kappa_mu, s_idx, res);
    gsl_vector_set(kappa_inf, s_idx, res_inf*K_d/gsl_pow_2(K_d+Ca));
    gsl_vector_set(kappa_sup, s_idx, res_sup*K_d/gsl_pow_2(K_d+Ca));
  }				
  
  /*
  The next code block does 3 linear regressions
  of tau on mean / min / max kappa_Fura
  */
  // Set up RNG for MC
  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng_env_setup();
  gsl_rng * rng = gsl_rng_alloc (T);
  
  if (do_robust) {
    // Fit and print tau vs mean kappa value
    robust_fit(kappa_mu,fit_res,n_stim,filename,output,rng,0,stim);
  
    // Fit and print tau vs min kappa value
    robust_fit(kappa_inf,fit_res,n_stim,filename,output,rng,1,stim);
  
    // Fit and print tau vs max kappa value
    robust_fit(kappa_sup,fit_res,n_stim,filename,output,rng,2,stim);
  } else {
    // Fit and print tau vs mean kappa value
    wls_fit(kappa_mu,fit_res,n_stim,filename,output,rng,0,stim);
  
    // Fit and print tau vs min kappa value
    wls_fit(kappa_inf,fit_res,n_stim,filename,output,rng,1,stim);
  
    // Fit and print tau vs max kappa value
    wls_fit(kappa_sup,fit_res,n_stim,filename,output,rng,2,stim);  
  }
  // Free RNG
  gsl_rng_free (rng);
  

  
  if (g_script) {
    char file4gp[512];
    //define a macro generating the gnuplot script
    #define mk_gp_script(p_name)\
      strcpy(file4gp,output);\
      strcat(file4gp,"_tau_vs_" #p_name "_kappa");\
      strcpy(out,output);\
      strcat(out,"_tau_vs_" #p_name "_kappa.gp");\
      fp = fopen(out,"w");\
      fprintf(fp,\
        "unset key\n"\
        "set grid\n"\
        "set ylabel 'τ (s)'\n"\
        "set xlabel '" #p_name " κ_Fura' noenhanced\n"\
        "set title 'Data from %s' noenhanced\n"\
        "set xrange [-1000<*:*]\n"\
        "plot '%s' index 0 using 1:2:(1.96*$3) linecolor rgb 'red' with yerrorbars,\\\n",\
    	  filename,file4gp);\
      fprintf(fp,\
    	  "'' index 1 using 1:2 with lines linecolor rgb 'black',\\\n");\
      fprintf(fp,\
    	  "'' index 1 using 1:3 with lines linecolor rgb 'blue' lt 'dotted',\\\n");\
      fprintf(fp,\
    	  "'' index 1 using 1:4 with lines linecolor rgb 'blue' lt 'dotted'");\
      fclose(fp);\
    
    #define mk_gp_script_robust(p_name)\
      strcpy(file4gp,output);\
      strcat(file4gp,"_tau_vs_" #p_name "_kappa_robust");\
      strcpy(out,output);\
      strcat(out,"_tau_vs_" #p_name "_kappa_robust.gp");\
      fp = fopen(out,"w");\
      fprintf(fp,\
        "unset key\n"\
        "set grid\n"\
        "set ylabel 'τ (s)'\n"\
        "set xlabel '" #p_name " κ_Fura' noenhanced\n"\
        "set title 'Data from %s' noenhanced\n"\
        "plot '%s' index 0 using 1:2:(1.96*$3) linecolor rgb 'red' with yerrorbars,\\\n",\
    	  filename,file4gp);\
      fprintf(fp,\
    	  "'' index 1 using 1:2 with lines linecolor rgb 'black',\\\n");\
      fprintf(fp,\
    	  "'' index 1 using 1:3 with lines linecolor rgb 'black' lt 'dotted',\\\n");\
      fprintf(fp,\
    	  "'' index 1 using 1:4 with lines linecolor rgb 'black' lt 'dotted',\\\n");\
      fprintf(fp,\
    	  "'' index 1 using 1:5 with lines linecolor rgb 'blue',\\\n");\
      fprintf(fp,\
    	  "'' index 1 using 1:6 with lines linecolor rgb 'blue' lt 'dotted',\\\n");\
      fprintf(fp,\
    	  "'' index 1 using 1:7 with lines linecolor rgb 'blue' lt 'dotted'");\
      fclose(fp);\
    
    if (do_robust) {
    mk_gp_script_robust(mean)
    
    mk_gp_script_robust(min)
    
    mk_gp_script_robust(max)
    } else {
    mk_gp_script(mean)
    
    mk_gp_script(min)
    
    mk_gp_script(max)
    }
    
  }

  gsl_vector_free(kappa_mu);
  gsl_vector_free(kappa_inf);
  gsl_vector_free(kappa_sup);
  aba_free(paba);
  ts_vector_free(pfura);
  free(stim); 
  return 0;
}

#define BETA0(c) gsl_vector_get((c),0)

#define BETA0_SE(cov) sqrt(gsl_matrix_get((cov),0,0))

#define BETA0_VAR(cov) gsl_matrix_get((cov),0,0)

#define BETA1(c) gsl_vector_get((c),1)

#define BETA1_SE(cov) sqrt(gsl_matrix_get((cov),1,1))

#define BETA1_VAR(cov) gsl_matrix_get((cov),1,1)

#define gamma_over_v_hat(c) 1.0/BETA1((c))

#define gamma_over_v_se(c,cov) BETA1_SE((cov))/	\
  gsl_pow_2(BETA1((c)))

#define kappa_s(c) BETA0((c))/BETA1((c))-1.0

#define kappa_s_se(c,cov) sqrt(BETA0_VAR((cov))/gsl_pow_2(BETA1((c)))+\
			       BETA1_VAR((cov))*gsl_pow_2(BETA0((c)))/gsl_pow_4(BETA1((c))))

/// @brief Performs a weighted linear regression and prints results
///        to `stderr` and to a file whose prefix is built from
///        `output`.
///
/// @returns 0 if everything is fine
int wls_fit(const gsl_vector * kappa, ///< [in] vector of kappa values.  
	    mono_exp_fit_res * fit_res, ///< [in] vector of mono_exp_fit_res.
	    size_t n_stim, ///< [in] the number of stimulation.
	    char * filename, ///< Name of file containing the data
	    char * output, ///< Prefix of file name where results are written
	    gsl_rng * rng, ///< Pointer to a gsl_rng
	    int what, ///< What is processed? 0 = mean, 1 = min, 2 = max
	    size_t *stim ///< [in] Vector of indexes of used stimulation
  ) 
{
  fprintf(stderr,"******************************************\n");
  if (what==0)
    fprintf(stderr,"* Doing tau vs mean kappa_Fura regression *\n");
  if (what==1)
    fprintf(stderr,"* Doing tau vs min kappa_Fura regression *\n");
  if (what==2)
    fprintf(stderr,"* Doing tau vs max kappa_Fura regression *\n");
  fprintf(stderr,"******************************************\n");

  gsl_vector * c = gsl_vector_alloc(2);
  gsl_matrix * cov = gsl_matrix_alloc(2,2);
  double chisq,tss;
  // Prepare data for fit
  gsl_vector * tau = gsl_vector_alloc(n_stim); // The dependent variable
  gsl_vector * w = gsl_vector_alloc(n_stim); // The weight of the above
  for (size_t s_idx=0; s_idx<n_stim; s_idx++) {
    gsl_vector_set(tau, s_idx, (fit_res[s_idx]).tau);
    gsl_vector_set(w, s_idx, 1.0/gsl_pow_2((fit_res[s_idx]).tau_se));
  }
  // Get total sum of squares (TSS)
  tss = gsl_stats_wtss(w->data,1,tau->data,1,tau->size);
  // Allocate design matrix
  gsl_matrix * X = gsl_matrix_alloc(n_stim,2);
  // Fill design matrix
  for (size_t s_idx=0; s_idx<n_stim; s_idx++) {
    gsl_matrix_set(X,s_idx,0,1.0);
    gsl_matrix_set(X,s_idx,1,gsl_vector_get(kappa,s_idx));
  }
  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n_stim, 2);
  gsl_multifit_wlinear (X, w, tau, c, cov, &chisq, work);
  gsl_multifit_linear_free (work);
  

  double kappa_S_MC[NREP];
  double beta0 = BETA0(c);
  double beta1 = BETA1(c);
  double cov00 = BETA0_VAR(cov);
  double cov11 = BETA1_VAR(cov);
  double beta0_se = BETA0_SE(cov);
  double kappa_S = kappa_s(c);
  double kappa_S_se = kappa_s_se(c,cov);
  double se11=BETA1_SE(cov);
  double rho=gsl_matrix_get(cov,0,1)/beta0_se/se11;
  for (size_t rep_idx=0; rep_idx<NREP; rep_idx++) {
    double num,denom;
    gsl_ran_bivariate_gaussian(rng,beta0_se,se11,rho,&num,&denom);
    num+=beta0;
    denom+=beta1;
    //double num=(beta0+gsl_ran_gaussian_ziggurat(rng,beta0_se));
    //double denom=(beta1+gsl_ran_gaussian_ziggurat(rng,se11));
    kappa_S_MC[rep_idx] = num/denom-1.0;
  }
  gsl_sort(kappa_S_MC,1,NREP);
  double kappa_S_l95=kappa_S_MC[(size_t)(0.025*NREP-1)];
  double kappa_S_u95=kappa_S_MC[(size_t)(0.975*NREP-1)];
  double kappa_S_l99=kappa_S_MC[(size_t)(0.005*NREP-1)];
  double kappa_S_u99=kappa_S_MC[(size_t)(0.995*NREP-1)];

  fprintf(stderr,"Best fit: tau = %g + %g kappa_Fura\n",beta0 , beta1 );
  fprintf(stderr,"Covariance matrix:\n");
  fprintf(stderr,"[ %+.5e, %+.5e  \n", cov00, gsl_matrix_get(cov,0,1));
  fprintf(stderr,"  %+.5e, %+.5e  ]\n", gsl_matrix_get(cov,1,0), cov11);
  fprintf(stderr,"Total sum of squares (TSS) = %g\n", tss);
  fprintf(stderr,"chisq (Residual sum of squares, RSS) = %g\n", chisq);
  double p_larger = gsl_cdf_chisq_Q(chisq,n_stim-2);
  fprintf(stderr,"Probability of observing a larger of equal RSS per DOF under the null hypothesis: %g\n",p_larger);
  fprintf(stderr,"R squared (1-RSS/TSS) = %g\n", 1 - chisq / tss);
  fprintf(stderr,"Estimated gamma/v with standard error: %g +/- %g\n",
	 gamma_over_v_hat(c), gamma_over_v_se(c,cov));
  fprintf(stderr,"Estimated kappa_S with standard error (using error propagation): %g +/- %g\n",
	  kappa_S, kappa_S_se);
  fprintf(stderr,"kappa_S confidence intervals based on parametric bootstrap\n");
  fprintf(stderr,"0.95 CI for kappa_S: [%g,%g]\n",kappa_S_l95,kappa_S_u95);
  fprintf(stderr,"0.99 CI for kappa_S: [%g,%g]\n",kappa_S_l99,kappa_S_u99);
  fprintf(stderr,"******************************************\n");
  if (what==0)
    fprintf(stderr,"* tau vs mean kappa_Fura regression done *\n");
  if (what==1)
    fprintf(stderr,"* tau vs min kappa_Fura regression done  *\n");
  if (what==2)
    fprintf(stderr,"* tau vs max kappa_Fura regression done  *\n");
  fprintf(stderr,"******************************************\n");

  // Print to file
  char out[512];
  strcpy(out,output);
  if (what==0)
    strcat(out,"_tau_vs_mean_kappa");
  if (what==1)
    strcat(out,"_tau_vs_min_kappa");
  if (what==2)
    strcat(out,"_tau_vs_max_kappa");
  FILE *fp = fopen(out,"w");
  if (what==0)
    fprintf(fp,"# τ vs mean κ linear regression for data set %s\n",filename);
  if (what==1)
    fprintf(fp,"# τ vs min κ linear regression for data set %s\n",filename);
  if (what==2)
    fprintf(fp,"# τ vs max κ linear regression for data set %s\n",filename);
  fprintf(fp,"# Using stimulation: %d",(int) stim[0]);
  for (size_t s_idx=1; s_idx<n_stim; s_idx++)
    fprintf(fp,", %d", (int) stim[s_idx]);
  fprintf(fp,"\n");
  fprintf(fp,"# Best fit: tau = %g + %g kappa_Fura\n",beta0 , beta1 );
  fprintf(fp,"# Covariance matrix:\n");
  fprintf(fp,"# [ %+.5e, %+.5e  \n", cov00, gsl_matrix_get(cov,0,1));
  fprintf(fp,"#   %+.5e, %+.5e  ]\n", gsl_matrix_get(cov,1,0), cov11);
  fprintf(fp,"# Total sum of squares (TSS) = %g\n", tss);
  fprintf(fp,"# chisq (Residual sum of squares, RSS) = %g\n", chisq);
  fprintf(fp,"# Probability of observing a larger of equal RSS per DOF under the null hypothesis: %g\n",p_larger);
  fprintf(fp,"# R squared (1-RSS/TSS) = %g\n", 1 - chisq / tss);
  fprintf(fp,"# Estimated gamma/v with standard error: %g +/- %g\n",
	  gamma_over_v_hat(c), gamma_over_v_se(c,cov));
  fprintf(fp,"# Estimates kappa_S with standard error (using error propagation): %g +/- %g\n",
	  kappa_S, kappa_S_se);
  fprintf(fp,"# kappa_S confidence intervals based on parametric bootstrap\n");
  fprintf(fp,"# 0.95 CI for kappa_S: [%g,%g]\n",kappa_S_l95,kappa_S_u95);
  fprintf(fp,"# 0.99 CI for kappa_S: [%g,%g]\n",kappa_S_l99,kappa_S_u99);
  fprintf(fp,"\n\n");
  fprintf(fp,"# The data\n"
	  "# κ_Fura  τ   τ_SE\n");		
  for (size_t s_idx=0; s_idx<n_stim; s_idx++)	
    fprintf(fp,"%g %g %g\n",
	    gsl_vector_get(kappa,s_idx),	
	    gsl_vector_get(tau,s_idx),		
	    (fit_res[s_idx]).tau_se);		
  fprintf(fp,"\n\n");
  fprintf(fp,
	  "# The fitted data\n"
	  "# Preticted values with lower and upper bounds of 0.95 CI\n"
	  "# κ_Fura  τ   τ_lwr  τ_upr\n");
  double kappa_bd = 0.0;
  if (kappa_S > 0.0)
    kappa_bd = 1.25*kappa_S+1;
  // double kappa_range = 1.05*gsl_vector_get(kappa,n_stim-1)+kappa_S_u99+1; 
  double kappa_range = 1.05*gsl_vector_get(kappa,n_stim-1)+kappa_bd; 
  double delta_kappa = kappa_range/250;

  gsl_vector * kappa_v = gsl_vector_alloc(2);
  gsl_vector_set(kappa_v,0,1.0);
  for (size_t i=0; i<250; i++) {
    double kappaP = -kappa_bd+delta_kappa*i;
    gsl_vector_set(kappa_v,1,kappaP);
    double tauP,tauP_err;
    gsl_multifit_linear_est(kappa_v, c, cov, &tauP, &tauP_err);
    fprintf(fp,
	    "%g %g %g %g\n",
	    kappaP, tauP,
	    tauP-1.96*tauP_err,
	    tauP+1.96*tauP_err);
  }
  gsl_vector_free(kappa_v);
  fclose(fp);
  
  gsl_vector_free(tau);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  
  return 0;
}

/// @brief Performs a robust linear regression and prints results
///        to `stderr` and to a file whose prefix is built from
///        `output`.
///
/// @returns 0 if everything is fine
int robust_fit(const gsl_vector * kappa, ///< [in] vector of kappa values.  
	       mono_exp_fit_res * fit_res, ///< [in] vector of mono_exp_fit_res.
	       size_t n_stim, ///< [in] the number of stimulation.
	       char * filename, ///< Name of file containing the data
	       char * output, ///< Prefix of file name where results are written
	       gsl_rng * rng, ///< Pointer to a gsl_rng
	       int what, ///< What is processed? 0 = mean, 1 = min, 2 = max
	       size_t *stim ///< [in] Vector of indexes of used stimulation
  ) 
{
  fprintf(stderr,"\n***********************************************************\n");
  if (what==0)
    fprintf(stderr,"* Doing tau vs mean kappa_Fura robust bi-square regression *\n");
  if (what==1)
    fprintf(stderr,"* Doing tau vs min kappa_Fura robust bi-square regression *\n");
  if (what==2)
    fprintf(stderr,"* Doing tau vs max kappa_Fura robust bi-square regression *\n");
  fprintf(stderr,"***********************************************************\n");

  gsl_vector * c = gsl_vector_alloc(2);
  gsl_vector * c_ols = gsl_vector_alloc(2);
  gsl_matrix * cov = gsl_matrix_alloc(2,2);
  gsl_matrix * cov_ols = gsl_matrix_alloc(2,2);
  double tss;
  // Find out tau estimation with largest variance
  double max_var = 0.0;
  for (size_t s_idx=0; s_idx<n_stim; s_idx++) {
    double var = gsl_pow_2((fit_res[s_idx]).tau_se);
    if (var > max_var)
      max_var = var;
  }
  // Get the number of pseudo observations
  size_t pseudo_n_stim = 0;
  for (size_t s_idx=0; s_idx<n_stim; s_idx++) {
    double var = gsl_pow_2((fit_res[s_idx]).tau_se);
    pseudo_n_stim += (size_t) round(max_var/var);
  }
  // Prepare pseudo data for fit
  // Pseudo data are obtained from actual data by
  // replicating some observations in order to have
  // the right standard error ratio
  gsl_vector * tau = gsl_vector_alloc(pseudo_n_stim);
  gsl_vector * pkappa = gsl_vector_alloc(pseudo_n_stim);
  size_t p_idx = 0;
  for (size_t s_idx=0; s_idx<n_stim; s_idx++) {
    double var = gsl_pow_2((fit_res[s_idx]).tau_se);
    size_t ncp = (size_t) round(max_var/var);
    double tau_value = (fit_res[s_idx]).tau;
    double kappa_value = gsl_vector_get(kappa,s_idx);  
    for (size_t c_idx=0; c_idx < ncp; c_idx++) {
      gsl_vector_set(tau,p_idx,tau_value);
      gsl_vector_set(pkappa,p_idx,kappa_value);
      p_idx += 1;
    }
  }
  // Get total sum of squares (TSS)
  tss = gsl_stats_tss(tau->data,1,tau->size);
  // Allocate design matrix
  gsl_matrix * X = gsl_matrix_alloc(pseudo_n_stim,2);
  // Fill design matrix
  for (size_t s_idx=0; s_idx<pseudo_n_stim; s_idx++) {
    gsl_matrix_set(X,s_idx,0,1.0);
    gsl_matrix_set(X,s_idx,1,gsl_vector_get(pkappa,s_idx));
  }
  // Do robust with bi-square weight function
  gsl_multifit_robust_workspace * work = gsl_multifit_robust_alloc (gsl_multifit_robust_bisquare,
								    pseudo_n_stim, 2);
  gsl_multifit_robust (X, tau, c, cov, work);
  gsl_multifit_robust_stats rstats = gsl_multifit_robust_statistics(work);
  double chisq = rstats.sse;
  double Rsq = rstats.Rsq;
  gsl_multifit_robust_free (work);

  // Do OLS fit
  gsl_multifit_robust_workspace * work_ols = gsl_multifit_robust_alloc (gsl_multifit_robust_ols,
									pseudo_n_stim, 2);
  gsl_multifit_robust (X, tau, c_ols, cov_ols, work_ols);
  gsl_multifit_robust_stats rstats_ols = gsl_multifit_robust_statistics(work_ols);
  double chisq_ols = rstats_ols.sse;
  double Rsq_ols = rstats_ols.Rsq;
  gsl_multifit_robust_free (work_ols);
  

  double kappa_S_MC[NREP];
  double beta0 = BETA0(c);
  double beta1 = BETA1(c);
  double cov00 = BETA0_VAR(cov);
  double cov11 = BETA1_VAR(cov);
  double beta0_se = BETA0_SE(cov);
  double kappa_S = kappa_s(c);
  double kappa_S_se = kappa_s_se(c,cov);
  double se11=BETA1_SE(cov);
  for (size_t rep_idx=0; rep_idx<NREP; rep_idx++) {
    double num=(beta0+gsl_ran_gaussian_ziggurat(rng,beta0_se));
    double denom=(beta1+gsl_ran_gaussian_ziggurat(rng,se11));
    kappa_S_MC[rep_idx] = num/denom-1.0;
  }
  gsl_sort(kappa_S_MC,1,NREP);
  double kappa_S_l95=kappa_S_MC[(size_t)(0.025*NREP-1)];
  double kappa_S_u95=kappa_S_MC[(size_t)(0.975*NREP-1)];
  double kappa_S_l99=kappa_S_MC[(size_t)(0.005*NREP-1)];
  double kappa_S_u99=kappa_S_MC[(size_t)(0.995*NREP-1)];

  fprintf(stderr,"Best fit: tau = %g + %g kappa_Fura\n",beta0 , beta1 );
  fprintf(stderr,"Covariance matrix:\n");
  fprintf(stderr,"[ %+.5e, %+.5e  \n", cov00, gsl_matrix_get(cov,0,1));
  fprintf(stderr,"  %+.5e, %+.5e  ]\n", gsl_matrix_get(cov,1,0), cov11);
  fprintf(stderr,"Total sum of squares (TSS) = %g\n", tss);
  fprintf(stderr,"chisq (Residual sum of squares, RSS) = %g\n", chisq);
  fprintf(stderr,"R squared (1-RSS/TSS) = %g\n", Rsq);
  fprintf(stderr,"Estimated gamma/v with standard error: %g +/- %g\n",
	  gamma_over_v_hat(c), gamma_over_v_se(c,cov));
  fprintf(stderr,"Estimated kappa_S with standard error (using error propagation): %g +/- %g\n",
	  kappa_S, kappa_S_se);
  fprintf(stderr,"kappa_S confidence intervals based on parametric bootstrap\n");
  fprintf(stderr,"0.95 CI for kappa_S: [%g,%g]\n",kappa_S_l95,kappa_S_u95);
  fprintf(stderr,"0.99 CI for kappa_S: [%g,%g]\n",kappa_S_l99,kappa_S_u99);
  fprintf(stderr,"***********************************************************\n");
  if (what==0)
    fprintf(stderr,"* tau vs mean kappa_Fura robust bi-square regression done *\n");
  if (what==1)
    fprintf(stderr,"* tau vs min kappa_Fura robust bi-square regression done  *\n");
  if (what==2)
    fprintf(stderr,"* tau vs max kappa_Fura robust bi-square regression done  *\n");
  fprintf(stderr,"***********************************************************\n");

  // Print OLS results to stderr
  fprintf(stderr,"\n****************************************************\n");
  if (what==0)
    fprintf(stderr,"* Doing tau vs mean kappa_Fura ordinary regression *\n");
  if (what==1)
    fprintf(stderr,"* Doing tau vs min kappa_Fura ordinary regression *\n");
  if (what==2)
    fprintf(stderr,"* Doing tau vs max kappa_Fura ordinary regression *\n");
  fprintf(stderr,"*****************************************************\n");

  double beta0_ols = BETA0(c_ols);
  double beta1_ols = BETA1(c_ols);
  double cov00_ols = BETA0_VAR(cov_ols);
  double cov11_ols = BETA1_VAR(cov_ols);
  double beta0_ols_se = BETA0_SE(cov_ols);
  double kappa_S_ols = kappa_s(c_ols);
  double kappa_S_ols_se = kappa_s_se(c_ols,cov_ols);
  double se11_ols=sqrt(cov11_ols);
  for (size_t rep_idx=0; rep_idx<NREP; rep_idx++) {
    double num=(beta0_ols+gsl_ran_gaussian_ziggurat(rng,beta0_ols_se));
    double denom=(beta1_ols+gsl_ran_gaussian_ziggurat(rng,se11_ols));
    kappa_S_MC[rep_idx] = num/denom-1.0;
  }
  gsl_sort(kappa_S_MC,1,NREP);
  double kappa_S_ols_l95=kappa_S_MC[(size_t)(0.025*NREP-1)];
  double kappa_S_ols_u95=kappa_S_MC[(size_t)(0.975*NREP-1)];
  double kappa_S_ols_l99=kappa_S_MC[(size_t)(0.005*NREP-1)];
  double kappa_S_ols_u99=kappa_S_MC[(size_t)(0.995*NREP-1)];

  fprintf(stderr,"Best fit: tau = %g + %g kappa_Fura\n",beta0_ols , beta1_ols );
  fprintf(stderr,"Covariance matrix:\n");
  fprintf(stderr,"[ %+.5e, %+.5e  \n", cov00_ols, gsl_matrix_get(cov_ols,0,1));
  fprintf(stderr,"  %+.5e, %+.5e  ]\n", gsl_matrix_get(cov_ols,1,0), cov11_ols);
  fprintf(stderr,"Total sum of squares (TSS) = %g\n", tss);
  fprintf(stderr,"chisq (Residual sum of squares, RSS) = %g\n", chisq_ols);
  fprintf(stderr,"R squared (1-RSS/TSS) = %g\n", Rsq_ols);
  fprintf(stderr,"Estimated gamma/v with standard error: %g +/- %g\n",
	  gamma_over_v_hat(c_ols), gamma_over_v_se(c_ols,cov_ols));
  fprintf(stderr,"Estimated kappa_S with standard error (using error propagation): %g +/- %g\n",
	  kappa_S_ols, kappa_S_ols_se);
  fprintf(stderr,"kappa_S confidence intervals based on parametric bootstrap\n");
  fprintf(stderr,"0.95 CI for kappa_S: [%g,%g]\n",kappa_S_ols_l95,kappa_S_ols_u95);
  fprintf(stderr,"0.99 CI for kappa_S: [%g,%g]\n",kappa_S_ols_l99,kappa_S_ols_u99);
  fprintf(stderr,"****************************************************\n");
  if (what==0)
    fprintf(stderr,"* tau vs mean kappa_Fura ordinary regression done *\n");
  if (what==1)
    fprintf(stderr,"* tau vs min kappa_Fura ordinary regression done  *\n");
  if (what==2)
    fprintf(stderr,"* tau vs max kappa_Fura ordinary regression done  *\n");
  fprintf(stderr,"****************************************************\n");
  
  // Print to file
  char out[512];
  strcpy(out,output);
  if (what==0)
    strcat(out,"_tau_vs_mean_kappa_robust");
  if (what==1)
    strcat(out,"_tau_vs_min_kappa_robust");
  if (what==2)
    strcat(out,"_tau_vs_max_kappa_robust");
  FILE *fp = fopen(out,"w");
  if (what==0)
    fprintf(fp,"# τ vs mean κ robust regression for data set %s\n",filename);
  if (what==1)
    fprintf(fp,"# τ vs min κ robust regression for data set %s\n",filename);
  if (what==2)
    fprintf(fp,"# τ vs max κ robust regression for data set %s\n",filename);
  fprintf(fp,"# Using stimulation: %d",(int) stim[0]);
  for (size_t s_idx=1; s_idx<n_stim; s_idx++)
    fprintf(fp,", %d", (int) stim[s_idx]);
  fprintf(fp,"\n");
  fprintf(fp,"# Results obtained with robust bi-square regression:");
  fprintf(fp,"# Best fit: tau = %g + %g kappa_Fura\n",beta0 , beta1 );
  fprintf(fp,"# Covariance matrix:\n");
  fprintf(fp,"# [ %+.5e, %+.5e  \n", cov00, gsl_matrix_get(cov,0,1));
  fprintf(fp,"#   %+.5e, %+.5e  ]\n", gsl_matrix_get(cov,1,0),cov11);
  fprintf(fp,"# Total sum of squares (TSS) = %g\n", tss);
  fprintf(fp,"# chisq (Residual sum of squares, RSS) = %g\n", chisq);
  fprintf(fp,"# R squared (1-RSS/TSS) = %g\n", Rsq);
  fprintf(fp,"# Estimated gamma/v with standard error: %g +/- %g\n",
	gamma_over_v_hat(c), gamma_over_v_se(c,cov));
  fprintf(fp,"# Estimates kappa_S with standard error (using error propagation): %g +/- %g\n",
	  kappa_S, kappa_S_se);
  fprintf(fp,"# kappa_S confidence intervals based on parametric bootstrap\n");
  fprintf(fp,"# 0.95 CI for kappa_S: [%g,%g]\n",kappa_S_l95,kappa_S_u95);
  fprintf(fp,"# 0.99 CI for kappa_S: [%g,%g]\n",kappa_S_l99,kappa_S_u99);
  fprintf(fp,"\n\n");
  fprintf(fp,"# Results obtained with ordinary regression:\n");
  fprintf(fp,"# Best fit: tau = %g + %g kappa_Fura\n",beta0_ols , beta1_ols );
  fprintf(fp,"# Covariance matrix:\n");
  fprintf(fp,"# [ %+.5e, %+.5e  \n", cov00_ols, gsl_matrix_get(cov_ols,0,1));
  fprintf(fp,"#   %+.5e, %+.5e  ]\n", gsl_matrix_get(cov_ols,1,0), cov11_ols);
  fprintf(fp,"# Total sum of squares (TSS) = %g\n", tss);
  fprintf(fp,"# chisq (Residual sum of squares, RSS) = %g\n", chisq_ols);
  fprintf(fp,"# R squared (1-RSS/TSS) = %g\n", Rsq_ols);
  fprintf(fp,"# Estimated gamma/v with standard error: %g +/- %g\n",
	  gamma_over_v_hat(c_ols), gamma_over_v_se(c_ols,cov_ols));
  fprintf(fp,"# Estimates kappa_S with standard error (using error propagation): %g +/- %g\n",
	  kappa_S_ols, kappa_S_ols_se);
  fprintf(fp,"# kappa_S confidence intervals based on parametric bootstrap\n");
  fprintf(fp,"# 0.95 CI for kappa_S: [%g,%g]\n",kappa_S_ols_l95,kappa_S_ols_u95);
  fprintf(fp,"# 0.99 CI for kappa_S: [%g,%g]\n",kappa_S_ols_l99,kappa_S_ols_u99);
  fprintf(fp,"\n\n");
  
  fprintf(fp,"# The data\n"
	  "# κ_Fura  τ   τ_SE\n");		
  for (size_t s_idx=0; s_idx<n_stim; s_idx++)	
    fprintf(fp,"%g %g %g\n",
	    gsl_vector_get(kappa,s_idx),	
	    (fit_res[s_idx]).tau,		
	    (fit_res[s_idx]).tau_se);		
  fprintf(fp,"\n\n");
  fprintf(fp,
	  "# The fitted data\n"
	  "# Preticted values with lower and upper bounds of 0.95 CI and the two fitting methods\n"
	  "# κ_Fura  τ   τ_lwr  τ_upr  τ_ols  τ_ols_lwr  τ_ols_upr\n");
  double kappa_range = 1.05*gsl_vector_get(pkappa,pseudo_n_stim-1)+kappa_S_u99+1; 
  double delta_kappa = kappa_range/250;

  gsl_vector * kappa_v = gsl_vector_alloc(2);
  gsl_vector_set(kappa_v,0,1.0);
  for (size_t i=0; i<250; i++) {
    double kappaP = -kappa_S_u99-1+delta_kappa*i;
    gsl_vector_set(kappa_v,1,kappaP);
    double tauP,tauP_err;
    gsl_multifit_robust_est(kappa_v, c, cov, &tauP, &tauP_err);
    double tauP_ols,tauP_ols_err;
    gsl_multifit_robust_est(kappa_v, c_ols, cov_ols, &tauP_ols, &tauP_ols_err);
    fprintf(fp,
	    "%g %g %g %g %g %g %g\n",
	    kappaP, tauP,
	    tauP-1.96*tauP_err,
	    tauP+1.96*tauP_err,
	    tauP_ols,
	    tauP_ols-1.96*tauP_ols_err,
	    tauP_ols+1.96*tauP_ols_err);
  }
  gsl_vector_free(kappa_v);
  fclose(fp);
  
  gsl_vector_free(tau);
  gsl_vector_free(pkappa);
  gsl_vector_free(c);
  gsl_vector_free(c_ols);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_matrix_free(cov_ols);
  
  return 0;
}
