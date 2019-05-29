#!/usr/bin/env python
# import modules
import os
import shutil
import subprocess
import random
import argparse

# parse program arguments
script_description = ("Runs aba_ratio in a systematic way on a given data set."
                      " All transients are analyzed first. The ones having a too "
                      "large residual sum of square (RSS) or a too large residuals "
                      "lag 1 correlation are detected and a new "
                      "aba_ratio run is performed with the other ones only. "
                      "The diagnostic figures are generaed in png format and "
                      "summary reports are generated in MarkDown and HTML "
                      "formats. These figure and reports, together with aba_ratio "
                      "output are stored in a directory called 'data_name_analysis' "
                      "created as a sub-directory of the one from which this "
                      "program is executed. The null disdribution of the residuals "
                      "lag 1 correlation is obtained by simulation (permutations "
                      "of the observed sequence of residuals).") 
parser = argparse.ArgumentParser(description=script_description)
parser.add_argument('-f', '--data-file-name',
                    dest='data_name',
                    help = ('The data file name without the ".h5" suffix.'),
                    required=True)
parser.add_argument('-d', '--data-directory',
                    dest='data_directory',
                    help = ('The directory name where the data are found.'),
                    required=True)
parser.add_argument('-p', '--path-to-code',
                    dest='path_to_aba_ratio',
                    help = ('The directory name where the code is found. '
                            'Only necessary if it not in the user PATH.'),
                    required=False)
# Get the baseline length
def _baseline_length(string):
    n = int(string)
    if n <= 0:
        msg = "The baseline length must be > 0"
        raise argparse.ArgumentTypeError(msg)
    return n
parser.add_argument('-l', '--baseline-length',
                    type=_baseline_length, dest='baseline_length',
                    help=('The baseline length for the fits '
                          '(default 7)'),
                    default=7)
# Set the (pseudo)random number generator seed
parser.add_argument('-s', '--rng-seed',
                    type=int, dest='seed',
                    help=('The random number generator seed '
                          '(default 18710305)'),
                    default=18710305)
# Get the number of permutations
def _nb_perm(string):
    n = int(string)
    if n <= 0:
        msg = "The number of permutations must be > 0"
        raise argparse.ArgumentTypeError(msg)
    return n
parser.add_argument('-r', '--number-of-permutations',
                    type=_nb_perm, dest='nperm',
                    help=('The number of permutations used for '
                          'the simulations (default 1000)'),
                    default=1000)
args = parser.parse_args()

data_name = args.data_name
data_dir = os.path.abspath(args.data_directory)
path = os.path.abspath(args.path_to_aba_ratio)
baseline_length = args.baseline_length
rng_seed = args.seed
random.seed(rng_seed)
nperm = args.nperm

# create directory where results will be kept
# Make directory to store results and copy data into it
new_dir = os.path.abspath(data_name+'_analysis')
if not os.path.exists(new_dir):  # Check if directory already exists
    os.makedirs(new_dir)  # if not, create it
os.chdir(new_dir)
shutil.copy(data_dir+'/'+ data_name +'.h5',data_name +'.h5')

# run the analysis of the whole set of transients
# run aba_ratio on the whole set of transients
print("Doing now analysis with all the transients.\n")
cmd = path + '/aba_ratio -i ' + data_name + '.h5 -g -b ' + str(baseline_length)
p = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
print(p.stderr.read().decode())

# get the total number of transients
# Figure out first the number of transients
with open(data_name + '_aba_tau_vs_mean_kappa','r') as fin:
    tau_vs_kappa = fin.read()

tau_vs_kappa.find('# Using stimulation:')
debut = tau_vs_kappa.find('# Using stimulation:')
fin = tau_vs_kappa.find('\n',debut)
working_string = tau_vs_kappa[(debut+len('# Using stimulation:')):fin]
nb_transients = int(working_string[-1])

# Get key statistics on each transient
# Next get: the RSS per dof for each transient
#           the key part of the fit info for each
#           as well as the standard error on tau
#           and statistics on the residuals correlation
# Create a list with the good transients        

def corr(lst: list) -> float:
    """Returns lag 1 correlation of its input.

    The elements of lst are assumed centered (no mean subtraction is
    done).
    """
    n = len(lst) - 1
    return sum([lst[i+1]*lst[i] for i in range(n)])/n


rss_per_dof = []  # a list of values
tau_se = []  # a list of values
res_corr_lag1 = []  # a list of tuples (correlation and Proba of a smaller value under the null 
fit_info = []  # a list strings
kept_transients = []  # a list of integers
for t_idx in range(1,nb_transients+1):
    # Get the residuals in order to study their lag 1 correlation
    residual = []
    for line in open(data_name + '_aba_RatioFit_s' + str(t_idx), 'r'):
        if ("#" in line) or ("\n" == line):
            continue
        else:
            residual.append(float(line.split(" ")[3]))
    residual_corr = corr(residual)  # lag 1 correlation of the residuals
    corr_null = [0]*nperm  # a list used for storing the corr. of the shuffled
                           # residuals
    for i in range(nperm):
        random.shuffle(residual)
        corr_null[i] = corr(residual)
    corr_null.sort()  # the corr. coef. of the shuffled residuals are now sorted
    for i in range(-1,-len(corr_null)-1,-1):
        # find out the index of the first corr. coef. smaller or equal than
        # the observed one
        if corr_null[i] <= residual_corr:
            break
    # Add to res_corr_lag1 a tuple with the observed value and its probability
    res_corr_lag1.append((residual_corr,1.0-min(1.0,(nperm+i+1)/nperm)))
    with open(data_name + '_aba_RatioFit_s' + str(t_idx),'r') as fin:
            RatioFit = fin.read()
    # Get the standard error on tau
    debut = RatioFit.find("# estimated tau ")
    debut = RatioFit.find(" and standard error ",debut)
    fin = RatioFit.find("\n",debut)
    debut = debut+len(" and standard error ")
    tau_se.append(float(RatioFit[debut:fin]))
    # Get the general information on the fit
    debut = RatioFit.find('# nobs = ')
    fin = RatioFit.find('# rss per degree of freedom: ')
    fit_info.append(RatioFit[debut:fin])
    # Get the RSS per DOF
    debut = fin + len('# rss per degree of freedom: ')
    fin = RatioFit.find('\n',debut)
    rss_per_dof.append(float(RatioFit[debut:fin]))  
    # Find out if the fit was good enough
    condition1 = RatioFit.find("WARNING: THE FIT IS NOT GOOD!") == -1
    condition2 = RatioFit.find(("Probability of observing a larger of equal "
                                "RSS per DOF under the null hypothesis: "
                                "-nan")) == -1
    condition3 = RatioFit.find(("Probability of observing a larger of equal "
                                "RSS per DOF under the null hypothesis: "
                                "nan")) == -1
    condition4 = res_corr_lag1[-1][1] > 0.01
    if condition1 and condition2 and condition3 and condition4:
        kept_transients.append(t_idx)

    
# If some transients are not well fitted redo the
# analysis using only the good ones
# if there are bad transients redo analysis with the good ones
if len(kept_transients) < nb_transients:
    if len(kept_transients) >= 3:
        print("Doing now analysis with the 'good' transients.\n")
        cmd += ' -s ' + str(kept_transients).strip('[]').replace(' ','')
        p = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
        stderr_msg = p.stderr.read().decode()


# Generates the figures and the output summary file
# in markdown (MD) format
fout = open(data_name+'_report.md','w')
fout.write("*Analysis of dataset " + data_name + "*\n")
fout.write("-----\n\n")
fout.write("[TOC]\n\n")

fout.write("The baseline length is: {0:d}.\n\n".format(baseline_length))
fout.write(("**When fitting tau against kappa_Fura only the transients "
            "for which the fit RSS and the lag 1 auto-correlation "
            "of the residuals were small enough, giving an overall "
            "probability of false negative of 0.02, were kept** (see "
            "the numerical summary associated with each transient).\n\n"))

fout.write("\nThe good transients are: " + str(kept_transients).strip('[]') + ".\n\n")
if len(kept_transients) < 3:
    fout.write("**Not enough good transients to keep going!**\n\n")

files = os.listdir()
files = [f for f in files if '.gp' in f]

def mkfig(f: str) -> str:
    """Makes figures from gnuplot script file f and returns file name of created figure.
    """
    gp_name = f
    gp_out = f.replace(".gp",".png")
    gp_cmd = ("set terminal pngcairo size 800,800 enhanced font 'Verdana,12';\n"
              "set o '" + gp_out + "';\n"
              "load '" + gp_name + "'\n")
    subprocess.run('gnuplot',input=gp_cmd,text=True)
    return gp_out

# Start with the loading curve
loading_name = [f for f in files if 'loading_curve' in f][0]
fout.write("# Loading curve\n")
fout.write(("The time at which the 'good' transients were recorded"
            " appear in red.\n\n"))
gp_out = mkfig(loading_name)
fout.write("!["+loading_name.strip(".gp")+"]("+gp_out+")\n\n")

# The transient
fout.write("# Transients \n")
fout.write(("On each graph, the residuals appear on top.\n"
            "**Under the null hypothesis**, if the monoexponential fit is correct "
            "**they should be centered on 0 and have a SD close to 1** (not "
            "exactly 1 since parameters were obtained through the fitting "
            "procedure form the data.\n\n"
            "The estimated [Ca2+] appears on the second row. The estimate "
            "is show in black together with pointwise 95% confidence intervals."
            " The fitted curve appears in red. **The whole transient is not "
            "fitted**, only a portion of it is: a portion of the baseline "
            "made of " + str(baseline_length) + " points and "
            "the decay phase starting at the time where the Delta[Ca2+] "
            "has reached 50% of its peak value.\n\n"
            "The time appearing on the abscissa is the time from the beginning "
            "of the experiment.\n\n"))
for t_idx in range(1,nb_transients+1):
    fout.write("## Transient " + str(t_idx) + "\n")
    transient_name = [f for f in files if 'RatioFit_s'+str(t_idx) in f][0]
    if t_idx in kept_transients:
        fout.write("**Transient "+str(t_idx)+" is 'good'.**\n\n")
    else:
        fout.write("**Transient "+str(t_idx)+" is a 'bad'.**\n\n")
    fout.write("### Fit graphical summary\n")
    gp_out = mkfig(transient_name)
    fout.write("!["+transient_name.strip(".gp")+"]("+gp_out+")\n\n")
    fout.write("### Fit numerical summary\n")
    #fout.write("\n\n~~~~~\n\n")
    fout.write("\n\n")
    fout.write(fit_info[t_idx-1].replace("#",">").replace("\n","\n\n"))
    fout.write("> Lag 1 residuals auto-correlation: {0:.3f}\n\n".\
               format(res_corr_lag1[t_idx-1][0]))
    fout.write("> Pr[Lag 1 auto-corr. > {0:.3f}] = {1:.3f}\n\n".\
               format(*res_corr_lag1[t_idx-1]))
    #fout.write("\n\n~~~~~\n\n")
    fout.write("\n\n")
if len(kept_transients) >= 3:
    # tau vs kappa figures
    fout.write("# tau vs kappa \n")
    fout.write(("Since the [Fura] changes during a transient (and it can change a lot "
                "during the early transients), the _unique_ value to use as '[Fura]' "
                "is not obvious. We therefore perform 3 fits: one using the minimal "
                "value, one using the mean and one using the maximal value.\n\n"
                "The observed tau (shown in red) are displayed with a 95% confidence "
                "interval that results from the fitting procedure and _is_ therefore "
                "_meaningful only if the fit is correct_!\n\n"
                "No serious attempt at quantifying the precision of [Fura] and "
                "therefore kappa_Fura has been made since the choice of which [Fura] "
                "to use has a larger effect and since the other dominating effect "
                "is often the certainty we can have that the saturating value (the "
                "[Fura] in the pipette) has been reached.\n\n"
                "The straight line in black is the result of a _weighted_ linear "
                "regression. The blue dotted lines correspond to the limits of "
                "_pointwise 95% confidence intervals_.\n\n"))
    suffix = ["min","mean","max"]
    for s in suffix:
        fout.write("## tau vs kappa  using the " + s + " [Fura] value\n")
        fout.write("### Fit graphical summary\n")
        tau_name = [f for f in files if 'tau_vs_'+s in f][0]
        gp_out = mkfig(tau_name)
        fout.write("!["+tau_name.strip(".gp")+"]("+gp_out+")\n\n")
        fout.write("### Fit numerical summary\n")
        with open(data_name + '_aba_tau_vs_' + s + '_kappa','r') as fin:
            tau_vs_kappa = fin.read()
            debut = tau_vs_kappa.find("# Best fit:")
            fin = tau_vs_kappa.find("\n\n# The data")
            #fout.write("\n\n~~~~~\n\n")
            fout.write("\n\n")
            fout.write(tau_vs_kappa[debut:fin].replace("#",">").replace("\n","\n\n"))
            #fout.write("\n\n~~~~~\n\n")
            fout.write("\n\n")
    
# write the RSS per dof for each good transient
rss_per_dof = [rss_per_dof[i] for i in range(nb_transients) if i+1 in kept_transients]
tau_se = [tau_se[i] for i in range(nb_transients) if i+1 in kept_transients]
res_corr_lag1 = [res_corr_lag1[i] for i in range(nb_transients) if i+1 in kept_transients]
corr_lag1 = [res_corr_lag1[i][0] for i in range(len(res_corr_lag1))]
corr_prob = [res_corr_lag1[i][1] for i in range(len(res_corr_lag1))]
fout.write(("# RSS per DOF, standard error of tau and lag 1 residual "
            "correlation for each 'good' tansient\n"))
fout.write("{0:d} out of {1:d} transients  were kept.\n\n".\
           format(len(kept_transients), nb_transients))
fout.write("sigma(tau): "+str(tau_se).strip('[]')+"\n\n")
fout.write("Residual correlation at lag 1: "+str(corr_lag1).strip('[]')+"\n\n")
fout.write(("Probablity of a correlation at lag 1 smaller or equal than "
            "observed: ")+str(corr_prob).strip('[]')+"\n\n")
fout.write("RSS/DOF: "+str(rss_per_dof).strip('[]')+"\n")


# Do some clean up
fout.close()
os.remove(data_name + ".h5")
    
# Try to generate an HTML version of the MD summary file
# Try generating the html output with the markdown module
# You can install it since it is not part of the standard library
# with: pip install markdown
# See https://python-markdown.github.io/
try:
    import markdown
    import codecs  # part of the standard library
    fin = codecs.open(data_name+"_report.md", mode="r", encoding="utf-8")
    md = fin.read()
    fin.close()
    html = markdown.markdown(md,extensions=["extra","toc"])
    fout = codecs.open(data_name+"_report.html", mode="w",
                       encoding="utf-8", errors="xmlcharrefreplace")
    fout.write(html)
    fout.close()
except ModuleNotFoundError:
    # check if pandoc is installed
    pandoc = subprocess.run("pandoc -v", shell=True)
    if pandoc.returncode == 0:
        # if yes, call pandoc
        pandoc_cmd = 'pandoc -s ' + data_name + '_report.md -o ' + data_name + '_report.html'
        subprocess.call(pandoc_cmd, shell=True)
    else:
        print(("Neither 'markdown' (Python module), nor 'pandoc' found,"
               " so no 'html' output was generated."))
