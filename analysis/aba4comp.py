#!/usr/bin/env python
# Imports required modules
import os
from multiprocessing import Pool
from statistics import mean, stdev 
from math import sqrt
from collections import OrderedDict
import markdown # Not part of standard lib. must be installed by user
import codecs  

# Makes sure that path are adapted to OS
DA_beta_dir = os.path.abspath('../data_paper/data_beta_escin')
DA_wc_dir = os.path.abspath('../data_paper/data_whole_cell')
path2aba_ratio = os.path.abspath('../code')
path2aba_boring = os.path.abspath('.')  

# beta-escin data set analysis
new_dir = os.path.abspath('DA-beta')
if not os.path.exists(new_dir):  # Check if directory already exists
    os.makedirs(new_dir)  # if not, create it
os.chdir(new_dir)
DAfiles = [file for file in os.listdir(DA_beta_dir) if file.endswith(".h5")]

def call_aba_boring(file_name):
    import subprocess
    da_cmd = (path2aba_boring + "/aba_boring.py "
              "-f " + file_name[:file_name.find(".h5")] + " "
              "-d " + DA_beta_dir + " "
              "-p " + path2aba_ratio)
    subprocess.call(da_cmd,shell=True)

with Pool() as p:
    p.map(call_aba_boring, DAfiles)

#  Go back to where we started
os.chdir('../')  

# whole-cell data set analysis
new_dir = os.path.abspath('DA-wc')
if not os.path.exists(new_dir):  # Check if directory already exists
    os.makedirs(new_dir)  # if not, create it
os.chdir(new_dir)

DAfiles = [file for file in os.listdir(DA_wc_dir) if file.endswith(".h5")]

def call_aba_boring(file_name):
    import subprocess
    da_cmd = (path2aba_boring + "/aba_boring.py "
              "-f " + file_name[:file_name.find(".h5")] + " "
              "-d " + DA_wc_dir + " "
              "-p " + path2aba_ratio)
    subprocess.call(da_cmd,shell=True)

with Pool() as p:
    p.map(call_aba_boring, DAfiles)

#  Go back to where we started
os.chdir('../')  

#  Define read_report function
def read_report(report_name: str) -> dict:
    """Read report in directory report_name and returns a dictionary with 
    relevant statitics.
    
    Parameters
    ----------
    report_name: the name of the directory containing the report (if
                 'DA_120906_E1_analysis' a file called 'DA_120906_E1_reprot.md'
                 will be looked for in that directory)
    
    Returns
    -------
    A dictionary with the following keys:
    nb_used: the number of 'usable' transients
    nb_total: the total number of transients
    setting: the best settinf (min, mean or max), ie, the setting yielding the
             best tau vs kappa (straight line) fit
    RSS: the residual sum of squares of the best fit
    Prob: the probability for the Chi2 distribution to exceed the RSS under 
          the null hypothesis
    CI: the 99% confidence interval for kappa_s for the bset fit (using a
        bootstrap).
    """
    from collections import OrderedDict
    import math  # get nan and inf
    fn = report_name.replace("_analysis","_report.md")
    with open(report_name+'/'+fn,'r') as fin:
        report = fin.read() 
    def mk_get_value(report: str):
        """Takes a strings and returns a function (more precisely a closure)
        making it easy to extract specific information contained in the strig.
        """
        import math
        def get_value(start: int,
                      bra: str,
                      cket: str = "\n",
                      what: str = "int") -> tuple:
            """Reads a value from report between bra and cket from start
            and returns a tuple containg the point where it stopped and
            the requested value.
            
            Pamaters
            --------
            start: a positive intege, the  position from which searched is done.
            bra: a string preceeding immediatly the value.
            cket: a string following the value (in other words we should have
                  bra value cket).
            what: a string, either 'int' or 'float' depending on the type of
                  value
    
            Returns
            -------
            A tuple with the position in report where cket was found and the
            value.
            """
            start = report.find(bra,start)+len(bra)
            end = report.find(cket,start)
            val = report[start:end]
            if what == "int":
                res = int(val)
            else:
                res = float(val)
            return end,res
        return get_value
    
    get_value = mk_get_value(report)
    get_value_par = {"start": 0,
                     "bra": ("# RSS per DOF, standard error of tau and lag 1 "
                             "residual correlation for each 'good' tansient\n"),
                     "cket": " out of ",
                     "what": "int"}
    start, nb_used = get_value(**get_value_par)
    get_value_par = {"start": start,
                     "bra": " out of ",
                     "cket": " transients",
                     "what": "int"}
    start, nb_total = get_value(**get_value_par)
        
    transients = OrderedDict()
    dict_lst = [{"bra": "> nobs = ",
                 "cket": "\n",
                 "what": "int"},  # n_obs
                {"bra": "> number of degrees of freedom = ",
                 "cket": "\n",
                 "what": "int"},  # n_dof
                {"bra": "> baseline length = ",
                 "cket": "\n",
                 "what": "int"},  # baseline_length
                {"bra": "> fit started from point ",
                 "cket": "\n",
                 "what": "int"},  # fit_start
                {"bra": "> estimated baseline ",
                 "cket": " and standard error ",
                 "what": "float"},  # Ca0_best
                {"bra": " and standard error ",
                 "cket": "\n",
                 "what": "float"},  # Ca0_se
                {"bra": "> estimated delta ",
                 "cket": " and standard error ",
                 "what": "float"},  # Delta_best
                {"bra": " and standard error ",
                 "cket": "\n",
                 "what": "float"},  # Delta_se
                {"bra": "> estimated tau ",
                 "cket": " and standard error ",
                 "what": "float"},  # Tau_best
                {"bra": " and standard error ",
                 "cket": "\n",
                 "what": "float"},  # Tau_se
                {"bra": "> residual sum of squares: ",
                 "cket": "\n",
                 "what": "float"},  # rss
                {"bra": "> RSS per degree of freedom: ",
                 "cket": "\n",
                 "what": "float"},  # rss_per_dof
                {"bra": ("> Probability of observing a larger of equal RSS "
                         "per DOF under the null hypothesis: "),
                 "cket": "\n",
                 "what": "float"},  # prob_rss
                {"bra": "> Lag 1 residuals auto-correlation: ",
                 "cket": "\n",
                 "what": "float"},  # lag1
                {"bra": "] = ",
                 "cket": "\n",
                 "what": "float"}]  # prob_lag1
                
    for t_idx in range(1,nb_total+1):
        start=report.find("## Transient "+str(t_idx))
        start=report.find("### Fit numerical summary",start)
        get_value_par = {"start": start}
        val = []
        for par in dict_lst:
            get_value_par.update(par)
            get_value_par["start"] = start
            start, res = get_value(**get_value_par)
            val.append(res)
        transients["Transient"+str(t_idx)]={"n_obs": val[0],
                                            "n_dof": val[1],
                                            "baseline_length": val[2],
                                            "fit_start": val[3],
                                            "Ca0_best": val[4],
                                            "Ca0_se": val[5],
                                            "Delta_best": val[6],
                                            "Delta_se": val[7],
                                            "tau_best": val[8],
                                            "tau_se": val[9],
                                            "rss": val[10],
                                            "rss_per_dof": val[11],
                                            "prob_rss": val[12],
                                            "lag1": val[13],
                                            "prob_lag1": val[14]}
    tau_vs_kappa = OrderedDict()
    dict_lst = [{"bra": "> chisq (Residual sum of squares, RSS) = ",
                 "cket": "\n",
                 "what": "float"},  # rss
                {"bra": ("> Probability of observing a larger of equal "
                         "RSS per DOF under the null hypothesis: ")  ,
                 "cket": "\n",
                 "what": "float"},  # prob
                {"bra": "> Estimated gamma/v with standard error: ",
                 "cket": " +/- ",
                 "what": "float"},  # gamma_best
                {"bra": " +/- ",
                 "cket": "\n",
                 "what": "float"},  # gamma_se
                {"bra": ("> Estimates kappa_S with standard "
                         "error (using error propagation): "),
                 "cket": " +/- ",
                 "what": "float"},  # kappa_S_best
                {"bra": " +/- ",
                 "cket": "\n",
                 "what": "float"},  # kappa_S_se
                {"bra": "> 0.99 CI for kappa_S: [",
                 "cket": ",",
                 "what": "float"},  # 99% CI lower
                {"bra": ",",
                 "cket": "]\n",
                 "what": "float"}]  # 99% CI upper
    settings = ["min","mean","max"]
    best_setting = "min"
    for suffix in settings:
        if nb_used < 3:
            tau_vs_kappa["tau_vs_kappa_"+suffix]={"rss": "NA",
                                                  "prob": "NA",
                                                  "gamma_best": "NA",
                                                  "gamma_se": "NA",
                                                  "kappa_S_best": "NA",
                                                  "kappa_S_se": "NA",
                                                  "CI_lower": "NA",
                                                  "CI_upper": "NA"}
        else:
            bra = "## tau vs kappa  using the "+suffix+" [Fura] value\n"
            start = report.find(bra)
            get_value_par = {"start": start}
            val = []
            for par in dict_lst:
                get_value_par.update(par)
                get_value_par["start"] = start
                start, res = get_value(**get_value_par)
                val.append(res)
            if suffix == "min":
                best_rss = val[0]
            elif val[0] < best_rss:
                best_rss = val[0]
                best_setting = suffix
            CI = "["+str(round(val[6]))+","+str(round(val[7]))+"]"
            tau_vs_kappa["tau_vs_kappa_"+suffix]={"RSS": val[0],
                                                  "Prob": val[1],
                                                  "gamma_best": val[2],
                                                  "gamma_se": val[3],
                                                  "kappa_S_best": val[4],
                                                  "kappa_S_se": val[5],
                                                  "CI_lower": round(val[6]),
                                                  "CI_upper": round(val[7]),
                                                  "CI": CI}
    
    return {"nb_used": nb_used,
            "nb_total": nb_total,
            "setting": best_setting,
            "transients": transients,
            "tau_vs_kappa": tau_vs_kappa}
  

# reads beta-escin analysis results
beta_dir = os.path.abspath('DA-beta')
os.chdir(beta_dir)
the_list = os.listdir(beta_dir)
the_list = [n for n in the_list if 'DA_' in n]
the_list.sort()
beta = OrderedDict([(dn.replace("_analysis",""), read_report(dn)) for dn in the_list])
#  Go back to where we started
os.chdir('../')  

# reads whole-cell analysis results
wc_dir = os.path.abspath('DA-wc')
os.chdir(wc_dir)
the_list = os.listdir(wc_dir)
the_list = [n for n in the_list if 'DA_' in n]
the_list.sort()
wc = OrderedDict([(dn.replace("_analysis",""), read_report(dn)) for dn in the_list])
#  Go back to where we started
os.chdir('../')  

# writes analysis summary to disk in md format
md="# Summaries for the beta-escin and whole-cell data sets\n"
md+="[TOC]\n\n"
md+="## The content of the tables\n"
md+=("In the next two tables, the columns have the following meaning:\n\n"
     "+ **Data set** gives the data set name.\n"
     "+ **used/total** displays the number of 'good' transients over the "
     "total number of transients.\n"
     "+ **setting** is the [Fura] concentration choice (among 'min', 'mean' "
     "and 'max') that yielded the best straight line fit.\n"
     "+ **RSS** is the residual sum of squares of the best fit (its magnitude "
     "depends on the number of 'good' transients in the dataset).\n"
     "+ **Prob** is the probability for the RSS to *exceed* the observed "
     "value under the null hypothesis.\n"
     "+ **99% CI** is the 99% confidence interval for the endogenous kappa "
     "(kappa_s) obtained with a bootstrap simulation.\n\n"
     "**Warning** if the data set yieleded less than 3 'good' transients, no "
     "straight line fit was performed (it would have been meaningless) and "
     "'NA' (Not Available) appears in the last 4 columns of the table.\n\n")

md+="## Analysis summary of beta-escin data sets\n\n"
md+="**tau vs kappa(Fura) straight line regression summary**\n\n"
tbl_head=["Data set","used/total","setting","RSS","Prob","99% CI"]
md+="~~~~~\n"
md+="{0:^25}|{1:^14}|{2:^11}|{3:^9}|{4:^9}|{5:^25}\n".format(*tbl_head)
md+="\n"
for data_name,value in beta.items():
    link="[{0}_report]({0}_analysis/{0}_report.html) ".format(data_name)
    ratio="{nb_used}/{nb_total}".format(**value)
    md+="{0:^26}{1:^15}".format(data_name,ratio)
    if value["nb_used"] < 3:
        md+="{0:^12}{0:^10}{0:^10}{0:^26}".format("NA")
    else:
        par = {"nb_used": value['nb_used'],
               "nb_total": value['nb_total'],
               "setting": value['setting']}
        good_tau = value['tau_vs_kappa']['tau_vs_kappa_'+value['setting']]
        par.update(good_tau)
        md+="{setting:^11s}{RSS:^10.2f}{Prob:^10.2f}{CI:^26s}".format(**par)
    md+="\n"

md+="~~~~~\n"

md+="## Analysis summary of whole-cell data sets\n\n"
md+="**tau vs kappa(Fura) straight line regression summary**\n\n"
tbl_head=["Data set","used/total","setting","RSS","Prob","99% CI"]
md+="~~~~~\n"
md+="{0:^25}|{1:^14}|{2:^11}|{3:^9}|{4:^9}|{5:^25}\n".format(*tbl_head)
md+="\n"
for data_name,value in wc.items():
    link="[{0}_report]({0}_analysis/{0}_report.html) ".format(data_name)
    ratio="{nb_used}/{nb_total}".format(**value)
    md+="{0:^26}{1:^15}".format(data_name,ratio)
    if value["nb_used"] < 3:
        md+="{0:^12}{0:^10}{0:^10}{0:^26}".format("NA")
    else:
        par = {"nb_used": value['nb_used'],
               "nb_total": value['nb_total'],
               "setting": value['setting']}
        good_tau = value['tau_vs_kappa']['tau_vs_kappa_'+value['setting']]
        par.update(good_tau)
        md+="{setting:^11s}{RSS:^10.2f}{Prob:^10.2f}{CI:^26s}".format(**par)
    md+="\n"

md+="~~~~~\n"


md+="# Transients numerical summaries of the beta-escin data sets\n"
dir_name = "DA-beta"
for data_name,value in beta.items():
    report_name = data_name+"_report"
    full_report_name = dir_name+"/"+data_name+"_analysis/"+report_name+".html"
    txt = "## Data set "+data_name+" numerical summary\n"
    txt += ("The full report with figures is located at ["
            "{report_name:s}]({full_report_name:s}).\n\n")
    txt = txt.format(**{"report_name": report_name,
                        "full_report_name": full_report_name})
    txt += ("This data set contains {nb_total:d} transients, "
            "{nb_used:d} of which were 'good' enough for the "
            "tau vs kappa fit.\n\n")
    txt = txt.format(**value)
    txt += ("The RSS per degree of freedom were (values in '()' correspond "
            "to 'bad' transients): ")
    nb_total = value['nb_total']
    for t_idx in range(1,nb_total+1):
        transient = value['transients']['Transient'+str(t_idx)]
        rss_per_dof = transient['rss_per_dof']
        if transient['prob_rss'] > 0.01:
            txt += "{0:.3f}".format(rss_per_dof)
        else:
            txt += "({0:.3f})".format(rss_per_dof)
        if t_idx < nb_total:
            txt += ", "
        else:
            txt += ".\n\n"
    txt += ("The lag 1 residual autocorrelation coefficient were "
            "(values in '()' correspond to 'bad' transients): ")
    for t_idx in range(1,nb_total+1):
        transient = value['transients']['Transient'+str(t_idx)]
        lag1 = transient['lag1']
        if transient['prob_rss'] > 0.01:
            txt += "{0:.3f}".format(lag1)
        else:
            txt += "({0:.3f})".format(lag1)
        if t_idx < nb_total:
            txt += ", "
        else:
            txt += ".\n\n"
    
    md += txt


md+="# Transients numerical summaries of the whole-cell data sets\n"
dir_name = "DA-wc"
for data_name,value in wc.items():
    report_name = data_name+"_report"
    full_report_name = dir_name+"/"+data_name+"_analysis/"+report_name+".html"
    txt = "## Data set "+data_name+" numerical summary\n"
    txt += ("The full report with figures is located at ["
            "{report_name:s}]({full_report_name:s}).\n\n")
    txt = txt.format(**{"report_name": report_name,
                        "full_report_name": full_report_name})
    txt += ("This data set contains {nb_total:d} transients, "
            "{nb_used:d} of which were 'good' enough for the "
            "tau vs kappa fit.\n\n")
    txt = txt.format(**value)
    txt += ("The RSS per degree of freedom were (values in '()' correspond "
            "to 'bad' transients): ")
    nb_total = value['nb_total']
    for t_idx in range(1,nb_total+1):
        transient = value['transients']['Transient'+str(t_idx)]
        rss_per_dof = transient['rss_per_dof']
        if transient['prob_rss'] > 0.01:
            txt += "{0:.3f}".format(rss_per_dof)
        else:
            txt += "({0:.3f})".format(rss_per_dof)
        if t_idx < nb_total:
            txt += ", "
        else:
            txt += ".\n\n"
    txt += ("The lag 1 residual autocorrelation coefficient were "
            "(values in '()' correspond to 'bad' transients): ")
    for t_idx in range(1,nb_total+1):
        transient = value['transients']['Transient'+str(t_idx)]
        lag1 = transient['lag1']
        if transient['prob_rss'] > 0.01:
            txt += "{0:.3f}".format(lag1)
        else:
            txt += "({0:.3f})".format(lag1)
        if t_idx < nb_total:
            txt += ", "
        else:
            txt += ".\n\n"
    
    md += txt
  

# converts md to html
html = markdown.markdown(md,extensions=["extra","toc"])
fout = codecs.open('beta_wc_comp.html',mode='w',
                   encoding='utf-8',errors='xmlcharrefreplace')
fout.write(html)
fout.close()
