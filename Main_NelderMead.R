# Core model: head and tail of glycolysis, lactate dehydrogenase, oxidative phosphorylation,
#             ATP hydrolysis
# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

# clear R workspace
rm(list=ls())
options(error = browser)

# load required libraries
# Package 'deSolve':  https://cran.r-project.org/web/packages/deSolve/index.html
library(deSolve)
# Package 'xlsx': https://cran.r-project.org/web/packages/xlsx/index.html
if (Sys.getenv("JAVA_HOME")!="") Sys.setenv(JAVA_HOME="")
library(xlsx)
# Package 'compiler': http://www.inside-r.org/r-doc/compiler/cmpfun
library(compiler)
library(XLConnect)
options( java.parameters = "-Xmx4g" )
xlcFreeMemory()

y_mem <<- vector()
cat("\nSTART\n")

# initialize character-variable for the filename
filename <- character(0)
# read filename from the following small file
# this file should contain a line stating : filename <- "your_own_choice"
source("filename.R")

# Read Parameter file
source("Parameters.R")

# several sourced files contain function definitions
# the statements in these function definitions are not immediately executed

# function definition to calculate several additional parameters from given parameters
source("Calc_NApars.R")
# load flux equations
source("Function_definitions.R")

# load initial start values concentration variables to be integrated based on the ODE system
source("xstart.R")
# plot routine
source("plotOutLsoda.R")

# initialize the parameter array
p2ck<-p2ck_default
# initialize the start values for the variables
xstart<-xstart_default

# create directory "Results" if it does not exist
mainDir<-"."
subDir <- "Results"

if (!file.exists(subDir)){
    cat("\nDirectory 'Results' does not exist and is created\n")
    dir.create("./Results") # Windows specific
}

# next statements make it possible to examine concentration vector after R stops a run
y_aux <- xstart_default
names(y_aux) <- names(xstart_default)

# The model function with the differential equations and evaluation of kinetic equations is loaded
# It will be executed by the integration routine
source("Model_function.R")
# the model function is compiled for faster execution
lmodel <- cmpfun(model)

# The model function for the runin part, before glucose is added
source("Model_function_runin.R")
# the model function is compiled for faster execution
lmodel_runin <- cmpfun(model_runin)

# experimental data read from Excel files
dataCurve <- read.xlsx("Data_1.xlsx",1)
dataCurve2 <- read.xlsx("Data_2.xlsx",1)
# read maximal FBP accumulation at various levels of glucose addition from Data3.xlsxQ
data3 <- read.xlsx("Data_3.xlsx",1)
# read lactate and FBP accumulation after 5 and 10 sec
data4 <- read.xlsx("Data_4.xlsx",1)
# matrices to store the differences (model result minus experimental results) are initialized
diffCurve <- matrix(nrow=dim(dataCurve)[1],ncol=5)
diffCurve2 <- matrix(nrow=dim(dataCurve2)[1],ncol=6)
diffCurve3 <- matrix(nrow=dim(data3)[1],ncol=5)
diffCurve4 <- matrix(nrow=dim(data4)[1],ncol=18)
# the columns of the calculated differences are given appropriate names
colnames(diffCurve) <- c("time","lactate","FBP","ATP","glucose")
colnames(diffCurve2) <- c("time","lactate","FBP","ATP","glucose","VO2")
colnames(diffCurve3) <- c("cell_concentration","glucose_concentration","FBP_exp","FBP_model","difference")
colnames(diffCurve4) <- c("cell_concentration","glucose_concentration","space1","lac_exp_5sec","lac_model_5sec","lac_diff_5sec",
                          "space2","lac_exp_10sec","lac_model_10sec","lac_diff_10sec","space3","FBP_exp_5sec","FBP_model_5sec",
                          "FBP_diff_5sec","space4","FBP_exp_10sec","FBP_model_10sec","FBP_diff_10sec")
# the first column, time,  set to the time points of the measurements
diffCurve[,1] <- dataCurve[,1]
diffCurve2[,1] <- dataCurve2[,1]

optpars <- read.xlsx("Optimized_parameters.xlsx",1)

# the names of the variables are stored in an array of characters to be used for output.
# The names are the same as for array y.
namesOutputVariables <- names(xstart_default)

# OPTIMIZATION EXECUTION SECTION

# specific data for simulation different from default read in
source(paste(filename,".R",sep=""))

uncoupler <- 0.  # is set to 1 only at the time when uncoupling of the mitochondria
                 # in the simulation protocol is required

names_xstart <- names(xstart)

####  flush.console() # works only under MS Windows. May increase running time very substantially.
                      # Uncomment if regular updates of results written to console are desirable.

end_time_run = p2ck["tbegin"]+p2ck["tduration1"]+p2ck["tduration2"]+p2ck["tduration3"]+p2ck["tduration4"]

# this section is used only if plots are made for each run during the ssq_priors call
# for plot for final calculation the framework is coded farther below
   if(p2ck["plot"] > 0.) {
    if (p2ck["plot_flux_pos_metabolites_neg"] >= 0) {
      plot(0,0,xlim=c(p2ck["tbegin"],end_time_run),
           xlab="Time (sec)", ylab="Net ATP synthesis (uM/s)",type="n",
           ylim=c(p2ck["y_axis_min_flux"],p2ck["y_axis_max_flux"]), col="blue")
           legend("topright", c("glycolysis", "oxidative phosphorylation"), fill=c("red", "blue"))
     }
    if (p2ck["plot_flux_pos_metabolites_neg"] < 0. ) {
      plot(0,0,xlim=c(p2ck["tbegin"],end_time_run),
           xlab="Time (sec)", ylab="Content (micromol/ml cell)",type="n",
           ylim=c(p2ck["y_axis_min_met"],p2ck["y_axis_max_met"]), col="blue")
           legend("topright", c("FBP","Del lactate","5 * ATP"), fill=c("red","blue"))
     }
    }

# TIMERS FOR THE SIMULATIONS
    out1_times <- p2ck["tbegin"]:0
    out2_times <- t(dataCurve[,"time"])
    out4_times <- t(dataCurve2[,"time"])
    out6_times <- c(seq(from=0,to=1,by=0.1),seq(from=2,to=60,by=1),seq(from=65,to=295,by=5),seq(from=300,to=3600,by=60)) #to=7200

index_AdNsynt <- -1
for (i in 1:length(out6_times)) {
  if (out6_times[i] == p2ck["index_time_AdNsynt"]) index_AdNsynt <- i 
}
if (index_AdNsynt < 0) stop("\n\nStopped becase index_time_AdNsynt is not in out6_times array")  

   out6_times_output <- c(seq(from=0,to=1,by=0.1),seq(from=2,to=60,by=1),seq(from=65,to=115,by=5),seq(from=120,to=3600,by=60)) #to=7200
                                          # time points for long period with high sugar

index_AdNsynt_out <- -1
for (i in 1:length(out6_times_output)) {
  if (out6_times_output[i] == p2ck["index_time_AdNsynt"]) index_AdNsynt_out <- i 
}
if (index_AdNsynt_out < 0) stop("\n\nStopped becase index_time_AdNsynt is not in out6_times_output array")  

direct_par_optim <- optpars[,"par_optim"]
names(direct_par_optim) <- optpars[,"names"]
#in the Optimized Parameters.xlsx file some parameters are chosen to be fixed
#some are chosen to be optimized
selected_par_optim <- as.logical(optpars[,"optimized"])
direct_par_optim_sel <- direct_par_optim[selected_par_optim]
direct_par_optim_fix <- direct_par_optim[!selected_par_optim]
par_optim <- log10(direct_par_optim_sel)

names_p2ck <- names(p2ck)
names_par <- names(par_optim)
names_par_fix <- names(direct_par_optim_fix)

for (par_name in names_par_fix) {
  p2ck[par_name] <- direct_par_optim_fix[par_name]
} 

ATPconc_start <- dataCurve2[1,"ATP"] / p2ck["Vcyt"] * 1000. /
                   ((1.-p2ck["DW_ratio"]) * p2ck["dens"])

lb <- log10(optpars[,"lb"])
lb <- lb[selected_par_optim]
names(lb) <- names_par
ub <- log10(optpars[,"ub"])
ub <- ub[selected_par_optim]
names(ub) <- names_par
invweight_lb <- optpars[,"invweight_lb"]
names(invweight_lb) <- names_par
invweight_ub <- optpars[,"invweight_ub"]
names(invweight_ub) <- names_par

# variables will be set to the initial glucose and lactate content, to be used repeatedly 
# in function ssq_priors()
# presently: no visible global binding. local variable by same name in ssq_priors gets assigned
glc_per_ml_initial <- 0.
lac_per_ml_initial <- 0.

# function to run the simulation over a number of desired periods and interventions
# provides among others log target density (ltd) which needs to be called if function spBayes is used
source("ssq_priors.R")

  lssq_priors <- cmpfun(ssq_priors)

  source(paste(filename,".R",sep=""))
  xstart<-xstart_default
  ssq_series <- numeric(0)
  # the number of calculated measures in par_optim_series, see below
  number_measures <- 13
  # the parameters tried in the intermediary steps of the Nelder-Mead series will be added to
  # par_optim_series in the ssq_priors function
  par_optim_series <- matrix(nrow=1,ncol=(length(par_optim)+number_measures))
  # to count the steps in the Nelder-Mead series
  count <- 0
  min_ssq_out2_4 <- 1e10
  par_min_ssq <- numeric()

  date_aux<-format(Sys.time(),"%d%b%Y_%H%M")
  # variable to print precise time of individual run on title Excel file and on plot. 
  # This variable is updated in function ssq_priors()
  date_aux1<-format(Sys.time(),"_%H.%M.%S")

# text output to screen is saved in file
  sink(file=paste(".\\Results\\run_output_",date_aux,".txt"),type="output",split=TRUE)

# Here the optimization run is done
  if(p2ck["ensemble_run"] > 0.) {
  # setting array.addition to TRUE means that the parameters for each step are added to par_optim_series by ssq_priors()
     array.addition <- TRUE

# OPTIMIZATION NELDER-MEAD
  # the log10 values of the parameters are given to the Nelder-Mead routine for better scaling.
  # This also prevents impossible negative parameter values
  # the parameter output is still log10 transformed
     optimized_par_Nelder_Mead <- optim(par_optim, fn = lssq_priors, method = "Nelder-Mead",
                                           control = list(maxit = p2ck["maxit_NelderMead"]))
     cat("optimized_par_Nelder_Mead\n",10^optimized_par_Nelder_Mead$par,"\n")

  # the first row of par_optim_series is a dummy and is removed here
     par_optim_series<-par_optim_series[-1,]
     rownames(par_optim_series) <- NULL

  # write the series of parameter estimates produced by the Nelder-Mead algorithm to an Excel file
  # the parameters have been written by the ssq_prior function in par_optim_series in nonlog values
     write.xlsx(x = par_optim_series,file = paste(".\\Results\\series_result_",filename,"_","_",
                          date_aux,".xlsx"),sheetName = "series",
                          append = FALSE,row.names = FALSE)

  # COPY OPTIMAL PARAMETERS FROM SERIES
  # find the result in the series with the minimum value for ssq in first column of par_optim_series
     ind <- which(par_optim_series[,1] == min(par_optim_series[,1],na.rm=TRUE))
     if (length(ind) > 1) {
       cat("\nWARNING There is more than 1 run with same minimum value of ssq - first one chosen")
       opt_par_optim <- par_optim_series[ind[1],]
     } else { opt_par_optim <- par_optim_series[ind,] }

     cat("\nIndex of optimal run ",ind,"\n")
     cat("\nOptimized parameter set\n",names(opt_par_optim),"\n",opt_par_optim,"\n")
     cat("\n\nThe full time course of the run with the lowest sqq is calculated and plotted")
   }

  # now plotting and Excel file writing is activated if it was not already in order to write final results
  p2ck["plot"] = 1. # set this to >0. to make plots
  p2ck["Excel_file"] = 1.  # set this to >0. to write Excel file for each run
  
  # for long plot to render sixth period in final plot
  if (p2ck["plot_long"] > 0.) { end_time_run <- p2ck["plot_long_endtime"] }

  # statement to initialize plot is repeated here for the case that only the end result is plotted
  if(p2ck["plot"] > 0.) {
    if (p2ck["plot_flux_pos_metabolites_neg"] >= 0) {
      plot(0,0,xlim=c(p2ck["tbegin"],end_time_run),
           xlab="Time (sec)", ylab="Net ATP synthesis (uM/s)",type="n",
           ylim=c(p2ck["y_axis_min_flux"],p2ck["y_axis_max_flux"]), col="blue")
           legend("topright", c("glycolysis", "oxidative phosphorylation"), fill=c("red", "blue"))
     }
    if (p2ck["plot_flux_pos_metabolites_neg"] < 0. ) {
      plot(0,0,xlim=c(p2ck["tbegin"],end_time_run),
           #xlab="Time (sec)", ylab="Concentration (microM)",type="n",
           xlab="Time (sec)", ylab="Content (micromol/ml cell)",type="n",
           #xlab="Time (sec)", ylab="Amount of Glucose (micromol/g wet weight)",type="n",
           ylim=c(p2ck["y_axis_min_met"],p2ck["y_axis_max_met"]), col="blue")
           legend("topright", c("FBP","Del lactate","5 * ATP"), fill=c("red","blue","green"))
     }
    }

  # setting array.addition to FALSE prevents that the parameters are written to par_optim_series by ssq_priors
  # during the recalculation of the optimal or calculation of the single value
  array.addition <- FALSE
  out1_times <- p2ck["tbegin"] : (p2ck["tbegin"]+p2ck["tduration1"])
if(p2ck["ensemble_run"] > 0.) {
  # get optimal parameters from the run
  par_optim2 <- opt_par_optim[(number_measures+1):length(opt_par_optim)]
  } else {
  # get the input parameters
      if(p2ck["load_opt_par_optim"] <= 0.) {
        par_optim2 <- 10^par_optim
      }
      else {
        load("opt_par_optim_saved")
        cat("\n\n WARNING: statement to assign loaded vector optimized values may need to be changed\n\n")
        par_optim2 <- opt_par_optim[(number_measures+1):length(opt_par_optim)] # 9:length(opt_par_optim [5:length(opt_par_optim)]
      }
}
# time points for long period with high sugar
  out6_times <- out6_times_output
  index_AdNsynt <- index_AdNsynt_out

  ssq_priors(log10(par_optim2))

  source("Rewrite Optimized Parameters.R")

  save.image(paste(".\\Results\\",filename,"_",date_aux,".RData",sep=""))

  title(paste("Final simulation ",date_aux," precise time ",date_aux1))

  cat("\nSave plot manually\n")

# Read data
data_model_high_glucose <- read.xlsx(file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),4)
data_exp_high_glucose <- dataCurve2 #read.xlsx(file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),21)
data_model_low_glucose <- read.xlsx(file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),2)
data_exp_low_glucose <- dataCurve #read.xlsx(file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),20)

source("Plot_Model_Exp_10panel.R")
plotA()

  #end writing R console text output to file
  sink()

  cat("\nEnd Program\n")