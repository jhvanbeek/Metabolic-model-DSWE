# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

# function to run the simulation over several desired periods and interventions
# the cost functions (sum of squares = ssq) are calculated
# provides the log target density (ltd) which is required if function spBayes is used
ssq_priors <- function(par_optim) {
                     date_aux1<<-format(Sys.time(),"_%H.%M.%S")
                     count <<- count + 1
                     cat("\ntime of this run = ",date_aux1,"  Count =",count,"\n")
                     names(par_optim) <- names_par
                     # the Nelder-Mead routine receives the log10 of the initial parameters from file
                     # and determines the next steps.
                     # It gives the parameters to the ssq_priors routine
                     # the model simulations require the nonlog values which are calculated here
                     nonlog_par_optim <- 10.^par_optim
                     names(nonlog_par_optim) <- names_par
                     # the present parameters that are evaluated are printed to the R console
                     print(nonlog_par_optim)
                     for (par_name in names_par) {
                            p2ck[par_name] <<- nonlog_par_optim[par_name]
                       } 
                     Vsynmax_penalty <- 0.
                     if (p2ck["Vsynmax_per_mg_mito"] < 1e-2) {
                       p2ck["Vsynmax_per_mg_mito"] <- 1e-2
                       Vsynmax_penalty <- 1e10
                      }
                                          
## SIMULATION FIRST PERIOD - THIS IS THE RUNIN FOR SECOND PERIOD

                     # initialize parameters and initial values for the first period
                     ssq_out2 = 0. 
  
                     p2ck["jump_to_glucose_per_ml_cell"] <<- p2ck["jump_to_glucose_per_ml_cell1"]
                     p2ck["vol_percent_cells"] <<- p2ck["vol_percent_cells1"]
                     p2ck["ATPconc"] <<- dataCurve[1,"ATP"] / p2ck["Vcyt"] * 1000  /
                                             ((1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     xstart["ATP"] <<- p2ck["ATPconc"]
                     p2ck <<- calculate_NAvalues_p2ck(p2ck,xstart)

  ## calculate initial values to calculate delta glc and delta lac, i.e. the change from t=0 for glucose and lactate, in model function
                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])

  # the concentration vector is calculated by integrating the ODE system over the first period
                     out1 <- as.data.frame(lsoda(xstart,out1_times, lmodel_runin, p2ck, hmax=p2ck["hmax"],
                                                  maxsteps=200000))
                     last_out1 <- length(out1[,1])
                     if(p2ck["Excel_file"]>0.) {
                       write.xlsx(x = out1, file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),
                                                                   sheetName = "out1", row.names = FALSE)
                      }
                     if(p2ck["plot"] > 0.) {plotOutLsoda(p2ck,out1)}

  ## the initial values for the second phase of the simulation are taken from the runin output
                     xinit <- as.numeric(out1[last_out1,2:(1+length(xstart))])
                     names(xinit)<-namesOutputVariables
  ## END OF THE FIRST PERIOD


  ## SIMULATION SECOND PERIOD - ADDITION LOW GLUCOSE CONCENTRATION

                     # initialize parameters and initial values for the second period
                     xinit["GLC"] <- p2ck["jump_to_Conc_glucose"]

                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])


                     # the concentration vector is calculated by integrating the ODE system over the second period
                     out2 <- as.data.frame(lsoda(xinit,out2_times, lmodel, p2ck, hmax=p2ck["hmax"],maxsteps=200000))

# calculate the sum of squares of deviations between model prediction and data
                     ssq_out2 = 0.
                     for (i in 1:length(out2_times)) {
                           ssq_out2 <- ssq_out2 + ((out2[i,"FBP_per_ml"] - dataCurve[i,"FBP"])/
                                                         dataCurve[i,"SD_FBP"])^2 +
                                          ((out2[i,"del_glc_per_ml"] -
                                               (dataCurve[i,"del_glucose_PCA_filtrates"] +
                                                dataCurve[i,"del_glucose_direct_suspension"])/2.)/
                                                   dataCurve[i,"SD_glc"])^2
                           if(!is.na(dataCurve[i,"ATP"])) {
                              ssq_out2 <- ssq_out2 + ((out2[i,"ATP_per_ml"] - dataCurve[i,"ATP"])/
                                                         dataCurve[i,"SD_ATP"])^2
                            }
                      }

                     if(p2ck["Excel_file"]>0.) {
                      write.xlsx(x = out2,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",
                                                         date_aux1,".xlsx"),
                                                  sheetName = "out2", append = TRUE, row.names = FALSE)
                      }

                     if(p2ck["plot"] > 0.) {plotOutLsoda(p2ck,out2)}
                     

  ## END OF THE SECOND PERIOD


## SIMULATION THIRD PERIOD - THIS IS THE RUNIN FOR FOURTH PERIOD

  ## calculate initial values to calculate delta glc and delta lac, i.e. the change from t=0 for glucose and lactate, in model function
                     p2ck["jump_to_glucose_per_ml_cell"] <<- p2ck["jump_to_glucose_per_ml_cell2"]
                     p2ck["vol_percent_cells"] <<- p2ck["vol_percent_cells2"]
                     p2ck <<- calculate_NAvalues_p2ck(p2ck,xstart)
                     xstart["ATP"] <<- ATPconc_start

                     # SIMULATION HIGH GLUCOSE DATA

  ## calculate initial values to calculate delta glc and delta lac in model functions
                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])


  # the concentration vector is calculated by integrating the ODE system over the third period
                     out3 <- as.data.frame(lsoda(xstart,out1_times, lmodel_runin, p2ck, hmax=p2ck["hmax"],
                                             maxsteps=200000))

                     last_out3 <- length(out3[,1])
                     if(p2ck["Excel_file"]>0.) {
                       write.xlsx(x = out3, file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),
                                                     sheetName = "out3", append = TRUE, row.names = FALSE)
                      }

                     if(p2ck["plot"] > 0.) {plotOutLsoda(p2ck,out3)}

  ## the initial values for the second phase of the simulation are taken from the runin output
                     xinit <- as.numeric(out3[last_out3,2:(1+length(xstart))])
                     names(xinit)<-namesOutputVariables
  ## END OF THE THIRD PERIOD


  ## SIMULATION FOURTH PERIOD - ADDITION HIGH GLUCOSE CONCENTRATION

                     xinit["GLC"] <- p2ck["jump_to_Conc_glucose"]

                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])


                     out4 <- as.data.frame(lsoda(xinit,out4_times, lmodel, p2ck, hmax=p2ck["hmax"],maxsteps=200000))
# calculate the sum of squares
                     ssq_out4 <- 0.
                     for (i in 1:length(out4_times)) {
                           ssq_out4 <- ssq_out4 + (( out4[i,"FBP_per_ml"] - dataCurve2[i,"FBP"])/
                                                     dataCurve2[i,"SD_FBP"])^2 +
                                   (( out4[i,"del_lac_per_ml"] - dataCurve2[i,"lactate"])/
                                                       dataCurve2[i,"SD_lac"])^2
                           if(!is.na(dataCurve2[i,"ATP"])) {
                              ssq_out4 <- ssq_out4 + ((out4[i,"ATP_per_ml"] - dataCurve2[i,"ATP"])/
                                                       dataCurve2[i,"SD_ATP"])^2
                            }
                           if(!is.na(dataCurve2[i,"glucose"])) {
                              ssq_out4 <- ssq_out4 + ((out4[i,"del_glc_per_ml"] - dataCurve2[i,"glucose"])/dataCurve2[i,"SD_glc"])^2
                            }
                           if(!is.na(dataCurve2[i,"VO2"])) {
                             ssq_out4 <- ssq_out4 + ((out4[i,"VO2"] - dataCurve2[i,"VO2"])/dataCurve2[i,"SD_VO2"])^2
                           }   
                     }
                     ssq_out4_pdp <- ssq_out4/31

                     # Penalties are applied if the present parameters during the optimization sequence
                     # are below the lower bound or above the upper bound set for the parameter range.
                     # This prevents that the parameters are given unrealistic values.
                     # In general these penalties will not add to the ssq value when this approaches 
                     # the optimal value,
                     # but they prevent the Nelder-Mead walk to stray off into unrealistic regions in
                     # parameter space where the integration algorithm may crash and local minima trap
                     # the optimization procedure
                     ssq_pars <- 0.
                     for (par_name in names_par) {
                       if (!is.na(lb[par_name]))   {
                           # par_optim and lb both contain the log10 values of the parameter
                           if (par_optim[par_name] < lb[par_name]) {
                                ssq_pars <- ssq_pars + ( (lb[par_name] - par_optim[par_name])/invweight_lb[par_name] )^2
                           }
                       }
                       if (!is.na(ub[par_name]))   {
                           # par_optim and ub both contain the log10 values of the parameter
                           if (par_optim[par_name] > ub[par_name]) {
                                ssq_pars <- ssq_pars +  ( (par_optim[par_name] - ub[par_name])/invweight_ub[par_name] )^2
                           }
                       }
                     }
                     ssq_pars <- ssq_pars + Vsynmax_penalty

                     if(p2ck["Excel_file"]>0.) {
                      write.xlsx(x = out4,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),
                                                              sheetName = "out4", append = TRUE, row.names = FALSE)
                      }

                     if(p2ck["plot"] > 0.) {plotOutLsoda(p2ck,out4)}

                     # the ssq_SD1 is calculated, which is the ssq with SDs set to 1
                     ssq_SD1 = 0.
                   if(p2ck["Excel_file"]>0.) {
                      for (i in 1:length(out2_times)) {
                          ssq_SD1 <- ssq_SD1 + ((out2[i,"FBP_per_ml"] - dataCurve[i,"FBP"]))^2 +
                                   (out2[i,"del_glc_per_ml"] - dataCurve[i,"del_glucose_PCA_filtrates"] )^2
                          if(!is.na(dataCurve[i,"ATP"])) {
                             ssq_SD1 <- ssq_SD1 + (out2[i,"ATP_per_ml"] - dataCurve[i,"ATP"])^2
                           }
                          }
                     }

                     ####  flush.console()  # works only under MS Windows, 
                                            # to see text output written to R console without delay
                     if(p2ck["Excel_file"]>0.) {
                       # calculate difference matrix model minus data
                        for (i in 1:length(out2_times)) {
                          # glucose decrease is counted positive in data and model
                          diffCurve[i,"lactate"] <- out2[i,"del_lac_per_ml"] - dataCurve[i,"lactate"]
                          diffCurve[i,"FBP"] <- out2[i,"FBP_per_ml"] - dataCurve[i,"FBP"]
                          diffCurve[i,"ATP"] <- out2[i,"ATP_per_ml"] - dataCurve[i,"ATP"]
                          diffCurve[i,"glucose"] <- out2[i,"del_glc_per_ml"] - (dataCurve[i,"del_glucose_PCA_filtrates"] + dataCurve[i,"del_glucose_direct_suspension"])/2.
                        }

                     for (i in 1:length(out4_times)) {
                          diffCurve2[i,"lactate"] <- out4[i,"del_lac_per_ml"] - dataCurve2[i,"lactate"]
                          diffCurve2[i,"FBP"] <- out4[i,"FBP_per_ml"] - dataCurve2[i,"FBP"]
                          diffCurve2[i,"ATP"] <- out4[i,"ATP_per_ml"] - dataCurve2[i,"ATP"]
                          diffCurve2[i,"glucose"] <- out4[i,"del_glc_per_ml"] - dataCurve2[i,"glucose"]
                          diffCurve2[i,"VO2"] <- out4[i,"VO2"] - dataCurve2[i,"VO2"]
                      }
                     }


#  SIMULATION LONG PERIOD - VERY HIGH GLUCOSE (11.1 mM)

                     p2ck["ATPconc"] <<- dataCurve2[1,"ATP"] / p2ck["Vcyt"] * 1000. / ((1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     p2ck["jump_to_glucose_per_ml_cell"] <<- 0.
                     p2ck["vol_percent_cells"] <<- 0.432
                     p2ck <<- calculate_NAvalues_p2ck(p2ck,xstart)
                     p2ck["jump_to_Conc_glucose"] <<- 11101.5

                     xstart["ATP"] <<- p2ck["ATPconc"]

  ## calculate initial values to calculate delta glc and delta lac in model functions
                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])

                     out5 <- as.data.frame(lsoda(xstart,out1_times, lmodel_runin, p2ck, hmax=p2ck["hmax"],
                                             maxsteps=200000))
                     last_out5 <- length(out5[,1])
                     if(p2ck["Excel_file"]>0.) {
                       write.xlsx(x = out5, file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),
                                            sheetName = "out5", append = TRUE, row.names = FALSE)
                      }

                     if(p2ck["plot"] > 0. && p2ck["plot_long"] > 0.) {plotOutLsoda(p2ck,out5)}

  ## the initial values for the second phase of the simulation are taken from the runin output
                     xinit <- as.numeric(out5[last_out5,2:(1+length(xstart))])
                     names(xinit)<-namesOutputVariables
  ## END OF THE SIMULATION OF THE INITIAL PERIOD WITHOUT GLUCOSE

                     xinit["GLC"] <- p2ck["jump_to_Conc_glucose"]

                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])


                     out6 <- as.data.frame(lsoda(xinit,out6_times, lmodel, p2ck, hmax=p2ck["hmax"],maxsteps=200000))

                     int_VO2 = 0. # integral of oxygen uptake
                     int_Vlac = 0. # integral of lactate production

                     # here a penalty for changes in FBP at late time points is imposed

                     FBPthreshold <- 0.5 * dataCurve2[13,"FBP"]
                     ssq_osc <- 0.
                     ssq_fbpthr <- 0.
                     ssq_atpthr <- 0.
                      for (i in 1:(length(out6_times)-1)) {
                           int_VO2 = int_VO2 + (out6[i+1,"time"] - out6[i,"time"])*(out6[i+1,"VO2"] + out6[i,"VO2"])/2.
                           int_Vlac = int_Vlac + (out6[i+1,"time"] - out6[i,"time"])*(out6[i+1,"Vlac"] + out6[i,"Vlac"])/2.
                           if (out6[i,"time"] > 600) {
                             ssq_osc <- ssq_osc + 10.*(out6[i-1,"FBP_per_ml"] - out6[i,"FBP_per_ml"])^2
                             if (out6[i,"ATP_per_ml"] > dataCurve2[1,"ATP"]) {
                               ssq_atpthr <- ssq_atpthr + (out6[i,"ATP_per_ml"] - dataCurve2[1,"ATP"])^2
                             }
                           if (out6[i,"FBP_per_ml"] < FBPthreshold) {
                             ssq_fbpthr <- ssq_fbpthr + (out6[i,"FBP_per_ml"] - FBPthreshold)^2
                             }                           
                           if (out6[i,"FBP_per_ml"] > dataCurve2[13,"FBP"]) {
                             ssq_fbpthr <- ssq_fbpthr + (out6[i,"FBP_per_ml"] - dataCurve2[13,"FBP"])^2
                           }                           
                           
                           }   
                     }
                     

                    #last_out6 <- length(out6[,1])

                     ssq_period_6 <- ssq_osc + ssq_atpthr + ssq_fbpthr 

                     av_VO2 = int_VO2/3600.
                     av_Vlac = int_Vlac/3600.
 
                     ssq_lac_VO2_rat <- (av_Vlac/av_VO2 - 3.)^2/0.2
                     ssq_out2_4_tot <- ssq_out4 + ssq_period_6 + ssq_pars + ssq_lac_VO2_rat + ssq_out2 # + ssq_Vstore


                     if(p2ck["plot"] > 0. && p2ck["plot_long"] > 0.) {plotOutLsoda(p2ck,out6)}

                     if(p2ck["Excel_file"]>0.) {
                        write.xlsx(x = out6,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),
                                       sheetName = "out6", append = TRUE, row.names = FALSE)
                     }


# Next a series of experiments is simulated where cells respiring on
# endogenous substrates were given different amounts of glucose at t=0 
# and the maximal FBP accumulation measured.
# In the simulations the maximal FBP concentration is calculated after giving a range of amounts of
# glucose.

## SIMULATION SEVENTH PERIOD - RUNIN FOR EIGTH PERIOD
                     ssq_FBP = NA #0.
               if(p2ck["Excel_file"]>0.) {
                     ssq_FBP = 0.
                  for (jj in 1:length(data3[,"cell_concentration"])) {
                     # initialize parameters and initial values for the first period
                     # here the maximal FBP concentration is calculated for a number
                     # of different amounts of added glucose in several experiments
                     p2ck["jump_to_glucose_per_ml_cell"] <<- data3[jj,"initial_glucose"]
                     p2ck["vol_percent_cells"] <<- data3[jj,"cell_concentration"]
                     p2ck["ATPconc"] <<- dataCurve[1,"ATP"] / p2ck["Vcyt"] * 1000  / ((1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     xstart["ATP"] <<- p2ck["ATPconc"]
                     p2ck <<- calculate_NAvalues_p2ck(p2ck,xstart)

  ## calculate initial values to calculate delta glc and delta lac, i.e. the change from t=0 for glucose and lactate, in model function
                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])

  # the concentration vector is calculated by integrating the ODE system over the first period
                     out7 <- as.data.frame(lsoda(xstart,out1_times, lmodel_runin, p2ck, hmax=p2ck["hmax"], maxsteps=200000))
                     last_out7 <- length(out7[,1])

  ## the initial values for the second phase of the simulation are taken from the runin output
                     xinit <- as.numeric(out7[last_out7,2:(1+length(xstart))])
                     names(xinit)<-namesOutputVariables
  ## END OF THE SEVENTH PERIOD


  ## SIMULATION EIGTH PERIOD - ADDITION LOW GLUCOSE CONCENTRATION

                     # initialize parameters and initial values for the second period
                     xinit["GLC"] <- p2ck["jump_to_Conc_glucose"]

                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])


                     # the concentration vector is calculated by integrating the ODE system over the second period
                     out8 <- as.data.frame(lsoda(xinit,out2_times, lmodel, p2ck, hmax=p2ck["hmax"],maxsteps=200000))

# calculate the sum of squares of deviations between model prediction and data
                     max_FBP <- max(out8[,"FBP_per_ml"])
                     diff_FBP <- max_FBP - data3[jj,"FBP_exp"]
                     ssq_FBP = ssq_FBP + diff_FBP^2

                     diffCurve3[jj,"cell_concentration"] <- data3[jj,"cell_concentration"]
                     diffCurve3[jj,"glucose_concentration"] <- data3[jj,"initial_glucose"]
                     diffCurve3[jj,"FBP_exp"] <- data3[jj,"FBP_exp"]
                     diffCurve3[jj,"FBP_model"] <- max_FBP
                     diffCurve3[jj,"difference"] <- diff_FBP

                     if(p2ck["Excel_file"]>0.) {
                        write.xlsx(x = out8,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),
                                     sheetName = paste("out8_",jj), append = TRUE, row.names = FALSE)
                        }

                      }
               }

## END OF THE EIGTH PERIOD
  
  #  Simulation of two experiments. Very low start concentrations of glucose were given at t=0 to 
  #  cells which respired on endogenous substrates before that time. The FBP 
  #  accumulated and the lactate produced after 5 and 10 sec was measured and calculated.
  ## SIMULATION NINTH PERIOD - RUNIN FOR TENTH PERIOD
                  ssq_out10 = NA #0.
                if(p2ck["Excel_file"]>0.) {
                  ssq_out10 = 0.
                  for (jj in 1:length(data4[,"cell_concentration"])) {
                     # initialize parameters and initial values for the first period
                     p2ck["jump_to_glucose_per_ml_cell"] <<- data4[jj,"initial_glucose"]
                     p2ck["vol_percent_cells"] <<- data4[jj,"cell_concentration"]
                     p2ck["ATPconc"] <<- dataCurve[1,"ATP"] / p2ck["Vcyt"] * 1000  / ((1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     xstart["ATP"] <<- p2ck["ATPconc"]

                     p2ck <<- calculate_NAvalues_p2ck(p2ck,xstart)

  ## calculate initial values to calculate delta glc and delta lac, i.e. the change from t=0 for glucose and lactate, in model function
                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xstart["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])

  # the concentration vector is calculated by integrating the ODE system over the first period
                     out9 <- as.data.frame(lsoda(xstart,out1_times, lmodel_runin, p2ck, hmax=p2ck["hmax"], maxsteps=200000))
                     last_out9 <- length(out9[,1])

  ## the initial values for the second phase of the simulation are taken from the runin output
                     xinit <- as.numeric(out9[last_out9,2:(1+length(xstart))])
                     names(xinit)<-namesOutputVariables
  ## END OF THE NINTH PERIOD


  ## SIMULATION TENTH PERIOD - ADDITION LOW GLUCOSE CONCENTRATION

                     # initialize parameters and initial values for the second period
                     xinit["GLC"] <- p2ck["jump_to_Conc_glucose"]

                     glc_per_ml_initial <<- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
                     lac_per_ml_initial <<- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) *
                                           xinit["LAC"])/ 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])


                     # the concentration vector is calculated by integrating the ODE system over the second period
                     out10 <- as.data.frame(lsoda(xinit,out4_times, lmodel, p2ck, hmax=p2ck["hmax"],maxsteps=200000))

                     ssq_out10 <- ssq_out10 + ( out10[2,"FBP_per_ml"] - data4[jj,"FBP_5sec"])^2. +
                                              ( out10[3,"FBP_per_ml"] - data4[jj,"FBP_10sec"])^2. +
                                              ( out10[2,"del_lac_per_ml"] - data4[jj,"lac_5sec"])^2. +
                                              ( out10[3,"del_lac_per_ml"] - data4[jj,"lac_10sec"])^2.
                     
                     diffCurve4[jj,"cell_concentration"] <- data4[jj,"cell_concentration"]
                     diffCurve4[jj,"glucose_concentration"] <- data4[jj,"initial_glucose"]
                     diffCurve4[jj,"lac_exp_5sec"] <- data4[jj,"lac_5sec"]
                     diffCurve4[jj,"lac_model_5sec"] <- out10[2,"del_lac_per_ml"]
                     diffCurve4[jj,"lac_diff_5sec"] <-  out10[2,"del_lac_per_ml"] - data4[jj,"lac_5sec"]
                     diffCurve4[jj,"lac_exp_10sec"] <-   data4[jj,"lac_10sec"]
                     diffCurve4[jj,"lac_model_10sec"] <- out10[3,"del_lac_per_ml"]
                     diffCurve4[jj,"lac_diff_10sec"] <-  out10[3,"del_lac_per_ml"] - data4[jj,"lac_10sec"]

                     diffCurve4[jj,"FBP_exp_5sec"] <- data4[jj,"FBP_5sec"]
                     diffCurve4[jj,"FBP_model_5sec"] <- out10[2,"FBP_per_ml"]
                     diffCurve4[jj,"FBP_diff_5sec"] <-  out10[2,"FBP_per_ml"] - data4[jj,"FBP_5sec"]
                     diffCurve4[jj,"FBP_exp_10sec"] <- data4[jj,"FBP_10sec"]
                     diffCurve4[jj,"FBP_model_10sec"] <- out10[3,"FBP_per_ml"]
                     diffCurve4[jj,"FBP_diff_10sec"] <-  out10[3,"FBP_per_ml"] - data4[jj,"FBP_10sec"]

                     if(p2ck["Excel_file"]>0.) {
                      write.xlsx(x = out10,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),
                                                  sheetName = paste("out10_",jj), append = TRUE, row.names = FALSE)
                      }

                     }
                }
  ## END OF THE TENTH PERIOD


                     ssq_FBP_pdp <- ssq_FBP/10. 
                     ssq_out10_pdp <- ssq_out10/4.
                     if(is.na(ssq_out2_4_tot)) ssq_out2_4_tot <- 100000000.

                     if (ssq_out2_4_tot <= min_ssq_out2_4) {
                        min_ssq_out2_4 <<- ssq_out2_4_tot
                     }

                     if(array.addition == TRUE) {    # variable array.addition is set in Main_NelderMead.R
                                                   # this statement makes it possible to use this function for the
                                                   # calculation of the full time course after the Nelder-Mead run is
                                                   # completed avoiding adding this run to par_optim_series
                     ssq_series <<- c(ssq_series,ssq_out2_4_tot)
                     par_optim_series <<- rbind(par_optim_series,c(ssq_out2_4_tot=ssq_out2_4_tot,
                                                                   min_ssq_out2_4_tot=min_ssq_out2_4,
                                                                   ssq_out4=ssq_out4,
                                                                   ssq_out4_pdp=ssq_out4_pdp,
                                                                   ssq_pars=ssq_pars, 
                                                                   ssq_osc=ssq_osc,
                                                                   ssq_atpthr=ssq_atpthr,
                                                                   ssq_fbpthr=ssq_fbpthr,
                                                                   ssq_period_6=ssq_period_6,
                                                                   ssq_out2=ssq_out2,
                                                                   ssq_lac_VO2_rat = ssq_lac_VO2_rat,
                                                                   av_VO2 = av_VO2,
                                                                   av_Vlac = av_Vlac,
                                                                   nonlog_par_optim))
                     }


                     if(p2ck["Excel_file"]>0.) {
                      write.xlsx(x = c(ssq_out4_pdp=ssq_out4_pdp,av_VO2 = av_VO2,av_Vlac = av_Vlac,
                                        p2ck,mitVsynmax=as.numeric(p2ck["Vsynmax"]),xstart),
                                     file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),
                                     sheetName = "parameters_startValues", append = TRUE, row.names = TRUE)
                      write.xlsx(x = dataCurve,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),sheetName = "dataCurve", append = TRUE, row.names = FALSE)
                      write.xlsx(x = dataCurve2,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),sheetName = "dataCurve2", append = TRUE, row.names = FALSE)
                      write.xlsx(x = diffCurve,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),sheetName = "diffCurve1", append = TRUE, row.names = FALSE)
                      write.xlsx(x = diffCurve2,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),sheetName = "diffCurve2", append = TRUE, row.names = FALSE)
                      write.xlsx(x = diffCurve3,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),sheetName = "diff_FBP", append = TRUE, row.names = FALSE)
                      write.xlsx(x = diffCurve4,file = paste(".\\Results\\out",filename,"_",date_aux,"_prec_",date_aux1,".xlsx"),sheetName = "diff_lact", append = TRUE, row.names = FALSE)
                      }
                     cat("ssq_out2_4_tot",ssq_out2_4_tot,
                         "min_ssq_out2_4_tot =",min_ssq_out2_4,
                         "s_osc",ssq_osc,
                         "ssq_fbpthr =",ssq_fbpthr,
                         "ssq_atpthr =",ssq_atpthr,
                         "\n",
                         "ssq_out4",ssq_out4,
                         "s_FBP =",ssq_FBP,
                         "s_p10 =",ssq_out10,
                         "ssq_pars =",ssq_pars,
                         "ssq_out4_pdp",ssq_out4_pdp,
                         "ssq_lac_VO2_rat",ssq_lac_VO2_rat,
                         "\n\n")
                     #a different version of ssq may be chosen
                     #ssq <- ssq_SD1
                     return(ssq_out2_4_tot)
                     }



