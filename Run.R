# here parameters more specific for calculation methods, plotting and timing of the different phases
# of the simulation are given

# NOTE: parameters which precede function definitions have already been calculated with default 
# parameters.
# in particular this applies to calculation of forward and backward max rates of enzymes based on dry-wet ratio etc. etc.

# adding parameters for calculation
p2ck["hmax"] <- 60.
p2ck["ensemble_run"] <- 1. # set to 1. for ensemble run
p2ck["load_opt_par_optim"] <- -1. # if > 0. and if there is no optimization run, parameters will be loaded from
                                  # object opt_par_optim from disk, but only if ensemble_run=0
                                  # parameters from Excel file will then be ignored
p2ck["maxit_NelderMead"] <- 5 #23400. # low number to test, high number for optimization

p2ck["index_time_AdNsynt"] <- 1800

# adding parameters for timing of simulations
p2ck["tstep"] <- 1.
p2ck["plotinterval"] <- 1.# gives plot interval and writing to console, but no effect on Excel output files
p2ck["tbegin"] <- -60.
p2ck["tduration1"] <- -p2ck["tbegin"] # determines runin time
p2ck["tduration2"] <- 300. # during second period. after runin
p2ck["pruning"] <- 1 #e.g., 10 : every tenth step of the ODE solver is written to Excel file

# adding parameters for plotting
p2ck["plot"] <- 0. # set this to >0., e.g. 1., to make plots
p2ck["plot_long"] <- 1. # set this to >0. if plot should include fifth and sixth period, i.e. long simulation
p2ck["plot_long_endtime"] <- 3600.
p2ck["plot_lines"] <- 0. # set this to >0 to plot lines, otherwise points are plotted
p2ck["plot_flux_pos_metabolites_neg"] <- -1 # >0 plot flux, <0 plot metabolites
p2ck["Excel_file"] <- 0.  # set this to >0. to write Excel file for each run

# parameters for plot axes
p2ck["y_axis_min_flux"] <- -200.
p2ck["y_axis_max_flux"] <- 200.
p2ck["y_axis_min_met"] <- 0.
p2ck["y_axis_max_met"] <- 16.

