# FUNCTION TO CALCULATE Not Available VARIABLES FOR PARAMETER ARRAY 'p2ck' WHICH DEPEND ON VALUES OTHER VARIABLES
# additional parameters are calculated for the p2ck array
# these were set to 'NA' values in the initialisation of the p2ck array

# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

calculate_NAvalues_p2ck <- function(p2ck,xstart) {

                           p2ck["AdNtot"] = xstart["ADP"] + p2ck["ATPconc"]
                           p2ck["NADtot"] = xstart["NAD"] + p2ck["NADHconc"]

                           # calculate Vmax of head and tail end glycolysis in microM/sec, 
                           # i.e. micromole per litre cell water per sec
                           # divide by 60 sec/min; start with value per mg protein;
                           # calculate based on one gram wet weight; multiply with mg protein in one gram wet weight;
                           # divide by ml of water per gram wet weight, taken for the whole cell
                           # it is assumed that density water is 1.0 g/ml
                           p2ck["Vmax_head"] <- p2ck["Vmax_head_per_mg"] * p2ck["Vcyt"] * p2ck["DW_ratio"] * p2ck["PD_ratio"] *
                                                                       1000. / ((1-p2ck["DW_ratio"]) * 1e-3)
                           p2ck["Vmax_tail"] <- p2ck["Vmax_tail_per_mg"] * p2ck["Vcyt"] * p2ck["DW_ratio"] * p2ck["PD_ratio"] *
                                              1000. / ((1-p2ck["DW_ratio"]) * 1e-3)
                           p2ck["Vldh_max_f"] <- p2ck["Vldh_max_f_per_mg"] * p2ck["Vcyt"] * p2ck["DW_ratio"] * p2ck["PD_ratio"] *
                                                                       1000. / ((1-p2ck["DW_ratio"]) * 1e-3)
                           # Vsynmax, maximal rate of mitochondrial ATP synthesis in microM / s
                           p2ck["Vsynmax"] <- p2ck["Vsynmax_per_mg_mito"] * p2ck["Vmat"] * p2ck["DW_ratio"] * p2ck["PD_ratio"] *
                                                                       1000. / ((1-p2ck["DW_ratio"]) * 1e-3)
                           if (!is.na(as.numeric(p2ck["wet_weight_cells_mg_per_ml"]))) {
                             ww <- p2ck["wet_weight_cells_mg_per_ml"]/1000.
                             p2ck["Vint"] <- (1.- ww/p2ck["dens"])/ ((1.-p2ck["DW_ratio"]) * ww)
                             }

                           if (!is.na(as.numeric(p2ck["vol_percent_cells"]))) {
                             ww <- p2ck["vol_percent_cells"]*10.*p2ck["dens"]/1000.
                             p2ck["Vint"] <- (1.- ww/p2ck["dens"])/ ((1.-p2ck["DW_ratio"]) * ww)
                             }

                           p2ck["jump_to_Conc_glucose"] <- p2ck["jump_to_glucose_per_ml_cell"]*1000./(p2ck["Vint"] + p2ck["Vcyt"])/
                                                                                         ((1.-p2ck["DW_ratio"]) * p2ck["dens"])
                           p2ck
                           }
