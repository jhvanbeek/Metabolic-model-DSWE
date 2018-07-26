# MODEL FUNCTION
# model function to be integrated by lsoda routine
# returns the time derivatives of state variables plus calculated results for diffusive and reactive fluxes

# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

# concentrations are given in micromole per liter of water in the compartment
# for ADP, ATP, NAD and FBP this means the water volume in the cytosolic space is used
# for glucose, lactate and pyruvate the sum of the water volume in the cytosol and extracellular space is used

# reaction and exchange rates are given in micromole per sec per liter of total intracellular water

model <- function(t, y, p2ck) {
  ydot <- vector()

   ####  flush.console()  # works only under MS Windows. Speeds up display of text to R console, but costs lots of CPU time

  ybelowzero <- y < 0.
  sumbelowzero <- sum(ybelowzero)
  if (sumbelowzero > 0) {
            y[ybelowzero] <- 0.
            p2ck["neg_flag"] <<- 1
     }

  y_calc_NADH <- as.numeric(p2ck["NADtot"] - y["NAD"])


#CALCULATED FLUXES - function definitions are given in file RM_function_definitions.R
  decrease_AdN <- xstart["ATP"] + xstart["ADP"] - y["ATP"] - y["ADP"]
  Vhyd <- as.numeric(lVhyd(t,decrease_AdN = decrease_AdN,                              
                         k_hyd_AdN=p2ck["k_hyd_AdN"],interc=p2ck["interc_ATP_hyd"]))
  Vmit_array <- as.numeric(lVmit(ADP=y["ADP"],Pi=p2ck["Pi_mit"],PYR=y["PYR"],NADH=y_calc_NADH,O2=p2ck["O2"],Kadp=p2ck["Kadp_mit"],
                                  Kpi=p2ck["Kpi_mit"],Kmit_pyr=p2ck["Kmit_pyr"],Kmit_nadh=p2ck["Kmit_nadh"],
                                  Kmit_O2=p2ck["Kmit_O2"],Vsynmax=p2ck["Vsynmax"],
                                  PO2_rat=p2ck["PO2_rat"],uncoupler=uncoupler))
  Vmit <- as.numeric(Vmit_array[1])
  VO2  <- as.numeric(Vmit_array[2])
  Vpyr_mit <- as.numeric(Vmit_array[3])
  Vnadh_mit <- as.numeric(Vmit_array[4])
  fPYR_mit <- as.numeric(Vmit_array[5])
  fNADH_mit <- as.numeric(Vmit_array[6])

  Vlac <- as.numeric(lVlac(PYR=y["PYR"],NADH=y_calc_NADH,LAC=y["LAC"],NAD=y["NAD"],Kldh_pyr=p2ck["Kldh_pyr"],
                      Kldh_nadh=p2ck["Kldh_nadh"],Kldh_lac=p2ck["Kldh_lac"],Kldh_nad=p2ck["Kldh_nad"],
                      KeqLDH=p2ck["KeqLDH"],Vldh_max_f=p2ck["Vldh_max_f"]))

  Vhead <- as.numeric(lVhead(ATP=y["ATP"],GLC=y["GLC"],Fhead_unb = y["Fhead_unb"], 
                           Katp=p2ck["Katp_head"],
                           Kglc=p2ck["Kglc_head"],Vmax=p2ck["Vmax_head"]))

  Vtail <- as.numeric(lVtail(FBP=y["FBP"],NAD=y["NAD"],ADP=y["ADP"],NADH=y_calc_NADH,Kfbp=p2ck["Kfbp_tail"],
                       Knad=p2ck["Knad_tail"],Kadp=p2ck["Kadp"],Ki_tail_nadh=p2ck["Ki_tail_NADH"],
                       Vmax=p2ck["Vmax_tail"]))

  if (y["store_glc"] > 6000.) store_inh <- 1000./(1000. + (y["store_glc"] - 6000.))
  else store_inh <- 1.
  Vstore_for <- p2ck["k_store"] * y["FBP"]* store_inh 
  Vstore <- Vstore_for

  AdN_tot <- (y["ATP"] + y["ADP"])
  AdN_start <- xstart["ATP"] + xstart["ADP"]
  if ( AdN_tot < AdN_start ) {
             Vadp_synth <- p2ck["kADPsynthesis"] * (AdN_start - AdN_tot)
      } else Vadp_synth <- 0.

# adenine nucleotide breakdown.
   VadpBreakdown <- as.numeric(p2ck["kADPbreakd"]*y["ADP"]^2)
   glyc_ATP_net <- 2.*Vtail - 2.*Vhead # net glycolytic ATP production  
                                       # NOTE: Coefficient for Vtail is 2: for each glyceraldehyde-3-phosphate
                                       # metabolized two ATP molecules are synthesized
   Vbal <- glyc_ATP_net + Vmit - Vhyd

#concentration derivatives
  # Vint is the interstitial or extracellular space volume
  # Vcyt is the cytosolic space volume
  # glucose (GLC) is distributed in Vint and Vcyt
  ydot["GLC"] <-  -Vhead / (p2ck["Vint"] + p2ck["Vcyt"])
  # lactate (LAC) is distributed in Vint and Vcyt
  ydot["LAC"] <-  Vlac / (p2ck["Vint"] + p2ck["Vcyt"])
  ydot["PYR"] <-  (Vtail - Vlac - Vpyr_mit ) / (p2ck["Vint"] + p2ck["Vcyt"]) # mitochondria take up pyruvate,
  ydot["ADP"] <- (- Vbal - VadpBreakdown + Vadp_synth)/ (p2ck["Vcyt"])
  ydot["ATP"] <- Vbal / p2ck["Vcyt"]
  # p2ck["FBPmultiplier"]: the FBP pool in the model also represents other phosphorylated intermediates. The total of FBP
  # plus the other phosphorylated glycolytic intermediates is taken into account in the fluxes.
  # Here the strict FBP pool size is calculated because this is measured experimentally.
  ydot["FBP"] <- (Vhead - 0.5 * Vtail - Vstore) / p2ck["FBPmultiplier"] / p2ck["Vcyt"]
  # Fhead_unb is the fraction of head enzyme that is unbound with FBP
  # FBP inhibition of enzyme in head section of glycolysis has been reported at high Mg concentration by Gumaa & McLean BBRC 1969
  # 2 mM F16BP inhibited by 90% at Mg2+ = 7 mM
  # this inhibition represents negative feedback on the head section exerted by glycolytic intermediates such
  # as glucose-6-phosphate and fructose-1,6-bisphosphate
  ydot["Fhead_unb"] <- - p2ck["kFBP_f"]* y["FBP"]*y["Fhead_unb"] + p2ck["kFBP_b"]*(1.-y["Fhead_unb"])
  ydot["NAD"] <- (-Vtail + Vlac + Vnadh_mit) / p2ck["Vcyt"]

  # the store_glc is calculated as moles per liter cytosolic water
  ydot["store_glc"] <- Vstore  / p2ck["Vcyt"]

  if (sum(is.na(ydot))>0){
    cat("\n\nAt least one ydot is NA, meaning it cannot be calculated")
    cat("\n\nydot = ",t(ydot))
    recover()
  }

  # calculating the content in micromol per ml cell volume (composed of total of water and dry components)
  FBP_per_ml <- as.numeric( (p2ck["Vcyt"] * y["FBP"]) / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
  glc_per_ml <- as.numeric( (p2ck["Vint"] + p2ck["Vcyt"]) * y["GLC"] / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
  lac_per_ml <- as.numeric( ( (p2ck["Vint"] + p2ck["Vcyt"]) * y["LAC"]) / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
  store_glc_per_ml <- as.numeric( (p2ck["Vcyt"] * y["store_glc"]) / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
  pyr_mit_per_ml_glc_unit <- as.numeric( (p2ck["Vcyt"] * 0.5 * y["pyr_mit"]) / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
  # change from initial values at start of current phase
  del_glc_per_ml <-  glc_per_ml_initial - glc_per_ml
  del_lac_per_ml <-  lac_per_ml - lac_per_ml_initial
  ATP_per_ml <- as.numeric( (p2ck["Vcyt"] * y["ATP"]) / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])
  ADP_per_ml <- as.numeric( (p2ck["Vcyt"] * y["ADP"]) / 1000 * (1.-p2ck["DW_ratio"]) * p2ck["dens"])

  functionValues <- c(
    NADH = y_calc_NADH,
    Vhyd = Vhyd,
    Vmit = Vmit,
    Vhead = Vhead,
    Vtail = Vtail,
    Vlac = as.numeric(Vlac),
    Vstore = as.numeric(Vstore),
    VglycATPnet = as.numeric(glyc_ATP_net),
    VadpBreakdown = VadpBreakdown,
    Vadp_synth = Vadp_synth,
    V_ATPtot = as.numeric(Vmit + glyc_ATP_net),
    VO2 = VO2,
    Vpyr_mit = as.numeric(Vpyr_mit),
    Vnadh_mit = as.numeric(Vnadh_mit),
    glc_per_ml = glc_per_ml,
    lac_per_ml = lac_per_ml,
    FBP_per_ml = FBP_per_ml,
    ATP_per_ml = ATP_per_ml,
    ADP_per_ml = ADP_per_ml,
    ATP_plus_ADP = ATP_per_ml + ADP_per_ml,
    del_lac_per_ml = del_lac_per_ml,
    del_glc_per_ml = del_glc_per_ml,
    store_glc_per_ml = store_glc_per_ml,
    Vstore_for = as.numeric(Vstore_for)
)

      list(ydot,functionValues)
  }