# FUNCTIONS TO CALCULATE REACTION AND TRANSPORT RATES
# the rate functions to calculate enzyme and transport fluxes are given below.
# these are used to calculate the time derivatives necessary to integrate the array of
# concentration variables in function 'model'.
# equation for Vmit is described in Van Beek JHGM, Am J Physiol 2007

# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

# the reaction or diffusion rates are always given per unit of volume of total intracellular water

# functions are compiled here using R function cmpfun

# Rate function for head of glycolysis
Vhead <- function(ATP,GLC,Fhead_unb,Katp,Kglc,Vmax) { 
  den_head <- 1 + ATP/Katp + GLC/Kglc + (ATP*GLC)/(Katp*Kglc)
  Vhead <- Vmax * Fhead_unb * ATP/Katp * GLC/Kglc / den_head 
  if(is.na(Vhead)) browser()
  Vhead
}
lVhead <- cmpfun(Vhead)


Vtail <- function(FBP,NAD,ADP,NADH,Kfbp,Knad,Kadp,Vmax,Ki_tail_nadh){
             a_fbp <- FBP/(FBP+Kfbp)
             a_nad <- NAD/(NAD+Knad)
             a_adp <- ADP/(ADP+Kadp)
             i_nadh <- 1/(1+NADH/Ki_tail_nadh)
             Vtail <- Vmax * a_nad * a_adp * a_fbp * i_nadh
             Vtail
             }
lVtail <- cmpfun(Vtail)


Vmit <- function(ADP,Pi,PYR,NADH,O2,Kadp,Kpi,Kmit_pyr,Kmit_nadh,Kmit_O2,Vsynmax,PO2_rat,uncoupler) { 
	DenSyn <- 1+ ADP/Kadp + Pi/Kpi + ADP*Pi/(Kadp*Kpi)
	fPYR <- PYR/(PYR + Kmit_pyr)
	fNADH <- NADH/(NADH + Kmit_nadh)
        sum_aux <- 5*fPYR + fNADH # weighted for 10 electrons per pyruvate, 2 electrons per NADH
	fCOM <-  sum_aux/6.
	Vmit <- Vsynmax*ADP*Pi / (Kpi*Kadp*DenSyn) * fCOM * O2/(O2 + Kmit_O2)
        VO2 <- Vmit/PO2_rat
        if (uncoupler > 0.) {
          Vmit <- -0.2*Vsynmax
          VO2 <- (Vsynmax/PO2_rat)* fCOM * O2/(O2 + Kmit_O2)
           }
        Vmit_pyr <- 5*fPYR/sum_aux * VO2/2.5 # 10 electrons per pyruvate, 4 electrons per O2
        Vmit_nadh <- fNADH/sum_aux * VO2 * 2 # 2 electrons per NADH, 4 electrons per O2
        c(Vmit,VO2,Vmit_pyr,Vmit_nadh,fPYR,fNADH)
	}
lVmit <- cmpfun(Vmit)

Vlac <- function(PYR,NADH,LAC,NAD,Kldh_pyr,Kldh_nadh,Kldh_lac,Kldh_nad,KeqLDH,Vldh_max_f) {
                        # Lactate dehydrogenase from Lambeth and Kushmerick
                        # positive flux means net lactate formation
        DenSyn <- 1+PYR/Kldh_pyr+NADH/Kldh_nadh+(PYR/Kldh_pyr)*(NADH/Kldh_nadh)+LAC/Kldh_lac+NAD/Kldh_nad+
                                               (LAC/Kldh_lac)*(NAD/Kldh_nad)
        result <- Vldh_max_f * (PYR/Kldh_pyr)*(NADH/Kldh_nadh) - Vldh_max_f * Kldh_lac * Kldh_nad /(Kldh_pyr*Kldh_nadh*KeqLDH)*
                                               (LAC/Kldh_lac)*(NAD/Kldh_nad)
        result <- result/DenSyn
        return(result)
        }
lVlac <- cmpfun(Vlac)


# function to calculate ATP hydrolysis in cytosol.
Vhyd <- function(t,decrease_AdN,k_hyd_AdN,interc) {
  outVhyd = interc - k_hyd_AdN * decrease_AdN
  if(decrease_AdN < 0.) outVhyd = interc
  if (outVhyd < 0) outVhyd <- 0.  
  return(outVhyd)
}
lVhyd <- cmpfun(Vhyd)