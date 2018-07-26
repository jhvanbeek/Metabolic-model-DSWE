# DEFAULT VALUES OF THE PARAMETER ARRAY
# initialise the model parameter array p2ck with default values
# these can be modified before simulation to establish the desired test condition

# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

p2ck_default  <-  c(
      	        neg_flag = 0.,
      	        #concentrations for conserved moieties
      	        # total NADH + NAD concentration is 200 microM
                # here the NADH concentration is given. NAD is a variable in the differential equations and
                # NADH is calculated as NADHtot - NAD
                NADHconc = 0.13,
      	        #initial ATP concentration, calculated later from experimental data
                # intracellular Pi concentration near mitochondria
                Pi_mit = 6000.,
                O2 = 150.,
                # fractional volumes of cytosol, intermembrane space and mitochondrial matrix
      		      Vcyt = 0.90,
      		      Vmat = 0.09,  # mitochondrial matrix volume fraction
                Vims = 0.01,
      		      k_hyd_ATP	= 0.118173,
                Vint = NA,
                FBPmultiplier = NA,
                                     # This is the ratio of total phosphorylated glycolytic intermediates
                                     # to fructose-1,6-bisphosphate.
                                     # to take into account that the total phosporylated carbohydrate pool which
                                     # is greater than the FBP pool which is measured
                                     # other metabolites are G6P, F6P, DHAP, glyceraldehyde-3-P and 
                                     # even beyond GAPDH
            		wet_weight_cells_mg_per_ml = NA, #117., #134., #167., #250., #167,
		            vol_percent_cells = 0.432, #adapted later to correspond to simulated experiments
		            vol_percent_cells1 = 2.6, # adapted later to correspond to simulated experiments
                vol_percent_cells2 = 2.9, # adapted later to correspond to simulated experiments

                # Dry to wet weight ratio
                DW_ratio = 0.18, # dry wet weight ratio suspension for Ehrlich ascites tumour cells

                # protein weight to dry weight ratio
                PD_ratio = 0.75,

                # ATP per O2 ratio
                PO2_rat = 5.6, # given by Coe, E.L. Cancer Research 26:269-275, 1966

                # density tissue 1.06 g wet weight/ml
                dens = 1.06,

                Vmax_head_per_mg = 0.003797685, # will be estimated from data
                Vmax_tail_per_mg = 0.3 * 0.46/60., 
                Vldh_max_f_per_mg = 2.0/60.,
                Vsynmax_per_mg_mito = 0.024144,
          
                k_store = 0.01,

		# ATP synthesis
		            Kadp_mit = 25.0, # microM, as in Van Beek 2007
		            Kpi_mit = 800.0, # microM, as in Van Beek 2007; Beard, 2005, Fig. 2B
                Kmit_pyr = 0.1, # Halestrap, Biochemical Journal 148:85-96, 1975 : this dependency prevents negative pyruvate concentrations
                Kmit_O2 = 0.26, # Froese BBA 1962
                Kmit_nadh = 0.001,

                ## ATP hydrolysis
                Khyd_ATP_ADP = 0.0001, # ATP hydrolysis is made low if ATP in cytosol falls to micromolar values
                kADPbreakd = 0.0001,
                kADPsynthesis = 5.e-5,

		# Total contents

                # parameters head end glycolysis
                Katp_head = 1000.,
                Kglc_head = 47.,

    # parameters tail end glycolysis : glyceraldehyde phosphate + Pi + NAD + 2ADP = pyruvate + NADH + 2ATP
		# defined in Excel file Optimized_parameters.xslx

                # lactate dehydrogenase (LDH) parameters
                Kldh_pyr = 335., # microM
                Kldh_nadh = 2., 
                Kldh_lac = 17000.,
                Kldh_nad = 849.,
                KeqLDH = 16198,

                jump_to_Conc_glucose_total_susp =  25., # adapted later to correspond with simulated experiments
                jump_to_Conc_glucose_total_susp2 = 4000., # adapted later to correspond with simulated experiments
                jump_to_glucose_per_ml_cell = 3.5,
                jump_to_glucose_per_ml_cell1 = 3.5, # micromol/ml cells
                jump_to_glucose_per_ml_cell2 = 26.6 # micromol/ml cells

                )
            # end of initialisation p2ck array
