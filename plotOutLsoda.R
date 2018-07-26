# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

plotOutLsoda <- function(p2ck,out) {
   if (p2ck["plot_lines"] > 0.) {
        if (p2ck["plot_flux_pos_metabolites_neg"] >= 0) {
           lines(out[,"time"], out[,"VglycATPnet"], col="red", lwd=3)
           lines(out[,"time"], out[,"Vmit"], col="blue", lwd=3)
             }

        if (p2ck["plot_flux_pos_metabolites_neg"] < 0) {
           lines(out[,"time"], out[,"FBP_per_ml"], col="red", lwd=3)
           lines(out[,"time"], out[,"del_lac_per_ml"], col="blue", lwd=3)
           lines(out[,"time"], out[,"ATP_per_ml"]*5., col="green", lwd=3)           
        }

   } else {
           for (i in 1:length(out[,1])) {
             if (p2ck["plot_flux_pos_metabolites_neg"] >= 0) {
               points(out[i,"time"], out[i,"VglycATPnet"], pch=20, col="red")
               points(out[i,"time"], out[i,"Vmit"], pch=20, col="blue")
             }
             if (p2ck["plot_flux_pos_metabolites_neg"] < 0) {
               points(out[i,"time"], out[i,"FBP_per_ml"], pch=20, col="red")
               points(out[i,"time"], out[i,"del_lac_per_ml"], pch=20, col="blue")
               points(out[i,"time"], out[i,"ATP_per_ml"]*5., pch=20, col="green")
             }
           }
     }
  }