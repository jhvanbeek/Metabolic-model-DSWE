# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

# Plotting function
plotA = function(){
  pdf(file = paste("Results\\Figure2_sim_",date_aux,".pdf",sep=''), width = 10, height = 10)
  par(mfrow = c(5,2), oma = c(6.0, 8.0, 3.0, 3.0), mai = c(0.1, 0.1, 0.4, 0.4)) # mai = c(0.5, 0.5, 0.6, 0.3))
  CEX <- 1.
  par(cex.axis = 1.6,cex.lab = 1, cex.main=1.5)   # cex.main=1.6
  
  plot(data_exp_low_glucose[,'time'], data_exp_low_glucose[,'FBP'], xlim=c(0,300), ylim=c(0,3.8),
       cex = 2, pch = 20, xaxp=c(0,300,3),xlab='',ylab='')
  lines(data_model_low_glucose[,'time'], data_model_low_glucose[,'FBP_per_ml'],type="l")
  title("Low glucose experiment",line=1.0,cex.main=2.)
  mtext(expression(paste('Fructose-1,6-bisphoshate', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 4.8, adj=0.5)
  mtext(expression(paste('(', mu, 'mol/ml cell)', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 2.65, adj=0.5)
  plot(data_exp_high_glucose[,'time'], data_exp_high_glucose[,'FBP'], xlim=c(0,300), ylim=c(0,3.8),
       cex = 2, pch = 20, xaxp=c(0,300,3),xlab='',ylab='')
  lines(data_model_high_glucose[,'time'], data_model_high_glucose[,'FBP_per_ml'],type="l")
  title("High glucose experiment",line=1.0,cex.main=2.)
  plot(data_exp_low_glucose[,'time'], data_exp_low_glucose[,'ATP'], xlim=c(0,300), ylim=c(0,3.8),
       cex = 2, pch = 20, xaxp=c(0,300,3),xlab='',ylab='')
  lines(data_model_low_glucose[,'time'], data_model_low_glucose[,'ATP_per_ml'],type="l")
  mtext(expression(paste('ATP', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 4.8, adj=0.5)
  mtext(expression(paste('(', mu, 'mol/ml cell)', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 2.65, adj=0.5)
  plot(data_exp_high_glucose[,'time'], data_exp_high_glucose[,'ATP'], xlim=c(0,300), ylim=c(0,3.8),
       cex = 2, pch = 20, xaxp=c(0,300,3),xlab='',ylab='')
  lines(data_model_high_glucose[,'time'], data_model_high_glucose[,'ATP_per_ml'],type="l")
  
  # here follow glucose and lactate changes
  plot(data_exp_low_glucose[,'time'], (data_exp_low_glucose[,'del_glucose_PCA_filtrates']+
                                         data_exp_low_glucose[,'del_glucose_direct_suspension'])/2., 
       xlim=c(0,300), ylim=c(0,12),
       cex = 2, pch = 20, xaxp=c(0,300,3),yaxp=c(0,12,3),xlab='',ylab='')
  lines(data_model_low_glucose[,'time'], data_model_low_glucose[,'del_glc_per_ml'],type="l")
  mtext(expression(paste('Glucose taken up', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 4.8, adj=0.5)
  mtext(expression(paste('(', mu, 'mol/ml cell)', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 2.65, adj=0.5)
  
  plot(data_exp_high_glucose[,'time'], data_exp_high_glucose[,'glucose'], xlim=c(0,300), ylim=c(0,12.),
       cex = 2, pch = 20, xaxp=c(0,300,3),yaxp=c(0,12,3),xlab='',ylab='')
  lines(data_model_high_glucose[,'time'], data_model_high_glucose[,'del_glc_per_ml'],type="l")
  plot(NA, NA, xlim=c(0,300), ylim=c(0,12.),
       cex = 2, pch = 20, xaxp=c(0,300,3),yaxp=c(0,12,3),xlab='',ylab='')
  lines(data_model_low_glucose[,'time'], data_model_low_glucose[,'del_lac_per_ml'],type="l")
  mtext(expression(paste('lactate produced', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 4.8, adj=0.5)
  mtext(expression(paste('(', mu, 'mol/ml cell)', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 2.65, adj=0.5)
  plot(data_exp_high_glucose[,'time'], data_exp_high_glucose[,'lactate'], xlim=c(0,300), ylim=c(0,12.),
       cex = 2, pch = 20, xaxp=c(0,300,3),yaxp=c(0,12,3),xlab='',ylab='')
  lines(data_model_high_glucose[,'time'], data_model_high_glucose[,'del_lac_per_ml'],type="l")
  
  plot(NA, NA, xlim=c(0,300), ylim=c(0,50.),
       cex = 2, pch = 20, xaxp=c(0,300,3),yaxp=c(0,50,5),xlab='',ylab='')
  lines(data_model_low_glucose[,'time'], data_model_low_glucose[,'VO2'],type="l")
  mtext(expression(paste('Oxygen consumption', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 4.8, adj=0.5)
  mtext(expression(paste('(', mu, 'M/s)', sep = '')),
        side = 2, outer = FALSE, cex = CEX, cex.lab=2, cex.axis=1., cex.main=1., cex.sub=1.,
        line = 2.65, adj=0.5)
  plot(data_exp_high_glucose[,'time'], data_exp_high_glucose[,'VO2'], xlim=c(0,300), ylim=c(0,50.),
       cex = 2, pch = 20, xaxp=c(0,300,3),yaxp=c(0,50,5),xlab='',ylab='')
  lines(data_model_high_glucose[,'time'], data_model_high_glucose[,'VO2'],type="l")
  
  
  mtext(expression(paste('Time (sec)', sep = '')), side = 1, outer = TRUE, cex = 2.0, line = 4.2)
  
  dev.off()
}