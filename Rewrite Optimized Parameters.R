# Script to overwrite initial optimizable parameters with final optimized parameters and save Excel file
# overwrite optpas with new values in par_optim2 

# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

optpars_new <- optpars
names_paroptim <- names(par_optim2)
for (i in names_paroptim) {
  for (j in 1:length(optpars[,"par_optim"])) {
      if(optpars_new[j,"names"]==i) optpars_new[j,"par_optim"] <- par_optim2[i] 
  } 
}

# write to new Optimized Parameters.xls file
write.xlsx(x = optpars_new,file = paste(".\\Results\\Optimized_parameters",filename,"_",date_aux,".xlsx"), append = FALSE, row.names = FALSE)
