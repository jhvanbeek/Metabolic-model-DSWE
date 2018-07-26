# program and copyright: JHGM van Beek 2014-2018
# Nelder-Mead version optimization
# This code is licensed under the GNU General Public License v3.0
# Permissions of this strong copyleft license are conditioned on making available complete source code 
# of licensed works and modifications, which include larger works using a licensed work, under the same 
# license. Copyright and license notices must be preserved. Contributors provide an express grant of 
# patent rights.

To run the program all files (R scripts and Excel data files) should be in the R work directory.
The R packages loaded in Main_NelderMead.R should be installed in the R environment. References to the source of packages is given in the Main_NelderMead.R script.
The package which reads and writes Excel files depends on the Java programming language and works only under MS Windows.
Run the program by starting ("sourcing") Main_NelderMead.R which will call further R script files and read and write Excel workbooks.
The results are written to subfolder Results. This directory is made if it does not exist already.

By setting p2ck["ensemble_run"] <- 1. in the Run.R script a parameter optimization run can be done. This may take several hours. By setting p2ck["ensemble_run"] <- 0. a single run with a chosen parameter set can be done.

The program code as given here was tested with R version 3.4.1 in RStudio version 1.0.153 under Windows 10 and with R version 3.1.0 (32 bit version) not under RStudio.


