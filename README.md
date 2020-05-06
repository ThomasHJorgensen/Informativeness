# Informativeness of Estimation Moments

This is the replication package for the paper "The Informativeness of Estimation Moments"
by Bo E. Honoré, Thomas H. Jørgensen and Àureo de Paula to appear in Journal of Applied Econometrics.
Working paper version can be found here: https://www.cemmap.ac.uk/publication/id/14668

The replication package has three main folders: STATA, Matlab and Python.

STATA:
------
The STATA-folder contains all the .do-files used to construct the empirical measures from the BHSP. All BHPS data can be downloaded free of charge at https://www.iser.essex.ac.uk/bhps/acquiring-the-data

Running the <tt>`main.do</tt> file produces all output required. In that file the path to the raw BHPS data should be specified along with the path of the replication package.

Matlab:
------
The following Matlab files produce the structural empirical estimation results reported in the paper. All results are stored in the "Results" folder.
* <tt>main.m:</tt>			main file. Running this file will run estimation for both educational groups and run the post-estimation files to produce the output.
* <tt>run_estimate.m</tt>:		estimation file. You can run this file by it's own.	
* <tt>run_post_estimation.m</tt>:	post-estimation production of tables and figures. Loads the estimation results from the run_estimate.m file	
* <tt>fun.m</tt>:			static class used to call all functions to estimate the empirical retirement model
* <tt>mex_combi.mexw64</tt>: 	C++ compiled file calculating the values for all combinations of retirement alternatives

The following folders/files are used to estimate the two examples in the paper:
* <tt>ex1_Probit</tt>:		folder containing the Matlab files for generating the results for this example. Run the "main_probit.m" file.
* <tt>ex2_Weibul</tt>:		folder containing the Matlab files for generating the results for this example. Run the "main_weibul.m" file.

Python:
------
* <tt>Probit.ipynb</tt>:		Python notebook that replicates the first (Probit) example in the paper. The exact numbers differ slighty due to different random numbers drawn in thiw implementation relative to the Matlab implementation.
