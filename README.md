# Informativeness of Estimation Moments

This is the replication package for the paper "The Informativeness of Estimation Moments"
by Bo E. Honoré, Thomas H. Jørgensen and Àureo de Paula to appear in Journal of Applied Econometrics.
Working paper version can be found here: https://www.cemmap.ac.uk/publication/id/14668

The replication package has three main folders: STATA, Matlab and Python.

STATA:
------
The STATA-folder contains all the .do-files used to construct the empirical measures from the BHSP. All BHPS data can be downloaded free of charge at https://www.iser.essex.ac.uk/bhps/acquiring-the-data

Running the "main.do" file produces all output required. In that file the path to the raw BHPS data should be specified along with the path of the replication package.

Matlab:
------
The following Matlab files produce the structural empirical estimation results reported in the paper. All results are stored in the "Results" folder.
<tt>`r sort(d$product_name)`</tt>
main.m:			main file. Running this file will run estimation for both educational groups and run the post-estimation files to produce the output.

run_estimate.m:		estimation file. You can run this file by it's own.	

run_post_estimation.m:	post-estimation production of tables and figures. Loads the estimation results from the run_estimate.m file	

fun.m:			static class used to call all functions to estimate the empirical retirement model

mex_combi.mexw64: 	C++ compiled file calculating the values for all combinations of retirement alternatives

The following folders/files are used to estimate the two examples in the paper:

ex1_Probit:		folder containing the Matlab files for generating the results for this example. Run the "main_probit.m" file.

ex2_Weibul:		folder containing the Matlab files for generating the results for this example. Run the "main_weibul.m" file.

Python:
------
Probit.ipynb:		Python notebook that replicates the first (Probit) example in the paper. The exact numbers differ slighty due to different random numbers drawn in thiw implementation relative to the Matlab implementation.
