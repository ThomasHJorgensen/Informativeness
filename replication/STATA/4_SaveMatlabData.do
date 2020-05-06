clear all
set more off
use STATA\BHPSdata\BHPSselectedHouseholds, replace

global SAMPLE 		MatlabSample /*Use the sample where retirees are included*/
global MatlabData 	age ExpRetAge SPA gp_many White Highskilled pensPPP pensEPS lWealth Income retired ExpWorseHealth cohort /*list of variables to be saved and used in Matlab code*/

******************************************************
* Save matlab .txt data files
******************************************************
quietly {
	forvalues member = 1/2 {
		foreach var in $MatlabData {
			replace `var'`member' = -2 if missing(`var'`member') /* Prepare data: Replace missing with -2 such that importation of data in Matlab works */
			outfile `var'`member' if ${SAMPLE}==1 using "Matlab\data\\`var'`member'.txt", comma wide replace  /*Each variable is saved in a.txt file with individual information in the collumn*/

		}
	}
}

