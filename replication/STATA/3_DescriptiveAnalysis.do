clear all
set type double
set scheme s1mono
ssc install estout /*install this package if not on machine*/

set more off
use STATA\BHPSdata\BHPSselectedHouseholds, replace

label var cohort1 "Birth year, husband"
label var cohort2 "Birth year, wife"


global SAMPLE ExpSample /*Use here the sample where both spouses are NOT retired*/

* Share of households in which both are working:
sum ExpSample if MatlabSample==1

* table with descriptive statistics  
global VARS age1 age2 ExpRetAge1 ExpRetAge2 ExpRetYearDiff Highskilled1 Highskilled2 gp_many1 gp_many2 ExpWorseHealth1 ExpWorseHealth2 Income1 Income2 pensPPP1 pensPPP2  pensEPS1 pensEPS2
eststo desc: estpost sum $VARS if ${SAMPLE}==1
estout desc using "STATA\figures\DescriptiveTable.tex", style(tex) replace varlabel(age1 "\cmidrule(lr){2-6} Age, husband") cells("mean(fmt(3) label(Mean)) sd(fmt(2) label(Std.)) min(fmt(0) label(Min)) max(fmt(0) label(Max)) count(fmt(0) label(Obs.))") mlabels(none) label 

* plot Expected joint retirement
hist ExpRetYearDiff if ${SAMPLE}==1,disc name(JointHist,replace) graphregion(color(white)) color(gs11) lcolor(black) xlabel(-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30)
graph export "STATA\figures\JointHist.eps",replace
* look at age-differences
hist ExpRetYearDiff if ${SAMPLE}==1 & AgeDiff>=2,disc name(JointHist_Agediff2,replace) graphregion(color(white)) color(gs11) lcolor(black) xlabel(-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30)
graph export "STATA\figures\JointHist_agediff2.eps",replace

cap drop AgeDiff_neg w
g AgeDiff_neg = - AgeDiff
egen w = group(AgeDiff ExpRetYearDiff) if  ${SAMPLE}==1
twoway (scatter  ExpRetYearDiff AgeDiff if ${SAMPLE}==1) (line AgeDiff AgeDiff_neg if ${SAMPLE}==1), name(RetDiff_scatter,replace)

* Reform plot of SPA age of husband and wife
preserve
keep if ${SAMPLE}==1
collapse SPA* , by(cohortm2)
replace cohortm2 = cohortm2/12
replace SPAm2 = SPAm2/12
replace SPAm1 = SPAm1/12
twoway (line SPAm1 cohortm2 , lcolor(gs9) lp("-") lwidth(.5)) (line SPAm2 cohortm2 , lcolor(black) lwidth(.5)), ///
			legend(on order(1 "Men" 2 "Women") ring(0) pos(5) col(1)) name(SPA_combined,replace) xtitle("Birth year, wife") ytitle("State pension age (SPA)") graphregion(color(white)) ///
			ylabel(60(1)65) xlabel(1940(5)1970)
restore
graph export "STATA\figures\SPA_combined.eps",replace


