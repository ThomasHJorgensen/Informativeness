* This do file saves a data set on the household level
clear all

set more off
set memory 1g

global MinAge 40
global MaxAge 59


use STATA\BHPSdata\BHPSpanel, replace
xtset pid wave

global vars "lvlong hlsf10c dobm doby hl2gp pppen jbpen jbpenm age jbstat retex ageret agexrt doiy4 doim fiyrl Wealth lWealth mlstat isced race nkids"
cap drop Wealth
egen Wealth 	= rowtotal(IMP2bankk IMP2savek IMP2svack) /*Savings are based on several (mutually exclusive) variables*/
replace Wealth = . if missing(IMP2bankk) & missing(IMP2savek) & missing(IMP2svack)
g lWealth = l.Wealth

********************************************************************************************
* 1. construct household data set                                                             *
* *****************************************************************************************
keep if hgr2r==1|hgr2r==2|hgr2r==3  /* only keep reference person and spouse */

gen nn=1
bysort hid wave: egen nnhh=sum(nn)

tab nnhh                                                                        

keep if nnhh==2  /* only keep obs where there are two spouses in the household */

bysort pid: egen sex_temp = max(sex)
replace sex = sex_temp 	/* Let households with missing sex information in some years have the same sex as in others*/
drop sex_temp

sort hid wave sex

by hid wave: egen nsex=sum(sex)
tab nsex
keep if nsex==3 /* delete same sex couple */


**********************************************************************************
* 2. Construct household identifier                                             *
**********************************************************************************
keep $vars hid wave sex pid
reshape wide $vars pid, i(hid wave) j(sex)  

sort pid2 pid1 wave
by pid2 pid1: gen PwS =_n /*running from 1 within each couple*/

by pid2: g temp = 1 if PwS==1
g HID = sum(temp) /*construct the new household identifier based on these "new" couples*/

* Let the new HID be the household identifier
drop hid
ren HID hid
order hid, first

sort  hid wave
xtset hid wave


**********************************************************************************
* 3. Construct new individual variables 
**********************************************************************************
foreach var of varlist lvlong* hlsf10c* ageret* retex* fiyrl* agexrt* pppen* jbpenm* race* hl2gp* age* doim* doiy4* dobm* doby*{
	replace `var' = . if `var'<=0 /*missings are recoded from -9 to .*/
}

forvalues sex=1/2 {

* a. Retirement age and expectations
cap rename agexrt`sex' RetAge`sex'
cap drop retired`sex'
g retired`sex' = inlist(jbstat`sex',4) & !missing(jbstat`sex')

g ExpRetAge`sex' 	   = ageret`sex'
replace ExpRetAge`sex' = retex`sex' if missing(ExpRetAge`sex')
replace ExpRetAge`sex' = RetAge`sex' if retired`sex'==1 & !missing(retired`sex')  /*insert actual retirement age, if retired*/
* Update retirement to include negative durations
replace retired`sex' = 1 if ExpRetAge`sex'<age`sex' & !missing(ExpRetAge`sex') /* if expected retirement age is less than actual age, we take them as being retired*/ 

g ExpRetYear`sex' = doiy4`sex' + (ExpRetAge`sex'-age`sex')

* b. Financial variables, income, wealth and pension schemes
g Income`sex' 		 = fiyrl`sex'/1000 /*Labor market income*/
replace lWealth`sex' = lWealth`sex'/1000

g pensPPP`sex' = 0
replace pensPPP`sex' = 1 if pppen`sex' == 1 & !missing(pppen`sex')
g pensEPS`sex' = 0
replace pensEPS`sex' = 1 if jbpenm`sex' == 1 & !missing(jbpenm`sex')

* c. Background characteristics
* Educational level
g Highskilled`sex' 		 = 0
replace Highskilled`sex' = 1 if isced`sex'>=6 /*High skilled is "first degree" first or second stage of tertiary education (ISCED 5 and 6: http://www.uis.unesco.org/Library/Documents/isced97-en.pdf)*/
replace Highskilled`sex' = . if missing(isced`sex')

* RACE: White. The race is determined in the first wave and then only updated for first-time interviews
cap drop temp
g temp = race`sex'==1
replace temp = . if missing(race`sex') 
sort hid wave
cap drop White`sex'
by hid: egen White`sex' = max(temp)

* d. Health-related
*GP visits, many: 10+, moderate:3-10
g gp_many`sex' = 0
replace gp_many`sex' = 1 if hl2gp`sex'==5 
replace gp_many`sex' = . if missing(hl2gp`sex')

cap drop gp_moderate`sex'
g gp_moderate`sex' = 0
replace gp_moderate`sex' = 1 if inlist(hl2gp`sex',3,4) 
replace gp_moderate`sex' = . if missing(hl2gp`sex')

*Expectations about health
cap drop ExpWorseHealth`sex'
g ExpWorseHealth`sex' 		= l2.hlsf10c`sex' < 3 /*if "def. sure" or "mostly sure" that health will worsen 2 waves prior (when they were asked) */
replace ExpWorseHealth`sex' = . if missing(l2.hlsf10c`sex')

rename doby`sex' cohort`sex'
}
cap drop ageret* retex* pid* isced* temp race* hl2gp* lvlong* hlsf10c*

*************************
* State Pension Age (SPA)
g cohortm2 = (cohort2*12+dobm2)
g cohortm1 = (cohort1*12+dobm1) 
cap drop SPAm2
g SPAm2 = cohortm2 - (1950*12+3)
replace SPAm2 = 0 if SPAm2<0
replace SPAm2 = 5*12 if SPAm2>5*12
replace SPAm2 = 60*12 + SPAm2
g SPAm1       = 65*12

*replace cohort2 = 1949 if cohort2==1950 & dobm2<=3 /*shift women into the cohort 1949 if born before april 1950. This is just to ease the construction of the treatment*/

g SPA1 	= 65
g SPA2  = 60 + (cohort2-1950)
replace SPA2 = 60 if SPA2<60
replace SPA2 = 65 if SPA2>65

*g Cohort5054 = inlist(cohort2,1950,1951,1952,1953,1954) /*Thomas changed Agusut '17: includes some with very small change in their SPA*/
g Cohort5054 = inlist(cohort2,1951,1952,1953,1954)
g Cohort5558 = cohort2>1954 & !missing(cohort2)
g treatment  = 0
replace treatment = 1 if Cohort5054==1 & !missing(Cohort5054)
replace treatment = 2 if Cohort5558==1 & !missing(Cohort5558)

***************************
* Household level variables
g ExpRetYearDiff = ExpRetYear1 - ExpRetYear2
g ExpRetYearDiff_work = ExpRetYearDiff
replace ExpRetYearDiff_work = . if (retired1==1 | retired2==1)
g AgeDiff 	  	 = age1 - age2

cap drop Married
g Married = 0
replace Married = 1 if mlstat1 == 1 | mlstat2 == 1
replace Married = . if missing(mlstat1)
cap drop mlstat*

g nkids = nkids2
drop nkids1 nkids2

* Asign row-total to generate household variables                                     
foreach var in Wealth lWealth Income{
	egen `var' = rowtotal(`var'1 `var'2)
	replace `var' = . if missing(`var'1) & missing(`var'2)
}


********************************************************************************
* Restrict sample
**********************************************************************************
g remove = 0
forvalues sex=1/2 {
	replace remove = 1 if (ExpRetAge`sex'<50 | ExpRetAge`sex' > 70) & !missing(ExpRetAge`sex')  /*Planned retirement in 50,70*/
}

g MatlabSample 	 = remove==0 & (!missing(ExpRetAge1) & !missing(ExpRetAge2)) & (retired1==0 | retired2==0)  & (age1>=${MinAge} & age1<=${MaxAge} & age2>=${MinAge} & age2<=${MaxAge})  /* must have expectations for both spouses*/ 
*g MatlabSample 	 = remove==0 & (!missing(ExpRetAge1) & !missing(ExpRetAge2)) & (!inlist(jbstat1,4) | !inlist(jbstat2,4))  & (age1>=${MinAge} & age1<=${MaxAge} & age2>=${MinAge} & age2<=${MaxAge})  /* must have expectations for both spouses*/ 
g ExpSample 	 = MatlabSample==1 & retired1==0 & retired2==0 /*both has to be working*/

cap drop _merge PwS
drop pppen* doim* doiy* jbstat* jbpen* dobm* fiyrl*

label var nkids "Number of children at home"
label var age1 "Age, husband"
label var age2 "Age, wife"
label var White1 "White, husband"
label var White2 "White, wife"
label var ExpRetAge1 "Planned retirement age, husband"
label var ExpRetAge2 "Planned retirement age, wife"
label var Highskilled1 "High skilled, husband"
label var Highskilled2 "High skilled, wife"
label var pensPPP1 "Private pension, husband"
label var pensPPP2 "Private pension, wife"
label var pensEPS1 "Employer pension, husband"
label var pensEPS2 "Employer pension, wife"
label var Income1 "Labor income (£1,000), husband"
label var Income2 "Labor income (£1,000), wife"
label var lWealth "Household savings (£1,000)"
label var ExpRetYearDiff "Diff. in planned retirement year (husband-wife)"
label var ExpRetYearDiff_work "Diff. in planned retirement year (husband-wife), both work"
label var ExpRetAge1 "Planned retirement age, husband"
label var ExpRetAge2 "Planned retirement age, wife"
label var gp_moderate1 "3-10 GP visits, husband"
label var gp_moderate2 "3-10 GP visits, wife"
label var gp_many1 "10+ GP visits, husband"
label var gp_many2 "10+ GP visits, wife"
label var ExpWorseHealth1 "Expect worse health, husband"
label var ExpWorseHealth2 "Expect worse health, wife"
label var SPA1 "State Pension Age (SPA), husband"
label var SPA2 "State Pension Age (SPA), wife"
label var Cohort5054 "Wife born in 1951-1954"
label var Cohort5558 "Wife born after 1954"
label var treatment "Treatment, 0-2"
label var cohort1 "cohort, husband"
label var cohort2 "cohort, wife"

label define treatment_v 0 "wife born before 1951" 1 "wife born in 1951-1954" 2 "wife born after 1954"
label values treatment treatment_v  


save STATA\BHPSdata\BHPSselectedHouseholds, replace
