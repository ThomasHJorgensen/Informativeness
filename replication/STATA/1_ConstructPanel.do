* This do-file merges relevant datasets within each wave, stores the relevant varibles and subsequently combines all waves
clear all
set more off
set memory 2g

* Specify which variables to keep (indivual and hh) 
global VARSi 	lvlong hlsf10c dobm doby hl2gp pppen jbpen jbpenm agexrt retex ageret hid bwtag* doiy4 doim pno hid hgr2r plbornc hsownd hscost race yr2uk4 nmar nchild lprnt lnprnt fiyr* xpfood ncars carown isced fenow feend  lchnmor pay* bankk* debt* nvest* save* svack* sex mlstat jbstat age  nchild
global VARShh   nkids 

local wave = 0
foreach PRE in a b c d e f g h i j k l m n o p q r {
local wave = `wave' + 1
use "${BHPS_path}\\`PRE'indresp.dta", replace
qui destring, replace

**********************************
* Merge with household variables *
* Sort by household id (hid)
sort `PRE'hid
tempfile Individual
save `Individual'
use "${BHPS_path}\\`PRE'hhresp.dta",replace
sort `PRE'hid 
merge `PRE'hid using `Individual'

* Remove prefix and add wave variable 
renpfix `PRE'
gen wave = `wave'
* Rename the pid (now id in 2006)
if `wave'==16 {
	ren id pid
}

***************************************
*To save RAM: only keep some variables.
* The rather clumsy way of doing this is due to the fact that some variables are only in the sample for some years
loc VarsToKeep
foreach v in pid wave ${VARSi} ${VARShh} {
	loc j    	
	cap unab j: `v'
		loc VarsToKeep  `VarsToKeep' `j' 
}
keep `VarsToKeep'
***************************************

* Append with previous data
if `wave'>1 {
	append using `AllWaves', force
}

tempfile AllWaves
save `AllWaves'

}  /* End of year loop*/

sort pid wave
drop if missing(pid)


*************************************************************************************
* Use bracketed wealth information:  1: impute to lowest value 2: impute to midpoint *
* Recode many variables to missing if negative (inapplicable)
foreach var of varlist nkids save* bankk* nvest* debt* svack*{
	replace `var' = . if `var'<0 
}
* Find the "highest" value in the brackets. Since they are ordered 1000, 5000, 10000, 500, I rename the last bracket to zero rather than 4 and loop from 0-3
cap ren bankkb4 bankkb0
cap ren savekb4 savekb0
cap ren nvestc4 nvestc0 
cap ren debtc4 debtc0
cap ren svackb4 svackb0
* For some variables, there is five groups in 2005
cap ren nvestc5 nvestc4 
cap ren debtc5 debtc4
cap ren svackb5 svackb4
cap drop Bracket
g Bracket = .
cap drop BracketNv
g BracketNv = .
cap drop BracketSa
g BracketSa = .
cap drop BracketDe
g BracketDe = .
cap drop BracketSv
g BracketSv = .
forvalues br = 0/3 {
	replace Bracket 	= `br' if bankkb`br'==1 & !missing(bankkb`br')
	replace BracketSa 	= `br' if savekb`br'==1 & !missing(savekb`br')
}
forvalues br = 0/4 {
	replace BracketNv 	= `br' if nvestc`br'==1 & !missing(nvestc`br')
	replace BracketDe 	= `br' if debtc`br'==1 & !missing(debtc`br')
	replace BracketSv 	= `br' if svackb`br'==1 & !missing(svackb`br')
}
* If "no" to bracket 0-> then let bracket be -1 and use that below
replace Bracket 	= -1 if bankkb0==2 & !missing(bankkb0)
replace BracketNv 	= -1 if nvestc0==2 & !missing(nvestc0)
replace BracketSa 	= -1 if savekb0==2 & !missing(savekb0)
replace BracketDe 	= -1 if debtc0==2 & !missing(debtc0)
replace BracketSv 	= -1 if svackb0==2 & !missing(svackb0)

* Construct the imputed measures
* First: use the bottom (midpoint in lower bracket (-1)
cap drop IMPbankk
g IMPbankk = bankk
replace IMPbankk = (Bracket==-1)*250 + (Bracket==0)*500 + (Bracket==1)*1000 + (Bracket==2)*5000 + (Bracket==3)*10000 if missing(IMPbankk)
replace IMPbankk = . if missing(bankk) & missing(Bracket)
cap drop IMPnvestk
g IMPnvestk = nvestk /*nvestk: ALL investments, both sole and joint*/
replace IMPnvestk = (BracketNv==-1)*500 + (BracketNv==0)*1000 + (BracketNv==1)*5000 + (BracketNv==2)*15000 + (BracketNv==3)*50000 + (BracketNv==4)*100000 if missing(IMPnvestk)
replace IMPnvestk = . if missing(nvestk) & missing(BracketNv)
cap drop IMPsavek
g IMPsavek = savek
replace IMPsavek = (BracketSa==-1)*250 + (BracketSa==0)*500 + (BracketSa==1)*1000 + (BracketSa==2)*5000 + (BracketSa==3)*10000 if missing(IMPsavek)
replace IMPsavek = . if missing(savek) & missing(BracketSa)
cap drop IMPdebty
g IMPdebty = debty
replace IMPdebty = (BracketDe==-1)*50 + (BracketDe==0)*100 + (BracketDe==1)*500 + (BracketDe==2)*1500 + (BracketDe==3)*5000 + (BracketDe==4)*10000 if missing(IMPdebty)
replace IMPdebty = . if missing(debty) & missing(BracketDe)
cap drop IMPsvack
g IMPsvack = svack
replace IMPsvack = (BracketSv==-1)*250 + (BracketSv==0)*500 + (BracketSv==1)*1000 + (BracketSv==2)*5000 + (BracketSv==3)*10000 + (BracketSv==4)*20000 if missing(IMPsvack)
replace IMPsvack = . if missing(svack) & missing(BracketSv)

* Second: use the midpoint between brackets (and put in an arbitrary point in the end..
cap drop IMP2bankk
g IMP2bankk = bankk
replace IMP2bankk = (Bracket==-1)*250 + (Bracket==0)*(500+1000)/2 + (Bracket==1)*(1000+5000)/2 + (Bracket==2)*(5000+10000)/2 + (Bracket==3)*(10000+2500) if missing(IMP2bankk)
replace IMP2bankk = . if missing(bankk) & missing(Bracket)
cap drop IMP2nvestk
g IMP2nvestk = nvestk 
replace IMP2nvestk = (BracketNv==-1)*500 + (BracketNv==0)*(1000+5000)/2 + (BracketNv==1)*(5000+15000)/2 + (BracketNv==2)*(15000+50000)/2 + (BracketNv==3)*(50000+100000)/2 + (BracketNv==4)*(100000+25000) if missing(IMP2nvestk)
replace IMP2nvestk = . if missing(nvestk) & missing(BracketNv)
cap drop IMP2savek
g IMP2savek = savek
replace IMP2savek = (BracketSa==-1)*250 + (BracketSa==0)*(500+1000)/2 + (BracketSa==1)*(1000+5000)/2 + (BracketSa==2)*(5000+10000)/2 + (BracketSa==3)*(10000+2500) if missing(IMP2savek)
replace IMP2savek = . if missing(savek) & missing(BracketSa)
cap drop IMP2debty
g IMP2debty = debty
replace IMP2debty = (BracketDe==-1)*50 + (BracketDe==0)*(100+500)/2 + (BracketDe==1)*(500+1500)/2 + (BracketDe==2)*(1500+5000)/2 + (BracketDe==3)*(5000+10000)/2 + (BracketDe==4)*(10000+2500) if missing(IMP2debty)
replace IMP2debty = . if missing(debty) & missing(BracketDe)
cap drop IMP2svack
g IMP2svack = svack
replace IMP2svack = (BracketSv==-1)*250 + (BracketSv==0)*(500+1000)/2 + (BracketSv==1)*(1000+5000)/2 + (BracketSv==2)*(5000+10000)/2 + (BracketSv==3)*(10000+20000)/2 + (BracketSv==4)*(20000+5000) if missing(IMP2svack)
replace IMP2svack = . if missing(svack) & missing(BracketSv)

* Delete variables used to impute wealth-information
drop  nvest* debt* savek* bankk* Bracket* svack*

*************************************************************************************
sort pid wave
save STATA\BHPSdata\BHPSpanel, replace


