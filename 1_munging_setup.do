/*
Data Munging for
	Weakening of the cognition and height association from 1957 to 2018: 
	Findings from three British birth cohort studies 
Bann et al. (2021)

Cleans and derives variables for the 1946, 1958, and 1970 cohorts, separately.
Files are appended together in 2_mice.R, where variables are also renamed to account for different measurement ages
*/

/*
Set-Up
*/

clear
cls
set linesize 140

global bann 		"S:/IOEQSS_Main/Bann/crosscinequality"
global raw			"D:"
global clean 		"D:/Projects/Cognition and BMI/Data"
global code			"D:/Projects/Cognition and BMI/Code"
global cognition 	"D:/Projects/Cognition Measures/Data"

// ssc install egenmore
// ssc install wridit
// net install dm0004_1.pkg // zanthro


/*
Make Programmes
*/
capture program drop gen_bmi
program define gen_bmi
	syntax, Ages(string)
	
	foreach age of local ages{
	    gen bmi_`age' = weight_`age' / ( (height_`age' / 100) ^ 2)
	}
end

capture program drop gen_age
program define gen_age
	args new_var int_m int_y birth_my
	
	tempvar date
	gen `date' = ym(1900 + `int_y', `int_m') if `int_y' >= 0 & `int_m' >= 0
	gen `new_var' = (`date' - `birth_my') / 12
	drop `date'
end
	
capture program drop gen_residuals
program define gen_residuals
	syntax varlist [, covars(varlist)]
	
	foreach var of local varlist{
		local age = substr("`var'", -2, .)
		local stub = substr("`var'", 1, strpos("`var'", "_") - 1)
		
		regress `var' age_cog_`age' `covars'
		predict `stub'_resid_`age', residuals
	}
end

capture program drop label_class
program define label_class
	label define father_class 	///
		1 "V Unskilled" 2 "IV Partly Skilled" ///
		3 "III Skilled Manual" 4 "III Skilled Non-Manual" ///
		5 "II Intermediate" 6 "I Professional"
	label values father_class* father_class
end

capture program drop closer_wide
program define closer_wide
	args file id
	use `id' visitage wtself htself bmi ///
	using "${raw}/CLOSER/`file'_closer_wp1", clear
	tostring visitage, gen(age)
	replace age = "0" + age if visitage < 10
	drop visitage
	gen self_report_ = max(wtself, htself)
	label values self_report_ wtself
	rename bmi bmi_
	keep `id' age *_
	reshape wide *_, i(`id') j(age) string
end


/*
Run Do Files
*/
do "${code}/1a_munging_nshd.do"
do "${code}/1b_munging_ncds.do"
do "${code}/1c_munging_bcs70.do"