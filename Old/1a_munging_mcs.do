/*
2001 Cohort (MCS)
*/
global sweeps 0 3 5 7 11 14 17

capture program drop fix_cnum
program define fix_cnum
	qui foreach x in A B C D E F G{
		replace CNUM00 = `x'CNUM00 if !missing(`x'CNUM00) & missing(CNUM00)
	}
	qui replace CNUM00 = 1 if missing(CNUM00)
	qui foreach x in A B C D E F G{
		replace `x'CNUM00 = CNUM00 if missing(`x'CNUM00) & !missing(CNUM00)
	}
end

closer_wide mcs mcsid
save "${clean}/2001c_temp.dta", replace  

use MCSID COUNTRY GOVWT2 NOCMHH WEIGHTGB using ///
	"${raw}/MCS/xwave/mcs_longitudinal_family_file.dta", clear
expand 2, gen(XX)
gen CNUM00 = XX + 1
drop XX
foreach x in A B C D E F G{
	gen `x'CNUM00 = CNUM00
}

merge m:1 MCSID using "${bann}/mcs/0y/mcs1_derived_variables.dta", ///
	nogen keepusing(ADFINH00 ADMHGT00 APHGTM00 ADMBMI00 APDBMIA0)
fix_cnum

forval i = 1/4{
	local j: word `i' of `c(ALPHA)'
	local s: word `i' of ${sweeps}
	
	merge 1:1 MCSID `j'CNUM00 using "${raw}/mcs/`s'y/mcs`i'_cm_derived.dta", ///
		nogen keepusing(`j'DC06E00)
	fix_cnum
}

forval i = 1/7{
	local j: word `i' of `c(ALPHA)'
	local s: word `i' of ${sweeps}
	
	preserve
		tempfile sex
		use MCSID `j'CNUM00 *CSEX00* using ///
			"${raw}/mcs/`s'y/mcs`i'_hhgrid.dta", clear
		drop if !inrange(`j'CNUM00, 1, 2)
		save `sex', replace
	restore
	
	merge 1:1 MCSID `j'CNUM00 using "`sex'", nogen
	fix_cnum
}

merge 1:1 MCSID DCNUM00 using "${raw}/mcs/7y/mcs4_cm_interview.dta", ///
	nogen keepusing(DCMCS4AGE DCHTCM00 DCWTCM00)
fix_cnum

merge 1:1 MCSID ECNUM00 using "${raw}/mcs/11y/mcs5_cm_interview.dta", ///
	nogen keepusing(ECHTCMA0 ECWTCMA0 EMCS5AGE)
fix_cnum

merge 1:1 MCSID ECNUM00 using "${raw}/mcs/11y/mcs5_cm_derived.dta", ///
	nogen keepusing(EVSABIL EDCE0600)
fix_cnum

merge 1:1 MCSID FCNUM00 using "${raw}/mcs/14y/mcs6_cm_derived.dta", ///
	nogen keepusing(FCWRDSC FCMCS6AG FDCE0600)
fix_cnum

merge 1:1 MCSID FCNUM00 using "${raw}/mcs/14y/mcs6_cm_interview.dta", ///
	nogen keepusing(FCHTCM00 FCWTCM00)	
fix_cnum

merge 1:1 MCSID GCNUM00 using "${raw}/mcs/17y/mcs7_cm_interview.dta", ///
	nogen keepusing(GCHTCM00 GCWTCM00)	
fix_cnum

merge 1:1 MCSID GCNUM00 using "${raw}/mcs/17y/mcs7_cm_derived.dta", ///
	nogen keepusing(GCMCS7AG)	
fix_cnum

rename MCSID mcsid
merge m:1 mcsid using "${bann}/derived/closer_education/mcs_ParentEducationDerived.dta", ///
	nogen keepusing(mmled mfled mcs_mledg mcs_fledg)
fix_cnum

preserve
	tempfile ses
	use mcsid csep using "${bann}/mcs/derived/mcs_bmi_ses_closer v001.dta" ///
		if !missing(csep), clear
	replace mcsid = subinstr(mcsid, ".", "", .)
	save "`ses'", replace
restore
merge m:1 mcsid using "`ses'", nogen	// 261 new participants! Presumably MCSID has changed.
fix_cnum

merge m:1 mcsid using "${clean}/2001c_temp.dta", // nogen


// Sex
recode *CSEX00* (min/-1 = .) (3/max = .)
gen male = .
foreach var of varlist *CSEX00* {
	replace male = 2 - `var' if !missing(`var') & missing(male)
}


// BMI
gen bmi_17 = GCWTCM00 / ( (GCHTCM00 / 100) ^ 2) if GCWTCM00 > 0 & GCHTCM00 > 0
gen self_report_17 = 2 if !missing(bmi_17)
label values self_report_17 wtself


// Ethnic Group
local ethnic_vars	ADC06E00 BDC06E00 CDC06E00 DDC06E00 EDCE0600 FDCE0600
recode `ethnic_vars' (min/-1 = .)
gen ethnic_group = .
foreach var of local ethnic_vars{
	replace ethnic_group = `var' if !missing(`var') & missing(ethnic_group)
}
label copy ADC06E00 ethnic_group
label define ethnic_group -9 "" -8 "" -1 "" 6 "Other", modify
label values ethnic_group ethnic_group


// Parental Education
gen mother_edu_years = mmled if inrange(mmled, 1, 10)
gen father_edu_years = mfled if inrange(mfled, 1, 10)

gen mother_edu_level = mcs_mledg - 1 if inrange(mcs_mledg, 1, 2)
gen father_edu_level = mcs_fledg - 1 if inrange(mcs_fledg, 1, 2)
label define parent_edu_level 0 "Low" 1 "High"
label values *_edu_level parent_edu_level


// Father's Social Class
recode csep 	(5 = 1 "V Unskilled") (4 = 2 "IV Partly Skilled") ///
				(3 = 3 "III Skilled Manual") (2 = 4 "III Skilled Non-Manual") ///
				(1 = 5 "II Intermediate") (0 = 6 "I Professional"), ///
				gen(father_class)


// Unit Tests
sum age_*
tab ethnic_group ADC06E00, m
sum maths_* verbal_* vocab_*
corr maths_* verbal_* vocab_*
corr *resid*
tab mother_edu_years mmled, m
tab father_edu_years mfled, m
tab mother_edu_level mcs_mledg, m
tab father_edu_level mcs_fledg, m
tab mother_edu_years mother_edu_level, m
tab father_edu_years father_edu_level, m
tab father_class csep, m


// Format
keep if inrange(COUNTRY, 1, 3) // DROP NORTHERN IRELAND
keep if NOCMHH == 1 & CNUM00 == 1 // KEEP SINGLETON BIRTHS
keep if ethnic_group == 1 // KEEP WHITE INDIVIDUALS
rename mcsid id
gen survey_weight = WEIGHTGB if WEIGHTGB >= 0 //GOVWT2 if GOVWT2 >= 0
gen cohort = "2001c"
keep id cohort survey_weight male ///
	bmi_* self_report_* ///
	maths_* verbal_* vocab_* ///
	mother_* father_*
drop *_00 *_01
compress
save "${clean}/2001c_cleaned.dta", replace
capture rm "${clean}/2001c_temp.dta"
