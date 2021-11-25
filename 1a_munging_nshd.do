/*
1946 Cohort (NSHD)
*/
use "${raw}/NSHD/46c_ht_cog.dta", clear
merge 1:1 nshdid_db1120  using "${raw}/NSHD/magels_fagels.dta", nogen
	
rename nshdid_db1120 id
tostring id, replace
merge 1:1 id using "${cognition}/1946c_cognition.dta", ///
	nogen keepusing(*_11 age_15 arithmetic_15 reading_15)
drop general_11 iq_*
	
// Cognition
rename arithmetic_15 maths_15
rename reading_15 vocab_15
// rename verbal_11 verbal_11
rename arithmetic_11 maths_11

// Sex
gen male = 2 - sex

// Birth Weight
gen birth_weight = mbwtu / 1000 if mbwtu != 9999 & !missing(mbwtu)

gen father_height = fht52 * 2.54 / 100 if inrange(fht52, 0, 100)
replace father_height = . if father_height < 1.4

gen mother_height = mht52 * 2.54 / 100 if inrange(mht52, 0, 100)
replace mother_height = . if mother_height < 1.4

// Parental Weight, Puberty, Maternal Age

// BMI
rename bmi??* bmi??
foreach var of varlist bmi?? {
	local year = substr("`var'", -2, 2)
	local age = cond(`year' < 46, 54 + `year', `year' - 46)
	local age = cond(`age' < 10, "0`age'", "`age'")
	gen bmi_`age' = `var' if inrange(`var', 0, 100)
}
recode bmi_?? (min/13 = .) (70/max = .)


// Parental Education
foreach parent in mother father{
	local p = substr("`parent'", 1, 1)
	
	gen `parent'_edu_years = .
	replace `parent'_edu_years = `p'agels - 9 if inrange(`p'agels, 10, 19)
	replace `parent'_edu_years = 10 if inrange(`p'agels, 20, 30)
	
	gen `parent'_edu_level = inrange(`p'agels, 15, 30) ///
		if inrange(`p'agels, 10, 30)
}
label define parent_edu_level 0 "Low" 1 "High"
label values *_edu_level parent_edu_level


// Father's Social Class
foreach var of varlist fsc50 fsc57{
	local suffix = cond("`var'" == "fsc50", "_alt", "")
	
	recode `var' 	(50 = 1 "V Unskilled") (40 = 2 "IV Partly Skilled") ///
					(35 = 3 "III Skilled Manual") (30 = 4 "III Skilled Non-Manual") ///
					(20 = 5 "II Intermediate") (10 = 6 "I Professional") ///
					(2 = .) (60 = .), gen(father_class`suffix')
}


// Unit Tests



// Format Dataset
// rename nshdid_db1120 id
// tostring id, replace
egen total_weight = total(inf)
gen survey_weight = inf * _N / total_weight
gen cohort = "1946c"
keep	id cohort survey_weight male ///
		birth_weight /// maternal_age puberty_* 
		bmi_* father_* mother_* *_11 *_15 // self_* 
compress
save "${clean}/1946c_cleaned.dta", replace