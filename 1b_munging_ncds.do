/*
1958 Cohort (NCDS)
*/
closer_wide ncds ncdsid

merge 1:1 ncdsid using ///
	"${raw}/NCDS/0y-16y/ncds0123.dta", ///
	nogen keepusing(n622 n2001 n2003 n553 n574 ///
		n1196 n1199 n1202 n1205 n236 n1687 n1811)

merge 1:1 ncdsid using ///
	"${raw}/NCDS/Derived/ncds_ParentEducationDerived.dta", ///
	nogen keepusing(n16med n16fed n016nmed n716daded)

rename ncdsid NCDSID
merge 1:1 NCDSID using ///
	"${raw}/NCDS/55y/ncds_2013_flatfile.dta", ///
	nogen keepusing(N9WT* N9HT*)
	
merge 1:1 NCDSID using ///
	"${raw}/NCDS/55y/ncds_2013_derived.dta", ///
	nogen keepusing(ND9BMI)
	
rename NCDSID id
merge 1:1 id using "${cognition}/1958c_cognition.dta", ///
	nogen keepusing(*_11 *_16)
drop general_11
	
// Cognition
// rename maths_16 maths_16
rename comprehension_16 vocab_16
// rename verbal_11 verbal_11
rename arithmetic_11 maths_11


// Sex
gen male = 2 - n622  if inrange(n622, 1, 2)

gen puberty_boy = n2001 if male == 1 & inrange(n2001, 1, 4)
replace puberty_boy = 0 if male == 0

gen puberty_girl = n2003 if male == 0 & inrange(n2003, 9, 16)
replace puberty_girl = 0 if male == 1

gen maternal_age = n553 if inrange(n553, 14, 50)

gen birth_weight = n574 * 28.3495 / 1000 if n574 >= 0


// Parental BMI
gen father_height = n1199 * 2.54 / 100 if n1199 >= 0
replace father_height = . if father_height < 1.4

gen mother_height = n1205 * 2.54 / 100 if n1205 >= 0
replace mother_height = . if mother_height < 1.4

gen father_weight = n1196 * 7 * 0.453592 if n1196 >= 0
gen mother_weight = n1202 * 7 * 0.453592 if n1202 >= 0

gen father_bmi = father_weight/(father_height ^ 2)
replace father_bmi = . if !inrange(father_bmi, 13, 70)
gen mother_bmi = mother_weight/(mother_height ^ 2)
replace mother_bmi = . if !inrange(mother_bmi, 13, 70)

drop father_height father_weight mother_height mother_weight

	
// BMI
gen bmi_55 = ND9BMI if inrange(ND9BMI, 0, 100)

gen height_55 = ( (N9HTFEET * 12) + N9HTINES )*2.54 ///
	if N9HTFEET >= 0 & N9HTINES >= 0
replace height_55 = (N9HTMEES * 100) + N9HTCMS if N9HTMEES >= 0 & N9HTCMS >= 0

gen weight_55 = ( (N9WTSTE * 14) + N9WTPOD )*0.454 ///
	if N9WTSTE >= 0 & N9WTPOD >= 0
replace weight_55 = N9WTKIS if N9WTKIS >= 0

gen self_height_55 =1 if !missing(bmi_55)
replace self_height_55 = 2 if !missing(height_55)

gen self_weight_55 = 1 if !missing(bmi_55)
replace self_weight_55 = 2 if !missing(weight_55)

gen self_report_55 = max(self_weight_55, self_height_55)
label values self_weight_55 wtself
drop self_?eight_*

foreach var of varlist self_*{
	levelsof `var'
	if r(r) == 1	drop `var'
}


// Parental Education
gen mother_edu_years = n16med if inrange(n16med, 1, 10)
gen father_edu_years = n16fed if inrange(n16fed, 1, 10)

gen mother_edu_level = n016nmed if inrange(n016nmed, 0, 1)
gen father_edu_level = n716daded if inrange(n716daded, 0, 1)
label define parent_edu_level 0 "Low" 1 "High"
label values *_edu_level parent_edu_level


// Father's Social Class
gen father_class_alt = 8 - n236 if inrange(n236, 2, 7)
gen father_class = 7 - n1687 if inrange(n1687, 1, 6)

label_class


// Unit Tests


			
// Format Dataset
keep if n1811 == 0 // SINGLETON BIRTHS

// rename NCDSID id
gen cohort = "1958c"
gen survey_weight = 1
keep id cohort survey_weight male ///
	maternal_age birth_weight puberty_* ///
	bmi_* self_* father_* mother_* *_11 *_16
drop *_00
compress
save "${clean}/1958c_cleaned.dta", replace