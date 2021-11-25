/*
1970 Cohort (BCS70)
*/
closer_wide bcs70 bcsid // Age 38?
save "${clean}/1970c_temp.dta", replace  
	
use BCSID SEX COB MULTIPNO using ///
	"${raw}/BCS70/xwave/bcs70_response_1970-2016.dta", clear
	
merge 1:1 BCSID using ///
	"${raw}/BCS70/46y/bcs_age46_main.dta", ///
	nogen keepusing(BD10BMI BD10MBMI)

rename BCSID bcsid
merge 1:1 bcsid using "${raw}/BCS70/0y/bcs7072a.dta", ///
	nogen keepusing(a0005a a0278 a0014 a0248 a0197)

merge 1:1 bcsid using "${raw}/BCS70/10y/sn3723.dta", ///
	nogen keepusing(e1_1 e1_2 e2_1 e2_2 c3_4)

merge 1:1 bcsid using "${raw}/BCS70/16y/bcs70_16-year_arithmetic_data.dta", ///
	nogen keepusing(mathscore)

merge 1:1 bcsid using "${raw}/BCS70/16y/bcs7016x.dta", ///
	nogen keepusing(ha1_1 ha1_2 ha7 ha7_1)
	
merge 1:1 bcsid using ///
	"${raw}/BCS70/derived/bcs70_ParentEducationDerived.dta", ///
	nogen keepusing(b016mmed b016ffed b05med b5fed2)
	
merge 1:1 bcsid using "${clean}/1970c_temp.dta", nogen

rename bcsid id
merge 1:1 id using "${cognition}/1970c_cognition.dta", ///
	nogen keepusing(*_10 age_16 home_test_16 apu_arithmetic_16 apu_vocab_a_16)
	
// Cognition
rename apu_arithmetic_16 maths_16
rename apu_vocab_a_16 vocab_16
rename bas_sim_10 verbal_10
// rename maths_10 maths_10


// Sex
gen male = 2 - SEX if inrange(SEX, 1, 2)

gen maternal_age = a0005a if inrange(a0005a, 14, 52)

gen birth_weight = a0278 / 1000 if a0278 >= 0

gen puberty_girl =  ha7_1 + 10 if male == 0 & inrange(ha7_1, 1, 6)
replace puberty_girl = 0 if male == 1

// Parental BMI
gen father_height = e2_1 / 100 if e2_1 >= 0
replace father_height = . if father_height < 1.4
gen mother_height = e1_1 / 100 if e1_1 >= 0
replace mother_height = . if mother_height < 1.4
gen father_weight = e2_2 if e2_2 >= 0
gen mother_weight = e1_2 if e1_2 >= 0

gen father_bmi = father_weight/(father_height ^ 2)
replace father_bmi = . if !inrange(father_bmi, 13, 70)
gen mother_bmi = mother_weight/(mother_height ^ 2)
replace mother_bmi = . if !inrange(mother_bmi, 13, 70)

drop father_height father_weight mother_height mother_weight


// BMI
gen bmi_46 = BD10MBMI if BD10MBMI > 0
replace bmi_46 = BD10BMI if BD10BMI > 0 & missing(bmi_46)

gen self_report_46 = 1 if BD10MBMI > 0 & !missing(BD10MBMI)
replace self_report_46 = 2 if self_report_46 != 1 & !missing(bmi_46)
label values self_report_46 wtself

foreach var of varlist self_*{
	levelsof `var'
	if r(r) == 1	drop `var'
}


// Parental Education
gen mother_edu_years = b016mmed if inrange(b016mmed, 1, 10)
gen father_edu_years = b016ffed if inrange(b016ffed, 1, 10)

gen mother_edu_level = b05med if inrange(b05med, 0, 1) // WHY NATURAL?
gen father_edu_level = b5fed2 if inrange(b5fed2, 0, 1)
label define parent_edu_level 0 "Low" 1 "High"
label values *_edu_level parent_edu_level


// Father's Social Class
gen father_class_alt = 7 - a0014 if inrange(a0014, 1, 6)
gen father_class = 7 - c3_4 if inrange(c3_4, 1, 6)

label_class
 

// Unit Tests 
 

// Format Dataset
keep if inrange(COB, 1, 3) // Born in England, Scotland or Wales
keep if MULTIPNO == -1 | a0248 == 1 // Singleton births

// rename bcsid id
gen cohort = "1970c"
gen survey_weight = 1
keep id cohort survey_weight male ///
	maternal_age birth_weight puberty_* ///
	bmi_* self_* father_* mother_* *_10 *_16
drop *_00 *_05
compress 
save "${clean}/1970c_cleaned.dta", replace 
capture rm "${clean}/1970c_temp.dta"