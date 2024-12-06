
/*******************************************************************************
File: generate extended version of Table 2
Purpose: Estimate the Phillips Curve's slope for 3 periods
*******************************************************************************/
clear all
set more off
set matsize 800

local install_packages = "no"

if "`install_packages'" == "yes" {
ssc install reghdfe
ssc install ranktest
ssc install moremata
ssc install carryforward
ssc install ivreg2
ssc install ftools
ssc install ivreghdfe
ssc install binscatter
}


*****************************************************************************
*Step 0: Set the local name with your name
*****************************************************************************

local name = "marina"

* set your path under your name
if "`name'" == "marina" { 
cd "/Users/maperrupato/Documents/UCSD/Computation/QJE Replication"
}

*The dataset contains these variables year state quarter mean_une statecode date constant infl_reg rp qt_bartik_sa
use "data/data_reg.dta", clear

*****************************************************************************
*Step 1: Calculate present values 
* The point here is to compute the PV approach for unemployment and prices. We do this for a benchmark value of beta = 0.99, and robustness values beta = 0.9, 0.95. The variables rp_sum_* and u_sum_* contain the outcome.
*****************************************************************************

* Set value of beta
global beta 0.99


* Vary truncation length
foreach truncation_length in 10 20 30 40 {

  * Calculate present value of unemployment
  quietly capture drop u_sum_`truncation_length'
  quietly generate u_sum_`truncation_length' = mean_une

  forvalues i = 1/`truncation_length' {

  	quietly replace u_sum_`truncation_length' = u_sum_`truncation_length' + $beta^`i'*F`i'.mean_une

  }


  * Calculate present value of relative prices. Similar procedure than for unemployment.
  quietly capture drop rp_sum_`truncation_length'
  quietly generate rp_sum_`truncation_length' = rp

  forvalues i = 1/`truncation_length' {

  	quietly replace rp_sum_`truncation_length' = rp_sum_`truncation_length' + $beta^`i'*F`i'.rp

  }

}

* Calculate present values for beta = 0.90, 0.95
foreach beta_local in 90 95 {

	local beta_local_adj = `beta_local' / 100

	quietly capture drop u_sum_20_`beta_local'
	quietly generate u_sum_20_`beta_local' = mean_une

	forvalues i = 1/20 {

		quietly replace u_sum_20_`beta_local' = u_sum_20_`beta_local' + `beta_local_adj'^`i'*F`i'.mean_une

	}
	quietly capture drop rp_sum_20_`beta_local'
	quietly generate rp_sum_20_`beta_local' = rp

	forvalues i = 1/20 {

		quietly replace rp_sum_20_`beta_local' = rp_sum_20_`beta_local' + `beta_local_adj'^`i'*F`i'.rp

	}
}


xtset statecode date

*****************************************************************************
* Step 2: Run Regressions 
***********************************************************************
* Generate some useful variables. infl_reg_sign just changes the sign of inflation so that the coefficients have the sign as in the paper. infl_reg_time_agg divides by 4 because the model is written in quarterly terms and the data are 12-month inflation rates. d20_qt_bartik_sa takes 20 quarter differences of the seasonally adjusted bartik instrument. L4_* variables are 4 lags of the respective variables.

generate infl_reg_sign = -1 * infl_reg
generate infl_reg_time_agg = -1 * infl_reg / 4
generate d20_qt_bartik_sa = qt_bartik_sa - L20.qt_bartik_sa
generate L4_d20_qt_bartik_sa = L4.d20_qt_bartik_sa
generate L4_rp = L4.rp
generate L4_mean_une = L4.mean_une

** Extension of Table 2: OLS and IV estimates by subperiod, with and without time fixed effects. Pretty self-explanatory.

* First we compute the estimates for Psi. Which means a regression of Pi on unemployment and relative prices, with or without the tradeable demand instrument and fixed effects. Then we compute estimates for kappa. The difference is that the RHS variable is the PV of u and rp, which are instrumented by the lagged values and the tradeable demand instrument.

*Regressions
eststo clear
xtset statecode date
* Subsample Estimate of Psi
quietly {
	* OLS IV
	* Volcker Era:
	eststo: reghdfe infl_reg_sign L4.mean_une L4.rp if (year<=1990), absorb(i.statecode i.date) cluster(statecode)
    * Greenspan Er
	eststo: reghdfe infl_reg_sign L4.mean_une L4.rp if (1990< year & year<=2006), absorb(i.statecode i.date) cluster(statecode)
	* Bernanke Era
	eststo: reghdfe infl_reg_sign L4.mean_une L4.rp if (2007<=year), absorb(i.statecode i.date) cluster(statecode)
	
    * IV Estimation
    * Volcker Era
    eststo: ivreghdfe infl_reg_sign (L4.mean_une = L4.d20_qt_bartik_sa) L4.rp if (year<=1990), absorb(i.statecode i.date) cluster(statecode)
    * Greenspan Er
    eststo: ivreghdfe infl_reg_sign (L4.mean_une = L4.d20_qt_bartik_sa) L4.rp if (1990< year & year<=2006), absorb(i.statecode i.date) cluster(statecode)
    * Bernanke Era
    eststo: ivreghdfe infl_reg_sign (L4.mean_une = L4.d20_qt_bartik_sa) L4.rp if (2007<=year), absorb(i.statecode i.date) cluster(statecode)
}

* Generate table
esttab, se keep(L4.mean_une)

estout using "output/psi_time_varying_ext.tex", style(tex) keep(L4.mean_une) varlabels(L4.mean_une "$\psi$") cells(b(star fmt(%9.3f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(L4.mean_une)

eststo clear

* Subsample Estimate of Kappa
eststo clear

quietly {
    * OLS IV
	* Volcker Era
    eststo: ts2sls infl_reg_time_agg  i.date (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), ifsecond(year <= 1990) absorb(statecode) cluster(statecode)
    * Greenspan Era
    eststo: ts2sls infl_reg_time_agg  i.date (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), ifsecond(1990 <year & year <= 2006) absorb(statecode) cluster(statecode)
	* Bernanke Era
    eststo: ts2sls infl_reg_time_agg  i.date (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), ifsecond(2007 <= year) absorb(statecode) cluster(statecode)

	
    * IV Estimation: Replace L4_mean_une with L4_d20_qt_bartik_sa as instrument
    * Volcker Era
    eststo volker: ts2sls infl_reg_time_agg i.date (u_sum_20 rp_sum_20 = L4_d20_qt_bartik_sa L4_rp), ifsecond(year <= 1990) absorb(statecode) cluster(statecode)
	scalar volcker_coef = e(b)[1, 1]
	scalar volcker_se = sqrt(e(V)[1, 1])

    * Greenspan Era
    eststo greenspan: ts2sls infl_reg_time_agg i.date (u_sum_20 rp_sum_20 = L4_d20_qt_bartik_sa L4_rp), ifsecond(1990 < year & year <= 2006) absorb(statecode) cluster(statecode)
	scalar greenspan_coef = e(b)[1, 1]
	scalar greenspan_se = sqrt(e(V)[1, 1])

	* Bernanke Era
    eststo bernanke: ts2sls infl_reg_time_agg i.date (u_sum_20 rp_sum_20 = L4_d20_qt_bartik_sa L4_rp), ifsecond(2007 <= year) absorb(statecode) cluster(statecode)
	scalar bernanke_coef = e(b)[1, 1]
	scalar bernanke_se = sqrt(e(V)[1, 1])

}

* Generate table
esttab, se keep(u_sum_20 rp_sum_20)

estout using "output/kappa_time_varying_ext.tex", style(tex) keep(u_sum_20) varlabels(u_sum_20 "$\kappa$") cells(b(star fmt(%9.4f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(u_sum_20)

***********************************************************************
*Test if difference between estimations are significant
***********************************************************************

* Compute the difference and SE of the difference
scalar diff_v = volcker_coef - bernanke_coef
scalar se_diff_v = sqrt(volcker_se^2 + bernanke_se^2)

scalar diff = greenspan_coef - bernanke_coef
scalar se_diff = sqrt(greenspan_se^2 + bernanke_se^2)

* Calculate the t-statistic
scalar t_stat_v = diff_v / se_diff_v
scalar t_stat = diff / se_diff

* Redirect output to a log file
capture log close
log using "output/t_test_results.txt", text replace

* Display and log results
display "Volker vs Bernanke Difference in estimates: " diff_v
display "Volker vs Bernanke Standard error of the difference: " se_diff_v
display "Volker vs Bernanke t-statistic: " t_stat_v

display "Difference in estimates: " diff
display "Standard error of the difference: " se_diff
display "t-statistic: " t_stat

if abs(t_stat_v) > 1.96 {
    display "Volker vs Bernanke: The estimates are statistically different at the 5% level."
} 
else {
    display "Volker vs Bernanke: The estimates are not statistically different at the 5% level."
}

if abs(t_stat) > 1.96 {
    display "Greenspan vs Bernanke: The estimates are statistically different at the 5% level."
} 
else {
    display "Greenspan vs Bernanke: The estimates are not statistically different at the 5% level."
}

* Close the log file
log close


***********************************************************************
* ORIGINAL Table 2: OLS and IV estimates by subperiod
***********************************************************************
*Regressions
eststo clear
xtset statecode date

* Subsample OLS Estimate of Psi
eststo clear
quietly {
	eststo: reghdfe infl_reg_sign L4.mean_une L4.rp if year <= 1990, absorb(i.statecode) cluster(statecode)
	eststo: reghdfe infl_reg_sign L4.mean_une L4.rp if year > 1990, absorb(i.statecode) cluster(statecode)
	eststo: reghdfe infl_reg_sign L4.mean_une L4.rp if year <= 1990, absorb(i.statecode i.date) cluster(statecode)
	eststo: reghdfe infl_reg_sign L4.mean_une L4.rp if year > 1990, absorb(i.statecode i.date) cluster(statecode)
	eststo: ivreghdfe infl_reg_sign (L4.mean_une = L4.d20_qt_bartik_sa) L4.rp if year <= 1990, absorb(i.statecode i.date) cluster(statecode)
	eststo: ivreghdfe infl_reg_sign (L4.mean_une = L4.d20_qt_bartik_sa) L4.rp if year > 1990, absorb(i.statecode i.date) cluster(statecode)

}

esttab, se keep(L4.mean_une)

estout using "output/psi_time_varying.tex", style(tex) keep(L4.mean_une) varlabels(L4.mean_une "$\psi$") cells(b(star fmt(%9.3f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(L4.mean_une)

eststo clear


* Subsample OLS Estimate of Kappa
eststo clear
quietly {

	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), ifsecond(year <= 1990) absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), ifsecond(year > 1990) absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg i.date (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), ifsecond(year <= 1990) absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg i.date (u_sum_20 rp_sum_20 = L4_mean_une L4_rp), ifsecond(year > 1990) absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg i.date (u_sum_20 rp_sum_20 = L4_d20_qt_bartik_sa L4_rp), ifsecond(year <= 1990) absorb(statecode) cluster(statecode)
	eststo: ts2sls infl_reg_time_agg i.date (u_sum_20 rp_sum_20 = L4_d20_qt_bartik_sa L4_rp), ifsecond(year > 1990) absorb(statecode) cluster(statecode)

}
esttab, se keep(u_sum_20 rp_sum_20)

estout using "output/kappa_time_varying.tex", style(tex) keep(u_sum_20) varlabels(u_sum_20 "$\kappa$") cells(b(star fmt(%9.4f)) se(par)) stats( , fmt(%7.0f %7.1f %7.2f)) nolabel replace mlabels(none) collabels(none) stardrop(u_sum_20)
