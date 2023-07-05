

clear all

import excel "DANE Industrial Goods Expenditure_COL_AR1.xls", sheet("data") firstrow

keep if (year > 1990) & (year < 2008)

gen ln_absorb = log(absorb)
reg ln_absorb year
predict ln_abs_detr, resid

tsset year

twoway (tsline ln_abs_detr), ytitle(log real manuf. Spending) title(Log Col market size)
regress ln_abs_detr l.ln_abs_detr


