

clear all
cd  "../US AR1 data files"

log using US_AR1.txt, text replace

use US_AR1_data
keep if (year > 1990) & (year < 2008)

gen US_mkt   = peso_dollar_xr*(US_manuf_sales + imports - exports)/COL_PPI
gen lnUS_mkt = log(US_mkt)
reg lnUS_mkt year
predict lnUS_mkt_detr, resid

tsset year

twoway (tsline lnUS_mkt_detr), ytitle(log real manuf. Spending) title(Log US market size)

regress lnUS_mkt_detr l.lnUS_mkt_detr
