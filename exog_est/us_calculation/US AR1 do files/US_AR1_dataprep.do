//                          US_AR1_dataprep.do

/* This script:
  (1) converts the FRED monthly US manufacturing PPI to an annual series 
  (2) Aggregates daily peso-dollar exchange rate to annual averages. 
  (3) Constructs manufactured imports and exports
  (4) Pulls in implicit peso price deflator for Colombian manufactured goods
  (5) Merges all data series
*/

clear all
cd  "../US AR1 data files"

// US manufacturing producer price index
// Data source: FRED

use "FRED_manuf_PPI", replace
destring observation_date, ignore('-') replace
gen year = floor(observation_date/10000)

collapse pcuomfgomfg, by(year)
rename pcuomfgomfg US_manuf_PPI 

save "FRED_annual_PPI", replace

// nominal US manufacturing sales
// Data source: FRED

use "US_manuf_sales_FRED", replace
destring observation_date, ignore('-') replace
gen year = floor(observation_date/10000)

rename slmnto02usa189n US_manuf_sales
save "clean_US_manuf_sales", replace

// Colombian peso US dollar xchange rate
// Data source: Banco de la República

use "COL_US_exch_rate", replace

rename marketexchangeratecolombianpeso peso_dollar_xr 
collapse peso_dollar_xr, by(year)

save "annual_peso_dollar_xr", replace

// US trade data
// Data source: Comtrade

use "Comtrade_US_XM_data", replace
keep year tradeflowcode commoditycode tradevalueus 
keep if commoditycode > 26 // drop non-manufactured goods
preserve
  keep if tradeflowcode == 1
  gen imports = tradevalueus
  collapse (sum) imports, by(year)
  save Comtrade_US_manuf_imports, replace
 restore
 
  keep if tradeflowcode == 2
  gen exports = tradevalueus
  collapse (sum) exports, by(year)
  save Comtrade_US_manuf_exports, replace
  
  merge 1:1 year using Comtrade_US_manuf_imports

// merging all annualized time series

merge 1:1 year using annual_peso_dollar_xr, gen(_merge0) 
merge 1:1 year using clean_US_manuf_sales, gen(_merge1)
merge 1:1 year using FRED_annual_PPI, gen(_merge2)
merge 1:1 year using US_manuf_XM.dta, gen(_merge3)
merge 1:1 year using COL_PPI_Marcela.dta, gen(_merge4)
drop observation_date

save US_AR1_data, replace
