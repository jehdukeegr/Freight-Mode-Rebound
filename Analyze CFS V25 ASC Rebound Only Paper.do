* This file analyzes the 2012 public-use microdata file of the Commodity Flow Survey
* Creates the figures and tables for the paper
* Implements simulation of mode choices with random extreme value errors
* Implements the Truck CAFE scenarios
** Changes from previous versions: Incorporates 2017 CFS data
* Version 21 iterative eliminates simulations

drop _all
clear matrix
clear mata
set mem 16g
set seed 1977

* Pathways
* Mac
*global PATH = "/Users/Jon/Dropbox/2012 Commodity Flow Survey"
*cd "$PATH"
* PC
global PATH = "d:\Dropbox\2012 Commodity Flow Survey"
cd "$PATH"
  
* 5 improvement in truck fuel economy coming from CAFE
global CAFE = 0.05
* Fuel intensities in gallons per ton-mile
global gptmTruck = (1/85) 
global gptmRail = (1/450)
global gptmBarge = (1/600)
global gptmAir = (1/7.5)
* Fuel emissions intensity in metric tons per gallon (https://www.eia.gov/environment/emissions/co2_vol_mass.php)
global fuelMTpg = 10.16/1000
global fuelMTpgJet = 9.57/1000
* Pass Through rate - assume average full pass through per Marion and Muehlegger (2011)
global passThruRate = 1
* Switch for simplified mode-choice probability pictures
global prettyPics = 1
* Switch for distance "penalties" for mode switching
global distPenalty = 1
* Use 2017 data
global useboth = 1
* Switch to use more expansive inland river system definition
global moreRivers = 1

**** Summary statistic type tables ****
* Characteristics by mode
*use "Data/CFS2012.dta", clear
*if $useboth == 1{
*append using "Data/CFS2017.dta",
*}
*gen weight = shipmt_wght
*gen distance = shipmt_dist_routed
*tabstat value vpt distance weight tonmiles, by(mode_txt) labelwidth(16)
*sort mode_txt
*collapse (mean) value vpt distance weight tonmiles [iweight = wgt_factor], by(mode_txt)
*outsheet using "Tables/ModeChars1217.csv", comma replace

**** Calculate distance penalties across modes ****
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
*append using "Data/CFS2017.dta",
keep if mode == 4 | mode == 6 | mode == 8  | mode == 11 
*  Drop "Reminder of State" shipments
drop if orig_ma == 99999
drop if dest_ma == 99999
drop if orig_ma == dest_ma
sort sctg orig_cfs_area dest_cfs_area mode
collapse (mean) miles (first) orig_state_txt dest_state_txt orig_description dest_description [iweight = wgt_factor], by(sctg orig_cfs_area dest_cfs_area mode)
gen Route =  orig_cfs_area + "_to_" + dest_cfs_area
* Reshape to wide
reshape wide miles orig_state_txt dest_state_txt orig_description dest_description, i(sctg Route) j(mode)
* "To truck" switching
gen TrucktoTruck = miles4/miles4 
gen RailtoTruck = miles6/miles4 
gen WatertoTruck = miles8/miles4
gen AirtoTruck = miles1/miles4
* "To rail" switching
gen TrucktoRail = miles4/miles6 
gen RailtoRail = miles6/miles6 
gen WatertoRail = miles8/miles6
gen AirtoRail = miles1/miles6
* "To water" switching
gen TrucktoBarge = miles4/miles8 
gen RailtoBarge = miles6/miles8
gen WatertoBarge = miles8/miles8
gen AirtoBarge = miles1/miles8
* "To air" switching
gen TrucktoAir = miles4/miles11
gen RailtoAir = miles6/miles11
gen WatertoAir = miles8/miles11
gen AirtoAir = miles1/miles11
*** Summary stats, use median
tabstat RailtoTruck, stats(n mean p5 p25 p50 p75 p95)
tabstat WatertoTruck, stats(n mean p5 p25 p50 p75 p95)
tabstat AirtoTruck, stats(n mean p5 p25 p50 p75 p95)
** Hard code globals
* Penalties are median values from route by sctg mean values
global TruckRailP = 1.15
global TruckWaterP = 1.05
global TruckAirP = 0.95
*** Summary stats, use median
tabstat TrucktoRail, stats(n mean p5 p25 p50 p75 p95)
tabstat WatertoRail, stats(n mean p5 p25 p50 p75 p95)
tabstat AirtoRail, stats(n mean p5 p25 p50 p75 p95)
** Hard code globals
* Penalties are median values from route by sctg mean values
global RailTruckP =0.87
global RailWaterP = 0.93
global RailAirP = 0.83
*** Summary stats, use median
tabstat TrucktoBarge, stats(n mean p5 p25 p50 p75 p95)
tabstat RailtoBarge, stats(n mean p5 p25 p50 p75 p95)
tabstat AirtoBarge, stats(n mean p5 p25 p50 p75 p95)
** Hard code globals
* Penalties are median values from route by sctg mean values
global BargeTruckP =0.95
global BargeRailP = 1.07
global BargeAirP = 0.99
*** Summary stats, use median
tabstat TrucktoAir, stats(n mean p5 p25 p50 p75 p95)
tabstat RailtoAir, stats(n mean p5 p25 p50 p75 p95)
tabstat WatertoAir, stats(n mean p5 p25 p50 p75 p95)
** Hard code globals
* Penalties are median values from route by sctg mean values
global AirTruckP = 1.05
global AirRailP = 1.20
global AirWaterP = 1.01

** Appendix Table 1
* Unweighted
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
sort sctg_txt
gen Truck = mode == 4
gen Rail = mode == 6
gen Air = mode == 11
gen Water = mode == 8
gen Pipeline = mode == 12
gen ParcelCourier = mode == 14
rename value Value
rename shipmt_dist_routed Miles
replace shipmt_wght = shipmt_wght/2000
rename shipmt_wght Tons
rename sctg SCTG
format Value %9.2f
format Miles %9.2f
format Tons %9.2f
format Air %9.2f
format Pipeline %9.2f
format Rail %9.2f
format Truck %9.2f
format Water %9.2f
format ParcelCourier %9.2f
#delimit ;
collapse (mean) Value_mean=Value Miles_mean=Miles Tons_mean=Tons 
Air_mean=Air Pipeline_mean=Pipeline Rail_mean=Rail 
Truck_mean=Truck Water_mean=Water ParcelCourier_mean=ParcelCourier 
(sd) Value_stdev=Value Miles_stdev=Miles Tons_stdev=Tons 
Air_stdev=Air Pipeline_stdev=Pipeline Rail_stdev=Rail 
Truck_stdev=Truck Water_stdev=Water ParcelCourier_stdev=ParcelCourier 
(min) Value_min=Value Miles_min=Miles Tons_min=Tons 
Air_min=Air Pipeline_min=Pipeline Rail_min=Rail 
Truck_min=Truck Water_min=Water  ParcelCourier_minn=ParcelCourier 
(max) Value_max=Value Miles_max=Miles Tons_max=Tons 
Air_max=Air Pipeline_max=Pipeline Rail_max=Rail 
Truck_max=Truck Water_max=Water ParcelCourier_max=ParcelCourier 
[iweight = wgt_factor], ;
#delimit cr
xpose, clear varname
order _varname
format v1 %9.4f
outsheet using "Tables/SumStatsUnWgt.csv", comma replace
* Ton-mile weighted
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
replace wgt_factor = wgt_factor*tonmiles
sort sctg_txt
gen Truck = mode == 4
gen Rail = mode == 6
gen Air = mode == 11
gen Water = mode == 8
gen Pipeline = mode == 12
gen ParcelCourier = mode == 14
rename value Value
rename shipmt_dist_routed Miles
replace shipmt_wght = shipmt_wght/2000
rename shipmt_wght Tons
rename sctg SCTG
format Value %9.2f
format Miles %9.2f
format Tons %9.2f
format Air %9.2f
format Pipeline %9.2f
format Rail %9.2f
format Truck %9.2f
format Water %9.2f
format ParcelCourier %9.2f
#delimit
collapse (mean) Value_mean=Value Miles_mean=Miles Tons_mean=Tons 
Air_mean=Air Pipeline_mean=Pipeline Rail_mean=Rail 
Truck_mean=Truck Water_mean=Water ParcelCourier_mean=ParcelCourier
(sd) Value_stdev=Value Miles_stdev=Miles Tons_stdev=Tons 
Air_stdev=Air Pipeline_stdev=Pipeline Rail_stdev=Rail 
Truck_stdev=Truck Water_stdev=Water ParcelCourier_stdev=ParcelCourier
(min) Value_min=Value Miles_min=Miles Tons_min=Tons 
Air_min=Air Pipeline_min=Pipeline Rail_min=Rail 
Truck_min=Truck Water_min=Water ParcelCourier_min=ParcelCourier  
(max) Value_max=Value Miles_max=Miles Tons_max=Tons 
Air_max=Air Pipeline_max=Pipeline Rail_max=Rail 
Truck_max=Truck Water_max=Water ParcelCourier_max=ParcelCourier 
[iweight = wgt_factor], ;
#delimit cr
xpose, clear varname
order _varname
format v1 %9.4f
outsheet using "Tables/SumStatsTmWgt.csv", comma replace

*** Figure 1a
* By SCTG
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
replace wgt_factor = wgt_factor*tonmiles
sort sctg_txt
gen Truck = mode == 4
gen Rail = mode == 6
gen Air = mode == 11
gen Water = mode == 8
gen Pipeline = mode == 12
gen ParcelCourier = mode == 14
rename vpt Value
rename shipmt_dist_routed Miles
replace shipmt_wght = shipmt_wght/2000
rename shipmt_wght Tons
rename sctg SCTG
format Value %9.2f
format Miles %9.2f
format Tons %9.2f
format Air %9.2f
format Pipeline %9.2f
format Rail %9.2f
format Truck %9.2f
format Water %9.2f
format ParcelCourier %9.2f
collapse (mean) SCTG Value Miles Tons Air Pipeline Rail Truck Water ParcelCourier tonmiles [iweight = wgt_factor], by(sctg_txt) 
drop if SCTG == .
save "Data/SumStatsTonMile.dta", replace
outsheet using "Tables/SumStatsTonMile.csv", comma replace

*** Figure 1b
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
sort mode_txt
gen distance = shipmt_dist_routed
gen weight = shipmt_wght
collapse (mean) value vpt distance weight tonmiles tons [iweight = wgt_factor], by(mode_txt)
drop if mode_txt == "Great Lakes" | mode_txt == "Rail and water" | mode_txt == "Water" | mode_txt == "Deep Sea" 
* Constructed on mean values
gen logtonmiles = log(tonmiles)
replace tonmiles = tonmiles/1000
* Modes by value and size
#delimit  
twoway (scatter vpt tonmiles if  
mode_txt=="Air"|
mode_txt=="Rail" |
mode_txt=="Parcel, USPS, or courier" |
mode_txt=="Great Lakes" |
mode_txt=="Water" |
mode_txt=="Inland Water" |
mode_txt=="Rail and water" |
mode_txt=="Truck and water" |
mode_txt=="Truck and rail" |
mode_txt=="Deep Sea" |
mode_txt=="Pipeline" |
mode_txt=="Truck", mlabposition(1)
xtitle("Ton-miles (1000s)") scheme(s1mono) xscale(r(0 315)) xlabel(0 100 200 300 )
yscale(log) ylabel(1 10 100 1000) mcolor(purple%50) 
ytitle("Value ($/lb.)") mlabel(mode_txt) title("Mean Shipment Value and Size by Mode"));
#delimit cr
graph export "Figures/Mode_Value_Size.pdf", replace

*** Figure 1c
* Share pictures for individual commodities - Color Versions
* Loop over commodities
#delimit ;
foreach g in "Grain" { ;
#delimit cr
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
keep if sctg_txt == "`g'"
keep if mode == 4 | mode == 6 | mode == 8 
sort tonmiles
xtile tonmile_decile = tonmiles, nq(10)
collapse (sum) wgt_factor, by(tonmile_decile mode)
rename wgt_factor shipments
by tonmile_decile: egen tot_ship= sum(shipments)
gen modeshare = shipments/tot_ship 
keep mode tonmile_decile modeshare
* Bar graph
#delimit ;
graph bar (sum) modeshare,
percentages
over(mode)
over(tonmile_decile)
asyvars stack
legend(cols(3) size(small) colfirst)
ytitle("Mode Share (%)") scheme(s1mono) 
bar(1, color(purple%50)) bar(2, color(navy%50)) bar(3, color(emerald%50))
title("`g': Mode Shares by Ton Mile Decile")
;
#delimit cr
graph export "Figures/`g'ShareDeciles_color.pdf", replace
}

*** Figure 2 and Appendix Table A2
* Implementation using ASC logit
* Grain
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
keep if sctg_txt == "Grain"
keep if mode == 4 | mode == 6 | mode == 8

* Reshape data for conditional logit implementation of mixed logit
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt8 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmBarge*Pfuel*tonmiles/10^6 if mode2 == 8
qui: save "Data/MixedLogitTemp.dta", replace
*** Estimate Mixed Logit and Save Parameter estimates
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles Miss) base(6)
eststo GrainResults
* Table for parameter estimates
#delimit;
estout GrainResults using "Tables/GrainResults_ASC.txt", 
mlabels("Grain") cells(b(fmt(%9.3f)) se(par fmt(%9.6f)) ) style(tab)  
stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace ;
#delimit cr
* Figure
* Focus on 2012
keep if year == 2012
*** Estimate Mixed Logit without Miss dummy to improve exposition
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(6)
predict ASCprobs, pr
* Improve Truck CAFE
replace PdTM = PdTM*(1-$CAFE) if mode2 == 4
predict ASCprobsCAFE, pr
* Reshape
keep tonmiles mode2 shipmt_id ASCprobs ASCprobsCAFE
reshape wide ASCprobs ASCprobsCAFE, i(shipmt_id) j(mode2)
sort tonmiles
** Graphs
* Fits
#delimit ;
twoway 
(scatter ASCprobs4 tonmiles if tonmiles < 70000, msymbol(O) mcolor(red*.25) msize(tiny)) 
(scatter ASCprobs6 tonmiles if tonmiles < 70000, msymbol(D) mcolor(green*.25) msize(tiny)) 
(scatter ASCprobs8 tonmiles if tonmiles < 70000, msymbol(T) mcolor(orange_red*.25) msize(tiny))
(scatter ASCprobsCAFE4 tonmiles if tonmiles < 70000, msymbol(O) mcolor(red*.5) msize(tiny)) 
(scatter ASCprobsCAFE6 tonmiles if tonmiles < 70000, msymbol(D) mcolor(green*.5) msize(tiny)) 
(scatter ASCprobsCAFE8 tonmiles if tonmiles < 70000, msymbol(T) mcolor(orange_red*.5) msize(tiny))
(lowess ASCprobs4 tonmiles if tonmiles < 70000, lcolor(red) lpattern(dash) bwidth(.05)) 
(lowess ASCprobs6 tonmiles if tonmiles < 70000, lcolor(green) lpattern(dash) bwidth(.05)) 
(lowess ASCprobs8 tonmiles if tonmiles < 70000, lcolor(orange_red) lpattern(dash) bwidth(.05))
(lowess ASCprobsCAFE4 tonmiles if tonmiles < 70000, lcolor(red) bwidth(.05)) 
(lowess ASCprobsCAFE6 tonmiles if tonmiles < 70000, lcolor(green) bwidth(.05)) 
(lowess ASCprobsCAFE8 tonmiles if tonmiles < 70000, lcolor(orange_red) bwidth(.05) scheme(s1mono) 
title("Grain Mode Choice Probabilities") xtitle("Ton-Miles") ytitle("Pr[mode j]")  
ylabel(, angle(horizontal)) yscale(titlegap(*+7))
legend(order(10 "Truck" 11 "Rail" 12 "Inland water" )  row(1)));
#delimit cr
graph export "Figures/GrainPredProbCAFEScatter.pdf", replace
qui: gen sctg_txt = "Grain"
qui: save "Data/SimGraphGrain.dta", replace

* Coal
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
keep if sctg_txt == "Coal"
keep if mode == 4 | mode == 6 | mode == 8

* Reshape data for conditional logit implementation of mixed logit
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt8 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmBarge*Pfuel*tonmiles/10^6 if mode2 == 8
qui: save "Data/MixedLogitTemp.dta", replace
*** Estimate Mixed Logit and Save Parameter estimates
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles Miss) base(6)
eststo CoalResults
* Table for parameter estimates
#delimit;
estout CoalResults using "Tables/CoalResults_ASC.txt", 
mlabels("Coal") cells(b(fmt(%9.3f)) se(par fmt(%9.6f)) ) style(tab)  
stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace ;
#delimit cr
* Figure
* Focus on 2012
keep if year == 2012
*** Estimate Mixed Logit without Miss dummy to improve exposition
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles Miss) base(6)
predict ASCprobs, pr
* Improve Truck CAFE
replace PdTM = PdTM*(1-$CAFE) if mode2 == 4
predict ASCprobsCAFE, pr
* Reshape
keep tonmiles mode2 shipmt_id ASCprobs ASCprobsCAFE
reshape wide ASCprobs ASCprobsCAFE, i(shipmt_id) j(mode2)
sort tonmiles
* Rescale to millions of tonmiles
replace tonmiles = tonmiles/10^6
* Fits
#delimit ;
twoway 
(scatter ASCprobs4 tonmiles if tonmiles < 6, msymbol(O) mcolor(red*.25) msize(tiny)) 
(scatter ASCprobs6 tonmiles if tonmiles < 6, msymbol(D) mcolor(green*.25) msize(tiny)) 
(scatter ASCprobs8 tonmiles if tonmiles < 6, msymbol(T) mcolor(orange_red*.25) msize(tiny))
(scatter ASCprobsCAFE4 tonmiles if tonmiles < 6, msymbol(O) mcolor(red*.5) msize(tiny)) 
(scatter ASCprobsCAFE6 tonmiles if tonmiles < 6, msymbol(D) mcolor(green*.5) msize(tiny)) 
(scatter ASCprobsCAFE8 tonmiles if tonmiles < 6, msymbol(T) mcolor(orange_red*.5) msize(tiny))
(lowess ASCprobs4 tonmiles if tonmiles < 6, lcolor(red) lpattern(dash) bwidth(.2)) 
(lowess ASCprobs6 tonmiles if tonmiles < 6, lcolor(green) lpattern(dash) bwidth(.2)) 
(lowess ASCprobs8 tonmiles if tonmiles < 6, lcolor(orange_red) lpattern(dash) bwidth(.2))
(lowess ASCprobsCAFE4 tonmiles if tonmiles < 6, lcolor(red) bwidth(.2)) 
(lowess ASCprobsCAFE6 tonmiles if tonmiles < 6, lcolor(green) bwidth(.2)) 
(lowess ASCprobsCAFE8 tonmiles if tonmiles < 6, lcolor(orange_red) bwidth(.2) scheme(s1mono) 
title("Coal Mode Choice Probabilities") xtitle("Ton-Miles") ytitle("Pr[mode j]")  
ylabel(, angle(horizontal)) yscale(titlegap(*+7))
legend(order(10 "Truck" 11 "Rail" 12 "Inland water" )  row(1)));
#delimit cr
graph export "Figures/CoalPredProbCAFEScatter.pdf", replace
qui: gen sctg_txt = "Coal"
qui: save "Data/SimGraphCoal.dta", replace

* Alcohol
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
keep if sctg_txt == "Alcohol"
keep if mode == 4 | mode == 6 
* Reshape data for conditional logit implementation of mixed logit
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: save "Data/MixedLogitTemp.dta", replace
*** Estimate Mixed Logit and Save Parameter estimates
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles TempCont) base(6)
eststo AlcoholResults
* Table for parameter estimates
#delimit;
estout AlcoholResults using "Tables/AlcoholResults_ASC.txt", 
mlabels("Alcohol") cells(b(fmt(%9.3f)) se(par fmt(%9.6f)) ) style(tab)  
stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace ;
#delimit cr
* Figure
* Focus on 2012
keep if year == 2012
*** Estimate Mixed Logit without Temp Cont to improve exposition
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles TempCont) base(6)
predict ASCprobs, pr
* Improve Truck CAFE
replace PdTM = PdTM*(1-$CAFE) if mode2 == 4
predict ASCprobsCAFE, pr
* Reshape
keep tonmiles mode2 shipmt_id ASCprobs ASCprobsCAFE
reshape wide ASCprobs ASCprobsCAFE, i(shipmt_id) j(mode2)
sort tonmiles
* Fits
#delimit ;
twoway 
(scatter ASCprobs4 tonmiles if tonmiles < 100000, msymbol(O) mcolor(red*.25) msize(tiny)) 
(scatter ASCprobs6 tonmiles if tonmiles < 100000, msymbol(D) mcolor(green*.25) msize(tiny)) 
(scatter ASCprobsCAFE4 tonmiles if tonmiles < 100000, msymbol(O) mcolor(red*.5) msize(tiny)) 
(scatter ASCprobsCAFE6 tonmiles if tonmiles < 100000, msymbol(D) mcolor(green*.5) msize(tiny)) 
(lowess ASCprobs4 tonmiles if tonmiles < 100000, lcolor(red) lpattern(dash) bwidth(.2)) 
(lowess ASCprobs6 tonmiles if tonmiles < 100000, lcolor(green) lpattern(dash) bwidth(.2)) 
(lowess ASCprobsCAFE4 tonmiles if tonmiles < 100000, lcolor(red) bwidth(.2)) 
(lowess ASCprobsCAFE6 tonmiles if tonmiles < 100000, lcolor(green) bwidth(.2) scheme(s1mono) 
title("Alcohol Mode Choice Probabilities") xtitle("Ton-Miles") ytitle("Pr[mode j]")  
ylabel(, angle(horizontal)) yscale(titlegap(*+7))
legend(order(7 "Truck" 8 "Rail" )  row(1)));
#delimit cr
graph export "Figures/AlcoholPredProbCAFEScatter.pdf", replace
qui: gen sctg_txt = "Alcohol"
qui: save "Data/SimGraphAlcohol.dta", replace

* Precision Instruments
use "Data/CFS2012.dta", clear
if $useboth == 1{
append using "Data/CFS2017.dta",
}
keep if sctg_txt == "Precision Instruments"
keep if mode == 4 | mode == 11 
* Reshape data for conditional logit implementation of mixed logit
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
qui: save "Data/MixedLogitTemp.dta", replace
*** Estimate Mixed Logit and Save Parameter estimates
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(4)
eststo InstrumentsResults
* Table for parameter estimates
#delimit;
estout InstrumentsResults using "Tables/InstrumentsResults_ASC.txt", 
mlabels("Instruments") cells(b(fmt(%9.3f)) se(par fmt(%9.6f)) ) style(tab)  
stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace ;
#delimit cr
* Figure
* Focus on 2012
keep if year == 2012
*** Estimate Mixed Logit without Temp Cont to improve exposition
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(4)
predict ASCprobs, pr
* Improve Truck CAFE
replace PdTM = PdTM*(1-$CAFE) if mode2 == 4
predict ASCprobsCAFE, pr
* Reshape
keep tonmiles mode2 shipmt_id ASCprobs ASCprobsCAFE
reshape wide ASCprobs ASCprobsCAFE, i(shipmt_id) j(mode2)
sort tonmiles
* Fits
#delimit ;
twoway 
(scatter ASCprobs4 tonmiles if tonmiles < 10000, msymbol(O) mcolor(red*.25) msize(tiny)) 
(scatter ASCprobs11 tonmiles if tonmiles < 10000, msymbol(D) mcolor(blue*.25) msize(tiny)) 
(scatter ASCprobsCAFE4 tonmiles if tonmiles < 10000, msymbol(O) mcolor(red*.5) msize(tiny)) 
(scatter ASCprobsCAFE11 tonmiles if tonmiles < 10000, msymbol(D) mcolor(blue*.5) msize(tiny)) 
(lowess ASCprobs4 tonmiles if tonmiles < 10000, lcolor(red) lpattern(dash) bwidth(.2)) 
(lowess ASCprobs11 tonmiles if tonmiles < 10000, lcolor(blue) lpattern(dash) bwidth(.2)) 
(lowess ASCprobsCAFE4 tonmiles if tonmiles < 10000, lcolor(red) bwidth(.2)) 
(lowess ASCprobsCAFE11 tonmiles if tonmiles < 10000, lcolor(blue) bwidth(.2) scheme(s1mono) 
title("Precision Mode Choice Probabilities") xtitle("Ton-Miles") ytitle("Pr[mode j]")  
ylabel(, angle(horizontal)) yscale(titlegap(*+7))
legend(order(7 "Truck" 8 "Air" )  row(1)));
#delimit cr
graph export "Figures/InstrumentsPredProbCAFEScatter.pdf", replace
qui: gen sctg_txt = "Instruments"
qui: save "Data/SimGraphInstruments.dta", replace

* Combine graphs into one for better presentation
use "Data/SimGraphGrain.dta", clear
append using "Data/SimGraphCoal.dta",
append using "Data/SimGraphAlcohol.dta",
append using "Data/SimGraphInstruments.dta",
** Graphs
* Grain
#delimit ;
twoway 
(scatter ASCprobs4 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", msymbol(O) mcolor(purple%5) msize(tiny)) 
(scatter ASCprobs6 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", msymbol(D) mcolor(navy%5) msize(tiny)) 
(scatter ASCprobs8 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", msymbol(T) mcolor(emerald%5) msize(tiny))
(scatter ASCprobsCAFE4 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", msymbol(O) mcolor(purple%20) msize(tiny)) 
(scatter ASCprobsCAFE6 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", msymbol(D) mcolor(navy%20) msize(tiny)) 
(scatter ASCprobsCAFE8 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", msymbol(T) mcolor(emerald%20) msize(tiny))
(lowess ASCprobs4 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", lcolor(purple) lpattern(dash) lwidth(medthick) bwidth(.05)) 
(lowess ASCprobs6 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", lcolor(navy) lpattern(dash) lwidth(medthick) bwidth(.05)) 
(lowess ASCprobs8 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", lcolor(emerald) lpattern(dash) lwidth(medthick) bwidth(.05))
(lowess ASCprobsCAFE4 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", lcolor(purple) lpattern(solid) bwidth(.05)) 
(lowess ASCprobsCAFE6 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", lcolor(navy) lpattern(solid) bwidth(.05)) 
(lowess ASCprobsCAFE8 tonmiles if tonmiles < 70000 & sctg_txt == "Grain", lcolor(emerald) lpattern(solid) bwidth(.05) scheme(s1mono) 
title("Grain Mode Choice Probabilities") xtitle("Ton-Miles") ytitle("Pr[mode j]")  
ylabel(, angle(horizontal)) yscale(titlegap(*+7)) xscale(r(0 70000)) xlab(0(20000)70000) 
legend(order(10 "Truck" 11 "Rail" 12 "Inland Water" ) row(1)) saving(grain, replace));
#delimit cr
graph export "Figures/GrainPredProbCAFEScatterV2.pdf", replace
* Coal
#delimit ;
twoway 
(scatter ASCprobs4 tonmiles if tonmiles < 6 & sctg_txt == "Coal", msymbol(O) mcolor(purple%5) msize(tiny)) 
(scatter ASCprobs6 tonmiles if tonmiles < 6 & sctg_txt == "Coal", msymbol(D) mcolor(navy%5) msize(tiny)) 
(scatter ASCprobs8 tonmiles if tonmiles < 6 & sctg_txt == "Coal", msymbol(T) mcolor(emerald%5) msize(tiny))
(scatter ASCprobsCAFE4 tonmiles if tonmiles < 6 & sctg_txt == "Coal", msymbol(O) mcolor(purple%20) msize(tiny)) 
(scatter ASCprobsCAFE6 tonmiles if tonmiles < 6 & sctg_txt == "Coal", msymbol(D) mcolor(navy%20) msize(tiny)) 
(scatter ASCprobsCAFE8 tonmiles if tonmiles < 6 & sctg_txt == "Coal", msymbol(T) mcolor(emerald%20) msize(tiny))
(lowess ASCprobs4 tonmiles if tonmiles < 6 & sctg_txt == "Coal", lcolor(purple) lpattern(dash) lwidth(medthick) bwidth(.2)) 
(lowess ASCprobs6 tonmiles if tonmiles < 6 & sctg_txt == "Coal", lcolor(navy) lpattern(dash) lwidth(medthick) bwidth(.2)) 
(lowess ASCprobs8 tonmiles if tonmiles < 6 & sctg_txt == "Coal", lcolor(emerald) lpattern(dash) lwidth(medthick) bwidth(.2))
(lowess ASCprobsCAFE4 tonmiles if tonmiles < 6 & sctg_txt == "Coal", lcolor(purple) lpattern(solid) bwidth(.2)) 
(lowess ASCprobsCAFE6 tonmiles if tonmiles < 6 & sctg_txt == "Coal", lcolor(navy) lpattern(solid) bwidth(.2)) 
(lowess ASCprobsCAFE8 tonmiles if tonmiles < 6 & sctg_txt == "Coal", lcolor(emerald) lpattern(solid) bwidth(.2) scheme(s1mono) 
title("Coal Mode Choice Probabilities") xtitle("Ton-Miles (millions)") ytitle("Pr[mode j]")  
ylabel(, angle(horizontal)) yscale(titlegap(*+7))
legend(order(10 "Truck" 11 "Rail" 12 "Inland Water" ) row(1)) saving(coal, replace));
#delimit cr
graph export "Figures/CoalPredProbCAFEScatterV2.pdf", replace
* Alcohol
#delimit ;
twoway 
(scatter ASCprobs4 tonmiles if tonmiles < 100000 & sctg_txt == "Alcohol", msymbol(O) mcolor(purple%5) msize(tiny)) 
(scatter ASCprobs6 tonmiles if tonmiles < 100000 & sctg_txt == "Alcohol", msymbol(D) mcolor(navy%5) msize(tiny)) 
(scatter ASCprobsCAFE4 tonmiles if tonmiles < 100000 & sctg_txt == "Alcohol", msymbol(O) mcolor(purple%20) msize(tiny)) 
(scatter ASCprobsCAFE6 tonmiles if tonmiles < 100000 & sctg_txt == "Alcohol", msymbol(D) mcolor(navy%20) msize(tiny)) 
(lowess ASCprobs4 tonmiles if tonmiles < 100000 & sctg_txt == "Alcohol", lcolor(purple) lpattern(dash) lwidth(medthick) bwidth(.2)) 
(lowess ASCprobs6 tonmiles if tonmiles < 100000 & sctg_txt == "Alcohol", lcolor(navy) lpattern(dash) lwidth(medthick) bwidth(.2)) 
(lowess ASCprobsCAFE4 tonmiles if tonmiles < 100000 & sctg_txt == "Alcohol", lcolor(purple) lpattern(solid) bwidth(.2)) 
(lowess ASCprobsCAFE6 tonmiles if tonmiles < 100000 & sctg_txt == "Alcohol", lcolor(navy) lpattern(solid) bwidth(.2) scheme(s1mono) 
title("Alcohol Mode Choice Probabilities") xtitle("Ton-Miles") ytitle("Pr[mode j]")  
ylabel(, angle(horizontal)) yscale(titlegap(*+7))
legend(order(7 "Truck" 8 "Rail"  ) row(1)) saving(alcohol, replace));
#delimit cr
graph export "Figures/AlcoholPredProbCAFEScatterV2.pdf", replace
* Instruments
#delimit ;
twoway 
(scatter ASCprobs11 tonmiles if tonmiles < 10000 & sctg_txt == "Instruments", msymbol(D) mcolor(gold%5) msize(tiny)) 
(scatter ASCprobs4 tonmiles if tonmiles < 10000 & sctg_txt == "Instruments", msymbol(O) mcolor(purple%5) msize(tiny)) 
(scatter ASCprobsCAFE11 tonmiles if tonmiles < 10000 & sctg_txt == "Instruments", msymbol(D) mcolor(gold%20) msize(tiny)) 
(scatter ASCprobsCAFE4 tonmiles if tonmiles < 10000 & sctg_txt == "Instruments", msymbol(O) mcolor(purple%20) msize(tiny)) 
(lowess ASCprobs11 tonmiles if tonmiles < 10000 & sctg_txt == "Instruments", lcolor(gold) lpattern(dash) lwidth(medthick) bwidth(.2)) 
(lowess ASCprobs4 tonmiles if tonmiles < 10000 & sctg_txt == "Instruments", lcolor(purple) lpattern(dash) lwidth(medthick) bwidth(.2)) 
(lowess ASCprobsCAFE4 tonmiles if tonmiles < 10000 & sctg_txt == "Instruments", lcolor(purple) lpattern(solid) bwidth(.2)) 
(lowess ASCprobsCAFE11 tonmiles if tonmiles < 10000 & sctg_txt == "Instruments", lcolor(gold) lpattern(solid) bwidth(.2) scheme(s1mono) 
title("Precision Instruments Mode Choice Probabilities") xtitle("Ton-Miles") ytitle("Pr[mode j]")  
ylabel(, angle(horizontal)) yscale(titlegap(*+7))
legend(order(7 "Truck" 8 "Air" ) row(1)) saving(instruments, replace));
#delimit cr
graph export "Figures/InstrumentsPredProbCAFEScatterV2.pdf", replace

*****************************************************************************************************
*** Main simulation results
**** Mixed logit models for CAFE simulations
*** 1. Models with truck, rail and barge as main modes
** Temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Agricultural Products" "Basic Chemicals" 
"Fertilizers" "Other Coal and Petroleum" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 | mode == 8

if $moreRivers == 1 {
replace Miss = 1 if orig_state_txt == "AL"
replace Miss = 1 if orig_state_txt == "PA"
replace Miss = 1 if orig_state_txt == "WV"
replace Miss = 1 if orig_state_txt == "OH"
replace Miss = 1 if orig_state_txt == "IN"
}

** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt8 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmBarge*Pfuel*tonmiles/10^6 if mode2 == 8
* Create variables for predictions
qui: gen PBAU = 0
qui: gen PCAFE = 0

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles Miss TempCont) base(6) cluster(routeid)
qui: matrix B = e(b)
qui: eststo Results
* Table for parameter estimates
estout Results using "Tables/`g'Results_ASC.txt", mlabels("`g'") cells(b(fmt(%9.3f)) t(par fmt(%9.3f)) ) style(tab) stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace

** Calculate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAU = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5])/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9])) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAU = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9])) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PBAU = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9])/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9])) if mode2 == 8;
#delimit cr

** Truck CAFE Scenario
* Truck
#delimit ;
qui: replace PCAFE = exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9])) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PCAFE = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9])) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PCAFE = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9] )) if mode2 == 8;
#delimit cr

*** Truck CAFE Calculations
gen tonmilesADJ = tonmiles
** Find most probably mode
* BAU
sort id
qui: by id: egen tempBAU = max(PBAU)
qui: gen maxprobBAU = tempBAU == PBAU
* CAFE
qui: by id: egen tempCAFE = max(PCAFE)
qui: gen maxprobCAFE = tempCAFE == PCAFE
** Correction for distance penalties
if $distPenalty == 1 {
sort shipmt_id mode2
qui: gen temp1 = 0
qui: replace temp1 = mode2 if maxprobBAU== 1
qui: by shipmt_id: egen BAUmode = max(temp1)
qui: drop temp1
** Adjust ton-miles should shipment switch mode
* Truck
qui: replace tonmilesADJ = tonmilesADJ/$RailTruckP if (mode2 == 6 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$BargeTruckP if (mode2 == 8 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$AirTruckP if (mode2 == 11 & BAUmode == 4)
* Rail
qui: replace tonmilesADJ = tonmilesADJ/$TruckRailP if (mode2 == 4 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$BargeRailP if (mode2 == 8 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$AirRailP if (mode2 == 11 & BAUmode == 6)
* Water
qui: replace tonmilesADJ = tonmilesADJ/$TruckWaterP if (mode2 == 4 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$RailWaterP if (mode2 == 6 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$AirWaterP if (mode2 == 11 & BAUmode == 8)
* Air
qui: replace tonmilesADJ = tonmilesADJ/$TruckAirP if (mode2 == 4 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$RailAirP if (mode2 == 6 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$BargeAirP if (mode2 == 8 & BAUmode == 11)

}

* Initialize variables
gen FuelBAU = .
gen FuelCAFE = .
gen MMTco2BAU = .
gen MMTco2CAFE = .
* Calculate Ton-miles, Fuel consumption and Carbon emissions - BAU
qui: gen TonMilesBAU = PBAU*tonmilesADJ
qui: replace FuelBAU = PBAU*$gptmTruck*tonmilesADJ if mode2 == 4
qui: replace FuelBAU = PBAU*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelBAU = PBAU*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelBAU = PBAU*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 4
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 6
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 8
qui: replace MMTco2BAU  = $fuelMTpgJet*FuelBAU if mode2 == 11
*** Calculate Ton-miles, Fuel consumption and Carbon emissions - CAFE
qui: gen TonMilesCAFE = PCAFE*tonmilesADJ
qui: replace FuelCAFE = PCAFE*$gptmTruck*tonmilesADJ*(1-$CAFE) if mode2 == 4
qui: replace FuelCAFE = PCAFE*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelCAFE = PCAFE*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelCAFE = PCAFE*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 4
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 6
qui: replace MMTco2CAFE = $fuelMTpg*FuelCAFE if mode2 == 8
qui: replace MMTco2CAFE  = $fuelMTpgJet*FuelCAFE if mode2 == 11
* Aggregate by mode
qui: collapse (first) sctg_txt (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE [iweight = wgt_factor], by(mode2) 

* Save file
qui: save "Tables/TruckCAFE`g'CAFE.dta", replace
}

*** 2. Models with truck, rail and barge as main modes
** Non-temp controlled
* Loop over commodities
#delimit ;
foreach g in "Grain" "Coal" "Gravel" 
"Transportation Equipment, not elsewhere classified" "Waste and Scrap" { ;
#delimit cr

* Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 | mode == 8

if $moreRivers == 1 {
replace Miss = 1 if orig_state_txt == "AL"
replace Miss = 1 if orig_state_txt == "PA"
replace Miss = 1 if orig_state_txt == "WV"
replace Miss = 1 if orig_state_txt == "OH"
replace Miss = 1 if orig_state_txt == "IN"
}

* Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt8 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmBarge*Pfuel*tonmiles/10^6 if mode2 == 8
* Create variables for predictions
qui: gen PBAU = 0
qui: gen PCAFE = 0

*** Estimate Mixed Logit and Save Parameter estimates
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles Miss) base(6) cluster(routeid)
qui: matrix B = e(b)
qui: eststo Results
* Table for parameter estimates
estout Results using "Tables/`g'Results_ASC.txt", mlabels("`g'") cells(b(fmt(%9.3f)) t(par fmt(%9.3f)) ) style(tab) stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace

** Calculate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAU = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7])) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAU = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7])) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PBAU = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7] )) if mode2 == 8;
#delimit cr

* Truck CAFE Scenario
* Truck
#delimit ;
qui: replace PCAFE = exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1] ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7])) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PCAFE = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7])) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PCAFE = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1] ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7] )) if mode2 == 8;
#delimit cr

*** Truck CAFE Calculations
gen tonmilesADJ = tonmiles
** Find most probably mode
* BAU
sort id
qui: by id: egen tempBAU = max(PBAU)
qui: gen maxprobBAU = tempBAU == PBAU
* CAFE
qui: by id: egen tempCAFE = max(PCAFE)
qui: gen maxprobCAFE = tempCAFE == PCAFE
** Correction for distance penalties
if $distPenalty == 1 {
sort shipmt_id mode2
qui: gen temp1 = 0
qui: replace temp1 = mode2 if maxprobBAU== 1
qui: by shipmt_id: egen BAUmode = max(temp1)
qui: drop temp1
** Adjust ton-miles should shipment switch mode
* Truck
qui: replace tonmilesADJ = tonmilesADJ/$RailTruckP if (mode2 == 6 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$BargeTruckP if (mode2 == 8 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$AirTruckP if (mode2 == 11 & BAUmode == 4)
* Rail
qui: replace tonmilesADJ = tonmilesADJ/$TruckRailP if (mode2 == 4 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$BargeRailP if (mode2 == 8 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$AirRailP if (mode2 == 11 & BAUmode == 6)
* Water
qui: replace tonmilesADJ = tonmilesADJ/$TruckWaterP if (mode2 == 4 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$RailWaterP if (mode2 == 6 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$AirWaterP if (mode2 == 11 & BAUmode == 8)
* Air
qui: replace tonmilesADJ = tonmilesADJ/$TruckAirP if (mode2 == 4 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$RailAirP if (mode2 == 6 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$BargeAirP if (mode2 == 8 & BAUmode == 11)

}

* Initialize variables
gen FuelBAU = .
gen FuelCAFE = .
gen MMTco2BAU = .
gen MMTco2CAFE = .
* Calculate Ton-miles, Fuel consumption and Carbon emissions - BAU
qui: gen TonMilesBAU = PBAU*tonmilesADJ
qui: replace FuelBAU = PBAU*$gptmTruck*tonmilesADJ if mode2 == 4
qui: replace FuelBAU = PBAU*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelBAU = PBAU*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelBAU = PBAU*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 4
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 6
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 8
qui: replace MMTco2BAU  = $fuelMTpgJet*FuelBAU if mode2 == 11
*** Calculate Ton-miles, Fuel consumption and Carbon emissions - CAFE
qui: gen TonMilesCAFE = PCAFE*tonmilesADJ
qui: replace FuelCAFE = PCAFE*$gptmTruck*tonmilesADJ*(1-$CAFE) if mode2 == 4
qui: replace FuelCAFE = PCAFE*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelCAFE = PCAFE*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelCAFE = PCAFE*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 4
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 6
qui: replace MMTco2CAFE = $fuelMTpg*FuelCAFE if mode2 == 8
qui: replace MMTco2CAFE  = $fuelMTpgJet*FuelCAFE if mode2 == 11
* Aggregate by mode
qui: collapse (first) sctg_txt (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE [iweight = wgt_factor], by(mode2) 

* Save file
qui: save "Tables/TruckCAFE`g'CAFE.dta", replace
}

*** 3. Models with truck and rail as main modes
** Temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Alcohol" "Animal Feed" "Milled Grain" "Metallic Ores"
"Other Chemical Products" "Other Prepared Foodstuffs"  { ;
#delimit cr

* Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
* Create variables for predictions
qui: gen PBAU = 0
qui: gen PCAFE = 0

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles TempCont) base(6) cluster(routeid)
qui: matrix B = e(b)
qui: eststo Results
* Table for parameter estimates
estout Results using "Tables/`g'Results_ASC.txt", mlabels("`g'") cells(b(fmt(%9.3f)) t(par fmt(%9.3f)) ) style(tab) stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace

** Calculate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAU = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAU = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 6;
#delimit cr

** Truck CAFE Scenario
* Truck
#delimit ;
qui: replace PCAFE = exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PCAFE = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 6;
#delimit cr

*** Truck CAFE Calculations
gen tonmilesADJ = tonmiles
** Find most probably mode
* BAU
sort id
qui: by id: egen tempBAU = max(PBAU)
qui: gen maxprobBAU = tempBAU == PBAU
* CAFE
qui: by id: egen tempCAFE = max(PCAFE)
qui: gen maxprobCAFE = tempCAFE == PCAFE
** Correction for distance penalties
if $distPenalty == 1 {
sort shipmt_id mode2
qui: gen temp1 = 0
qui: replace temp1 = mode2 if maxprobBAU== 1
qui: by shipmt_id: egen BAUmode = max(temp1)
qui: drop temp1
** Adjust ton-miles should shipment switch mode
* Truck
qui: replace tonmilesADJ = tonmilesADJ/$RailTruckP if (mode2 == 6 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$BargeTruckP if (mode2 == 8 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$AirTruckP if (mode2 == 11 & BAUmode == 4)
* Rail
qui: replace tonmilesADJ = tonmilesADJ/$TruckRailP if (mode2 == 4 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$BargeRailP if (mode2 == 8 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$AirRailP if (mode2 == 11 & BAUmode == 6)
* Water
qui: replace tonmilesADJ = tonmilesADJ/$TruckWaterP if (mode2 == 4 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$RailWaterP if (mode2 == 6 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$AirWaterP if (mode2 == 11 & BAUmode == 8)
* Air
qui: replace tonmilesADJ = tonmilesADJ/$TruckAirP if (mode2 == 4 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$RailAirP if (mode2 == 6 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$BargeAirP if (mode2 == 8 & BAUmode == 11)

}

* Initialize variables
gen FuelBAU = .
gen FuelCAFE = .
gen MMTco2BAU = .
gen MMTco2CAFE = .
* Calculate Ton-miles, Fuel consumption and Carbon emissions - BAU
qui: gen TonMilesBAU = PBAU*tonmilesADJ
qui: replace FuelBAU = PBAU*$gptmTruck*tonmilesADJ if mode2 == 4
qui: replace FuelBAU = PBAU*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelBAU = PBAU*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelBAU = PBAU*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 4
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 6
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 8
qui: replace MMTco2BAU  = $fuelMTpgJet*FuelBAU if mode2 == 11
*** Calculate Ton-miles, Fuel consumption and Carbon emissions - CAFE
qui: gen TonMilesCAFE = PCAFE*tonmilesADJ
qui: replace FuelCAFE = PCAFE*$gptmTruck*tonmilesADJ*(1-$CAFE) if mode2 == 4
qui: replace FuelCAFE = PCAFE*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelCAFE = PCAFE*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelCAFE = PCAFE*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 4
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 6
qui: replace MMTco2CAFE = $fuelMTpg*FuelCAFE if mode2 == 8
qui: replace MMTco2CAFE  = $fuelMTpgJet*FuelCAFE if mode2 == 11
* Aggregate by mode
qui: collapse (first) sctg_txt (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE [iweight = wgt_factor], by(mode2) 

* Save file
qui: save "Tables/TruckCAFE`g'CAFE.dta", replace
}

*** 4. Models with truck and rail as main modes
** Non temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Articles of Base Metal" "Logs and Other Wood in the Rough" "Non-Metallic Mineral Products"
"Paper" "Plastics and Rubber" "Primary Base Metal" "Pulp, Newsprint, Paper, and Paperboard" 
"Sand" "Textiles" "Vehicles" "Wood Products" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
* Create variables for predictions
qui: gen PBAU = 0
qui: gen PCAFE = 0

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(6) cluster(routeid)
qui: matrix B = e(b)
qui: eststo Results
* Table for parameter estimates
estout Results using "Tables/`g'Results_ASC.txt", mlabels("`g'") cells(b(fmt(%9.3f)) t(par fmt(%9.3f)) ) style(tab) stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAU = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAU = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 6;
#delimit cr

** Truck CAFE Scenario
* Truck
#delimit ;
qui: replace PCAFE = exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PCAFE = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 6;
#delimit cr

*** Truck CAFE Calculations
gen tonmilesADJ = tonmiles
** Find most probably mode
* BAU
sort id
qui: by id: egen tempBAU = max(PBAU)
qui: gen maxprobBAU = tempBAU == PBAU
* CAFE
qui: by id: egen tempCAFE = max(PCAFE)
qui: gen maxprobCAFE = tempCAFE == PCAFE
** Correction for distance penalties
if $distPenalty == 1 {
sort shipmt_id mode2
qui: gen temp1 = 0
qui: replace temp1 = mode2 if maxprobBAU== 1
qui: by shipmt_id: egen BAUmode = max(temp1)
qui: drop temp1
** Adjust ton-miles should shipment switch mode
* Truck
qui: replace tonmilesADJ = tonmilesADJ/$RailTruckP if (mode2 == 6 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$BargeTruckP if (mode2 == 8 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$AirTruckP if (mode2 == 11 & BAUmode == 4)
* Rail
qui: replace tonmilesADJ = tonmilesADJ/$TruckRailP if (mode2 == 4 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$BargeRailP if (mode2 == 8 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$AirRailP if (mode2 == 11 & BAUmode == 6)
* Water
qui: replace tonmilesADJ = tonmilesADJ/$TruckWaterP if (mode2 == 4 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$RailWaterP if (mode2 == 6 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$AirWaterP if (mode2 == 11 & BAUmode == 8)
* Air
qui: replace tonmilesADJ = tonmilesADJ/$TruckAirP if (mode2 == 4 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$RailAirP if (mode2 == 6 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$BargeAirP if (mode2 == 8 & BAUmode == 11)

}

* Initialize variables
gen FuelBAU = .
gen FuelCAFE = .
gen MMTco2BAU = .
gen MMTco2CAFE = .
* Calculate Ton-miles, Fuel consumption and Carbon emissions - BAU
qui: gen TonMilesBAU = PBAU*tonmilesADJ
qui: replace FuelBAU = PBAU*$gptmTruck*tonmilesADJ if mode2 == 4
qui: replace FuelBAU = PBAU*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelBAU = PBAU*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelBAU = PBAU*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 4
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 6
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 8
qui: replace MMTco2BAU  = $fuelMTpgJet*FuelBAU if mode2 == 11
*** Calculate Ton-miles, Fuel consumption and Carbon emissions - CAFE
qui: gen TonMilesCAFE = PCAFE*tonmilesADJ
qui: replace FuelCAFE = PCAFE*$gptmTruck*tonmilesADJ*(1-$CAFE) if mode2 == 4
qui: replace FuelCAFE = PCAFE*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelCAFE = PCAFE*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelCAFE = PCAFE*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 4
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 6
qui: replace MMTco2CAFE = $fuelMTpg*FuelCAFE if mode2 == 8
qui: replace MMTco2CAFE  = $fuelMTpgJet*FuelCAFE if mode2 == 11
* Aggregate by mode
qui: collapse (first) sctg_txt (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE [iweight = wgt_factor], by(mode2) 
 
* Save file
qui: save "Tables/TruckCAFE`g'CAFE.dta", replace
}

*** 5. Models with Air and Truck
** Temp controlled
* Loop over commodities
#delimit ;
foreach g in  "Animals" "Pharmaceuticals" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 11
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
* Create variables for predictions
qui: gen PBAU = 0
qui: gen PCAFE = 0

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles TempCont) base(11) cluster(routeid)
qui: matrix B = e(b)
qui: eststo Results
* Table for parameter estimates
estout Results using "Tables/`g'Results_ASC.txt", mlabels("`g'") cells(b(fmt(%9.3f)) t(par fmt(%9.3f)) ) style(tab) stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace

** Calculate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAU = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PBAU = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 11;
#delimit cr

** Truck CAFE Scenario
* Truck
#delimit ;
qui: replace PCAFE = exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PCAFE = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 11;
#delimit cr

*** Truck CAFE Calculations
gen tonmilesADJ = tonmiles
** Find most probably mode in each scenario
* BAU
sort id
qui: by id: egen tempBAU = max(PBAU)
qui: gen maxprobBAU = tempBAU == PBAU
* CAFE
qui: by id: egen tempCAFE = max(PCAFE)
qui: gen maxprobCAFE = tempCAFE == PCAFE
** Correction for distance penalties
if $distPenalty == 1 {
sort shipmt_id mode2
qui: gen temp1 = 0
qui: replace temp1 = mode2 if maxprobBAU== 1
qui: by shipmt_id: egen BAUmode = max(temp1)
qui: drop temp1
** Adjust ton-miles should shipment switch mode
* Truck
qui: replace tonmilesADJ = tonmilesADJ/$RailTruckP if (mode2 == 6 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$BargeTruckP if (mode2 == 8 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$AirTruckP if (mode2 == 11 & BAUmode == 4)
* Rail
qui: replace tonmilesADJ = tonmilesADJ/$TruckRailP if (mode2 == 4 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$BargeRailP if (mode2 == 8 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$AirRailP if (mode2 == 11 & BAUmode == 6)
* Water
qui: replace tonmilesADJ = tonmilesADJ/$TruckWaterP if (mode2 == 4 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$RailWaterP if (mode2 == 6 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$AirWaterP if (mode2 == 11 & BAUmode == 8)
* Air
qui: replace tonmilesADJ = tonmilesADJ/$TruckAirP if (mode2 == 4 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$RailAirP if (mode2 == 6 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$BargeAirP if (mode2 == 8 & BAUmode == 11)

}

* Initialize variables
gen FuelBAU = .
gen FuelCAFE = .
gen MMTco2BAU = .
gen MMTco2CAFE = .
* Calculate Ton-miles, Fuel consumption and Carbon emissions - BAU
qui: gen TonMilesBAU = PBAU*tonmilesADJ
qui: replace FuelBAU = PBAU*$gptmTruck*tonmilesADJ if mode2 == 4
qui: replace FuelBAU = PBAU*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelBAU = PBAU*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelBAU = PBAU*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 4
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 6
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 8
qui: replace MMTco2BAU  = $fuelMTpgJet*FuelBAU if mode2 == 11
*** Calculate Ton-miles, Fuel consumption and Carbon emissions - CAFE
qui: gen TonMilesCAFE = PCAFE*tonmilesADJ
qui: replace FuelCAFE = PCAFE*$gptmTruck*tonmilesADJ*(1-$CAFE) if mode2 == 4
qui: replace FuelCAFE = PCAFE*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelCAFE = PCAFE*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelCAFE = PCAFE*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 4
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 6
qui: replace MMTco2CAFE = $fuelMTpg*FuelCAFE if mode2 == 8
qui: replace MMTco2CAFE  = $fuelMTpgJet*FuelCAFE if mode2 == 11
* Aggregate by mode
qui: collapse (first) sctg_txt (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE [iweight = wgt_factor], by(mode2) 

* Save file
qui: save "Tables/TruckCAFE`g'CAFE.dta", replace
}


*** 6. Models with air and truck 
** Not temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Precision Instruments" "Printed Products" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 11
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
* Create variables for predictions
qui: gen PBAU = 0
qui: gen PCAFE = 0

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(11) cluster(routeid)
qui: matrix B = e(b)
qui: eststo Results
* Table for parameter estimates
estout Results using "Tables/`g'Results_ASC.txt", mlabels("`g'") cells(b(fmt(%9.3f)) t(par fmt(%9.3f)) ) style(tab) stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAU = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PBAU = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 11;
#delimit cr

** Truck CAFE Scenario
* Truck
#delimit ;
qui: replace PCAFE = exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PCAFE = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1])) if mode2 == 11;
#delimit cr

*** Truck CAFE Calculations
gen tonmilesADJ = tonmiles
** Find most probably mode
* BAU
sort id
qui: by id: egen tempBAU = max(PBAU)
qui: gen maxprobBAU = tempBAU == PBAU
* CAFE
qui: by id: egen tempCAFE = max(PCAFE)
qui: gen maxprobCAFE = tempCAFE == PCAFE
** Correction for distance penalties
if $distPenalty == 1 {
sort shipmt_id mode2
qui: gen temp1 = 0
qui: replace temp1 = mode2 if maxprobBAU== 1
qui: by shipmt_id: egen BAUmode = max(temp1)
qui: drop temp1
** Adjust ton-miles should shipment switch mode
* Truck
qui: replace tonmilesADJ = tonmilesADJ/$RailTruckP if (mode2 == 6 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$BargeTruckP if (mode2 == 8 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$AirTruckP if (mode2 == 11 & BAUmode == 4)
* Rail
qui: replace tonmilesADJ = tonmilesADJ/$TruckRailP if (mode2 == 4 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$BargeRailP if (mode2 == 8 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$AirRailP if (mode2 == 11 & BAUmode == 6)
* Water
qui: replace tonmilesADJ = tonmilesADJ/$TruckWaterP if (mode2 == 4 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$RailWaterP if (mode2 == 6 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$AirWaterP if (mode2 == 11 & BAUmode == 8)
* Air
qui: replace tonmilesADJ = tonmilesADJ/$TruckAirP if (mode2 == 4 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$RailAirP if (mode2 == 6 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$BargeAirP if (mode2 == 8 & BAUmode == 11)

}

* Initialize variables
gen FuelBAU = .
gen FuelCAFE = .
gen MMTco2BAU = .
gen MMTco2CAFE = .
* Calculate Ton-miles, Fuel consumption and Carbon emissions - BAU
qui: gen TonMilesBAU = PBAU*tonmilesADJ
qui: replace FuelBAU = PBAU*$gptmTruck*tonmilesADJ if mode2 == 4
qui: replace FuelBAU = PBAU*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelBAU = PBAU*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelBAU = PBAU*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 4
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 6
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 8
qui: replace MMTco2BAU  = $fuelMTpgJet*FuelBAU if mode2 == 11
*** Calculate Ton-miles, Fuel consumption and Carbon emissions - CAFE
qui: gen TonMilesCAFE = PCAFE*tonmilesADJ
qui: replace FuelCAFE = PCAFE*$gptmTruck*tonmilesADJ*(1-$CAFE) if mode2 == 4
qui: replace FuelCAFE = PCAFE*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelCAFE = PCAFE*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelCAFE = PCAFE*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 4
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 6
qui: replace MMTco2CAFE = $fuelMTpg*FuelCAFE if mode2 == 8
qui: replace MMTco2CAFE  = $fuelMTpgJet*FuelCAFE if mode2 == 11
* Aggregate by mode
qui: collapse (first) sctg_txt (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE [iweight = wgt_factor], by(mode2) 

* Save file
qui: save "Tables/TruckCAFE`g'CAFE.dta", replace
}

*** 7. Models with truck, rail and air as main modes
** Temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Mixed Freight" { ;
#delimit cr
** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 | mode == 11
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
* Create variables for predictions
qui: gen PBAU = 0
qui: gen PCAFE = 0

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles TempCont) base(6) cluster(routeid)
qui: matrix B = e(b)
qui: eststo Results
* Table for parameter estimates
estout Results using "Tables/`g'Results_ASC.txt", mlabels("`g'") cells(b(fmt(%9.3f)) t(par fmt(%9.3f)) ) style(tab) stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace

** Calculate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAU = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7]) ) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAU = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7]) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PBAU = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7]) ) if mode2 == 11;
#delimit cr

** Truck CAFE Scenario
* Truck  
#delimit ;
qui: replace PCAFE = exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7]) ) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PCAFE = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7]) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PCAFE = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7]) ) if mode2 == 11;
#delimit cr

*** Truck CAFE Calculations
gen tonmilesADJ = tonmiles
** Find most probably mode
* BAU
sort id
qui: by id: egen tempBAU = max(PBAU)
qui: gen maxprobBAU = tempBAU == PBAU
* CAFE
qui: by id: egen tempCAFE = max(PCAFE)
qui: gen maxprobCAFE = tempCAFE == PCAFE
** Correction for distance penalties
if $distPenalty == 1 {
sort shipmt_id mode2
qui: gen temp1 = 0
qui: replace temp1 = mode2 if maxprobBAU== 1
qui: by shipmt_id: egen BAUmode = max(temp1)
qui: drop temp1
** Adjust ton-miles should shipment switch mode
* Truck
qui: replace tonmilesADJ = tonmilesADJ/$RailTruckP if (mode2 == 6 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$BargeTruckP if (mode2 == 8 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$AirTruckP if (mode2 == 11 & BAUmode == 4)
* Rail
qui: replace tonmilesADJ = tonmilesADJ/$TruckRailP if (mode2 == 4 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$BargeRailP if (mode2 == 8 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$AirRailP if (mode2 == 11 & BAUmode == 6)
* Water
qui: replace tonmilesADJ = tonmilesADJ/$TruckWaterP if (mode2 == 4 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$RailWaterP if (mode2 == 6 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$AirWaterP if (mode2 == 11 & BAUmode == 8)
* Air
qui: replace tonmilesADJ = tonmilesADJ/$TruckAirP if (mode2 == 4 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$RailAirP if (mode2 == 6 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$BargeAirP if (mode2 == 8 & BAUmode == 11)

}

* Initialize variables
gen FuelBAU = .
gen FuelCAFE = .
gen MMTco2BAU = .
gen MMTco2CAFE = .
* Calculate Ton-miles, Fuel consumption and Carbon emissions - BAU
qui: gen TonMilesBAU = PBAU*tonmilesADJ
qui: replace FuelBAU = PBAU*$gptmTruck*tonmilesADJ if mode2 == 4
qui: replace FuelBAU = PBAU*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelBAU = PBAU*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelBAU = PBAU*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 4
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 6
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 8
qui: replace MMTco2BAU  = $fuelMTpgJet*FuelBAU if mode2 == 11
*** Calculate Ton-miles, Fuel consumption and Carbon emissions - CAFE
qui: gen TonMilesCAFE = PCAFE*tonmilesADJ
qui: replace FuelCAFE = PCAFE*$gptmTruck*tonmilesADJ*(1-$CAFE) if mode2 == 4
qui: replace FuelCAFE = PCAFE*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelCAFE = PCAFE*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelCAFE = PCAFE*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 4
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 6
qui: replace MMTco2CAFE = $fuelMTpg*FuelCAFE if mode2 == 8
qui: replace MMTco2CAFE  = $fuelMTpgJet*FuelCAFE if mode2 == 11
* Aggregate by mode
qui: collapse (first) sctg_txt (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE [iweight = wgt_factor], by(mode2) 

* Save file
qui: save "Tables/TruckCAFE`g'CAFE.dta", replace
}

*** 8. Models with truck, rail and air as main modes
** Not temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Machinery" "Miscellaneous Manufactured Products" { ;
#delimit cr
** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 | mode == 11
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
* Create variables for predictions
qui: gen PBAU = 0
qui: gen PCAFE = 0

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(6) cluster(routeid)
qui: matrix B = e(b)
qui: eststo Results
* Table for parameter estimates
estout Results using "Tables/`g'Results_ASC.txt", mlabels("`g'") cells(b(fmt(%9.3f)) t(par fmt(%9.3f)) ) style(tab) stats(N r2_p, fmt(%8.0f %9.6f)) starlevels(* 0.1 ** 0.05 *** 0.01) varlabels(_cons Constant) replace

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAU = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5]) ) if mode2 == 4;
#delimiti cr
* Rail
#delimit ;
qui: replace PBAU = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5]) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PBAU = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5])  /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5])  ) if mode2 == 11;
#delimit cr

** Truck CAFE Scenario  
* Truck
#delimit ;
qui: replace PCAFE = exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])/
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5]) ) if mode2 == 4;
#delimiti cr
* Rail
#delimit ;
qui: replace PCAFE = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]) /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5]) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PCAFE = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5])  /
(exp((1-$CAFE)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3])
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1])
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5])  ) if mode2 == 11;
#delimit cr

*** Truck CAFE Calculations
gen tonmilesADJ = tonmiles
** Find most probably mode in each scenario
* BAU
sort id
qui: by id: egen tempBAU = max(PBAU)
qui: gen maxprobBAU = tempBAU == PBAU
* CAFE
qui: by id: egen tempCAFE = max(PCAFE)
qui: gen maxprobCAFE = tempCAFE == PCAFE
** Correction for distance penalties
if $distPenalty == 1 {
sort shipmt_id mode2
qui: gen temp1 = 0
qui: replace temp1 = mode2 if maxprobBAU== 1
qui: by shipmt_id: egen BAUmode = max(temp1)
qui: drop temp1
** Adjust ton-miles should shipment switch mode
* Truck
qui: replace tonmilesADJ = tonmilesADJ/$RailTruckP if (mode2 == 6 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$BargeTruckP if (mode2 == 8 & BAUmode == 4)
qui: replace tonmilesADJ = tonmilesADJ/$AirTruckP if (mode2 == 11 & BAUmode == 4)
* Rail
qui: replace tonmilesADJ = tonmilesADJ/$TruckRailP if (mode2 == 4 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$BargeRailP if (mode2 == 8 & BAUmode == 6)
qui: replace tonmilesADJ = tonmilesADJ/$AirRailP if (mode2 == 11 & BAUmode == 6)
* Water
qui: replace tonmilesADJ = tonmilesADJ/$TruckWaterP if (mode2 == 4 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$RailWaterP if (mode2 == 6 & BAUmode == 8)
qui: replace tonmilesADJ = tonmilesADJ/$AirWaterP if (mode2 == 11 & BAUmode == 8)
* Air
qui: replace tonmilesADJ = tonmilesADJ/$TruckAirP if (mode2 == 4 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$RailAirP if (mode2 == 6 & BAUmode == 11)
qui: replace tonmilesADJ = tonmilesADJ/$BargeAirP if (mode2 == 8 & BAUmode == 11)

}

* Initialize variables
gen FuelBAU = .
gen FuelCAFE = .
gen MMTco2BAU = .
gen MMTco2CAFE = .
* Calculate Ton-miles, Fuel consumption and Carbon emissions - BAU
qui: gen TonMilesBAU = PBAU*tonmilesADJ
qui: replace FuelBAU = PBAU*$gptmTruck*tonmilesADJ if mode2 == 4
qui: replace FuelBAU = PBAU*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelBAU = PBAU*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelBAU = PBAU*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 4
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 6
qui: replace MMTco2BAU  = $fuelMTpg*FuelBAU if mode2 == 8
qui: replace MMTco2BAU  = $fuelMTpgJet*FuelBAU if mode2 == 11
*** Calculate Ton-miles, Fuel consumption and Carbon emissions - CAFE
qui: gen TonMilesCAFE = PCAFE*tonmilesADJ
qui: replace FuelCAFE = PCAFE*$gptmTruck*tonmilesADJ*(1-$CAFE) if mode2 == 4
qui: replace FuelCAFE = PCAFE*$gptmRail*tonmilesADJ if mode2 == 6
qui: replace FuelCAFE = PCAFE*$gptmBarge*tonmilesADJ if mode2 == 8
qui: replace FuelCAFE = PCAFE*$gptmAir*tonmilesADJ if mode2 == 11
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 4
qui: replace MMTco2CAFE  = $fuelMTpg*FuelCAFE if mode2 == 6
qui: replace MMTco2CAFE = $fuelMTpg*FuelCAFE if mode2 == 8
qui: replace MMTco2CAFE  = $fuelMTpgJet*FuelCAFE if mode2 == 11
* Aggregate by mode
qui: collapse (first) sctg_txt (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE [iweight = wgt_factor], by(mode2) 

* Save file
qui: save "Tables/TruckCAFE`g'CAFE.dta", replace
}

*** Table 1
* Combine CAFE files
drop _all
* 
#delimit ;
foreach g in "Agricultural Products" "Basic Chemicals" "Fertilizers" "Other Coal and Petroleum" 
"Coal" "Grain" "Gravel" "Transportation Equipment, not elsewhere classified" 
"Waste and Scrap" "Alcohol" "Animal Feed" "Milled Grain" 
"Other Chemical Products" "Other Prepared Foodstuffs" "Articles of Base Metal" 
"Logs and Other Wood in the Rough" "Non-Metallic Mineral Products" "Metallic Ores"
"Paper" "Plastics and Rubber" "Primary Base Metal" "Pulp, Newsprint, Paper, and Paperboard" 
"Sand" "Textiles" "Vehicles" "Wood Products" "Animals" "Pharmaceuticals" "Precision Instruments" 
"Printed Products" "Mixed Freight" "Machinery" "Miscellaneous Manufactured Products" { ;
#delimit cr
append using "Tables/TruckCAFE`g'CAFE.dta"
}
replace TonMilesBAU = TonMilesBAU/10^9
replace TonMilesCAFE = TonMilesCAFE/10^9
replace MMTco2BAU = MMTco2BAU/10^6
replace MMTco2CAFE = MMTco2CAFE/10^6
replace FuelBAU = FuelBAU/10^6
replace FuelCAFE = FuelCAFE/10^6
sort sctg_txt mode2
save "Tables/TruckCAFECombined.dta", replace
* Make duplicate file to model fit comparison
save "Tables/TruckCAFECombined2.dta", replace
* Make aggregate emissions and fuel file
use "Tables/TruckCAFECombined.dta", clear
sort mode2
collapse (sum) TonMilesBAU MMTco2BAU FuelBAU TonMilesCAFE MMTco2CAFE FuelCAFE, by(mode2)
save "Tables/TruckCAFECombinedTotals.dta", replace
gsort - mode2
outsheet using "Tables/FuelCO2AllCAFE.csv", comma replace
* Analyze model "fit" by comparing mode shares by good
use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
sort mode
collapse (sum) tonmiles [iweight = wgt_factor], by(mode sctg_txt)
 * Rename mode variable to do merge
rename mode mode2
sort sctg_txt mode2
merge 1:1 sctg_txt mode2 using "Tables/TruckCAFECombined2.dta", 
* We don't match the "minor" modes that were previously dropped for each SCTG, so drop them here
keep if _merge == 3
drop _merge
keep mode2 sctg_txt tonmiles TonMilesBAU
rename mode2 mode
replace tonmiles = tonmiles/10^9
* Create table
fillin mode sctg
order sctg_txt mode tonmiles TonMilesBAU
sort sctg_txt
replace tonmiles = 0 if tonmiles == .
replace TonMilesBAU = 0 if TonMilesBAU == .
drop _fillin
reshape wide tonmiles TonMilesBAU, i(sctg_txt) j(mode)
*order sctg_txt tonmiles* predtonmiles*
outsheet using "Tables/MLFit.csv", comma replace

*** Rebound calculations
use "Tables/TruckCAFECombined.dta", clear
sort sctg_txt
by sctg_txt: egen BAUemissions = sum(MMTco2BAU)
by sctg_txt: egen CAFEemissions = sum(MMTco2CAFE)
gen temp = MMTco2BAU
replace temp = (1-$CAFE)*temp if mode2 == 4
by sctg_txt: egen PredCAFEemissions = sum(temp)
drop temp
gen PredAbate = (BAUemissions-PredCAFEemissions)
gen ActAbate =  (BAUemissions-CAFEemissions)
gen Ratio = (ActAbate)/PredAbate
gen Rebound = 1-(ActAbate)/PredAbate
gen TonMileInc = (TonMilesCAFE-TonMilesBAU)/TonMilesBAU
order Ratio PredAbate ActAbate TonMileInc sctg_txt
* Clean up some names
replace sctg_txt = "Transportation Equipment, nec" if sctg_txt == "Transportation Equipment, not elsewhere classified"
replace sctg_txt = "Pulp, Paper, Newsprint" if sctg_txt == "Pulp, Newsprint, Paper, and Paperboard"
gen shortfall = CAFEemissions- PredCAFEemissions
keep Ratio Rebound sctg_txt BAUemissions shortfall PredAbate ActAbate
duplicates drop
gsort - Ratio
gen order = _n

*** Figure 4
** Version used in paper - 4 colors version
* The trick is to create duplicate observations so weights remain the same, but 3/4 of obs within each group have a missing
* x or y coordinate
* Create groups by available transportation modes
gen group = . 
replace group = 1 if sctg_txt == "Agricultural Products" 
replace group = 1 if sctg_txt == "Basic Chemicals" 
replace group = 1 if sctg_txt == "Fertilizers"
replace group = 1 if sctg_txt == "Other Coal and Petroleum"
replace group = 1 if sctg_txt == "Coal"
replace group = 1 if sctg_txt == "Grain"  
replace group = 1 if sctg_txt == "Gravel"
replace group = 1 if sctg_txt == "Transportation Equipment, nec"
replace group = 1 if sctg_txt == "Waste and Scrap" 
replace group = 2 if sctg_txt == "Alcohol"
replace group = 2 if sctg_txt == "Animal Feed"
replace group = 2 if sctg_txt == "Milled Grain"
replace group = 2 if sctg_txt == "Metallic Ores"
replace group = 2 if sctg_txt == "Other Chemical Products"
replace group = 2 if sctg_txt == "Other Prepared Foodstuffs"
replace group = 2 if sctg_txt == "Articles of Base Metal" 
replace group = 2 if sctg_txt == "Logs and Other Wood in the Rough"
replace group = 2 if sctg_txt == "Non-Metallic Mineral Products"
replace group = 2 if sctg_txt == "Paper"
replace group = 2 if sctg_txt == "Plastics and Rubber"
replace group = 2 if sctg_txt == "Primary Base Metal"
replace group = 2 if sctg_txt == "Pulp, Paper, Newsprint"
replace group = 2 if sctg_txt == "Sand"
replace group = 2 if sctg_txt == "Textiles"
replace group = 2 if sctg_txt == "Vehicles"
replace group = 2 if sctg_txt == "Wood Products"
replace group = 3 if sctg_txt == "Animals"
replace group = 3 if sctg_txt == "Pharmaceuticals"
replace group = 3 if sctg_txt == "Precision Instruments"
replace group = 3 if sctg_txt == "Printed Products"  
replace group = 4 if sctg_txt == "Mixed Freight" 
replace group = 4 if sctg_txt == "Machinery"
replace group = 4 if sctg_txt == "Miscellaneous Manufactured Products" 
* Create duplicate observations
expand 4, gen(created)
* Create missing x-var
replace Rebound = . if created == 1
* Create new groups
sort order created
replace group = 1 in 2
replace group = 2 in 3
replace group = 4 in 4
replace group = 1 in 6
replace group = 2 in 7
replace group = 4 in 8
replace group = 1 in 10
replace group = 2 in 11
replace group = 4 in 12
replace group = 1 in 14
replace group = 2 in 15
replace group = 4 in 16
replace group = 1 in 18
replace group = 2 in 19
replace group = 3 in 20
replace group = 2 in 22
replace group = 3 in 23
replace group = 4 in 24
replace group = 2 in 26
replace group = 3 in 27
replace group = 4 in 28
replace group = 1 in 30
replace group = 2 in 31
replace group = 3 in 32
replace group = 1 in 34
replace group = 3 in 35
replace group = 4 in 36
replace group = 1 in 38
replace group = 3 in 39
replace group = 4 in 40
replace group = 1 in 42
replace group = 3 in 43
replace group = 4 in 44
replace group = 1 in 46
replace group = 2 in 47
replace group = 3 in 48
replace group = 2 in 51
replace group = 3 in 52
replace group = 1 in 54
replace group = 3 in 55
replace group = 4 in 56
replace group = 1 in 58
replace group = 3 in 59
replace group = 4 in 60
replace group = 2 in 62
replace group = 3 in 63
replace group = 4 in 64
replace group = 2 in 66
replace group = 3 in 67
replace group = 4 in 68
replace group = 2 in 70
replace group = 3 in 71
replace group = 4 in 72
replace group = 1 in 74
replace group = 3 in 75
replace group = 4 in 76
replace group = 1 in 78
replace group = 3 in 79
replace group = 4 in 80
replace group = 1 in 82
replace group = 3 in 83
replace group = 4 in 84
replace group = 1 in 86
replace group = 3 in 87
replace group = 4 in 88
replace group = 1 in 90
replace group = 3 in 91
replace group = 4 in 92
replace group = 1 in 94
replace group = 3 in 95
replace group = 4 in 96
replace group = 1 in 98
replace group = 3 in 99
replace group = 4 in 100
replace group = 1 in 102
replace group = 3 in 103
replace group = 4 in 104
replace group = 2 in 106
replace group = 3 in 107
replace group = 4 in 108
replace group = 1 in 110
replace group = 3 in 111
replace group = 4 in 112
replace group = 1 in 114
replace group = 3 in 115
replace group = 4 in 116
replace group = 1 in 118
replace group = 3 in 119
replace group = 4 in 120
replace group = 2 in 122
replace group = 3 in 123
replace group = 4 in 124
replace group = 2 in 126
replace group = 3 in 127
replace group = 4 in 128
replace group = 1 in 130
replace group = 3 in 131
replace group = 4 in 132
* Draw graph
#delimit ;
twoway (scatter order Rebound [w= BAUemissions] if group == 1, xline(1) scheme(s1mono) mlcolor(gs6) msymbol(circle) mfcolor("68 1 84%50"))
(scatter order Rebound [w= BAUemissions] if group == 2, msymbol(circle) mlcolor(gs6) mfcolor("46 110 142%50"))
(scatter order Rebound [w= BAUemissions] if group == 4, msymbol(circle) mlcolor(gs6) mfcolor("45 178 125%50"))
(scatter order Rebound [w= BAUemissions] if group == 3, msymbol(circle) mlcolor(gs6) mfcolor("253 231 37%50") ytitle("") yscale(off)) 
(scatter order Rebound, ms(i) mlabel( sctg_txt ) xscale(range(-0.2 0.6)) xlabel(-0.2[0.2]0.6)
mlabposition(2) mlabangle(0) mlabsize(vsmall) xtitle("Substitution Rebound Effect")
yscale(range(1 33)) ylabel(1[1]33)
legend(row(1) order(1 "Barge, Rail, Truck" 2 "Rail, Truck" 3 "Air, Rail, Truck" 4 "Air, Truck") ));
#delimit cr
graph export "Figures/TruckCAFEBubble3_4color.pdf", replace

***** Elasticity Calculations
*** Assume a 1% increase in modal fuel intensity
*** Generate bar graphs for Rail-Truck and Air-Truck shipments

*** 1. Models with truck, rail and barge as main modes
** Temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Agricultural Products" "Basic Chemicals" 
"Fertilizers" "Other Coal and Petroleum" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 | mode == 8
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt8 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmBarge*Pfuel*tonmiles/10^6 if mode2 == 8
* Create some variables for error and predictions
qui: gen PBAUrand = 0
qui: gen PTruckrand_E = 0
qui: gen PRailrand_E = 0
qui: gen PBargerand_E = 0
qui: gen PAirrand_E = 0
qui: save "Data/MixedLogitTemp.dta", replace

*** Estimate Mixed Logit and Save Parameter estimates 
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles Miss TempCont) base(6)
qui: matrix B = e(b)

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAUrand = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAUrand = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PBAUrand = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 8;
#delimit cr

** Truck Elasticity Scenario
* Truck
#delimit ;
qui: replace PTruckrand_E = exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PTruckrand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PTruckrand_E = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 8;
#delimit cr

** Rail Elasticity Scenario
* Truck
#delimit ;
qui: replace PRailrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PRailrand_E = exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PRailrand_E = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 8;
#delimit cr

** Barge Elasticity Scenario
* Truck
#delimit ;
qui: replace PBargerand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp((1.01)*$gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBargerand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp((1.01)*$gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PBargerand_E = exp((1.01)*$gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + TempCont*B[1,4] + B[1,5] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp((1.01)*$gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,6] + Miss*B[1,7] + TempCont*B[1,8] + B[1,9]   )) if mode2 == 8;
#delimit cr
 
* Initialize variables
gen tonmilesBAU = .
gen tonmilesTruck_E = .
gen tonmilesRail_E = .
gen tonmilesBarge_E = .
gen tonmilesAir_E = .

* BAU
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 4
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 6
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 8
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 11
* Truck Elasticity Scenario
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 4
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 6
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 8
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 11
* Rail Elasticity Scenario
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 4
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 6
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 8
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 11
* Barge Elasticity Scenario
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 4
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 6
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 8
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 11
* Air Elasticity Scenario
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 4
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 6
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 8
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 11

qui: collapse (first) sctg_txt (sum) tonmilesBAU tonmilesTruck_E tonmilesRail_E tonmilesBarge_E tonmilesAir_E [iweight = wgt_factor], by(mode2) 

* Unit conversion
qui: replace tonmilesBAU = tonmilesBAU/10^9
qui: replace tonmilesTruck_E = tonmilesTruck_E/10^9
qui: replace tonmilesRail_E = tonmilesRail_E/10^9
qui: replace tonmilesBarge_E = tonmilesBarge_E/10^9
qui: replace tonmilesAir_E = tonmilesAir_E/10^9

qui: save "Tables/`g'Elast.dta", replace

}
* Create elasticity table for each commodity
#delimit ;
foreach g in "Agricultural Products" "Basic Chemicals" 
"Fertilizers" "Other Coal and Petroleum" { ;
#delimit cr
use "Tables/`g'Elast.dta", clear
collapse (mean) tonmiles* , by(sctg_txt mode2)
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PRailElast = ((tonmilesRail_E-tonmilesBAU)/tonmilesBAU)/.01
gen PBargeElast = ((tonmilesBarge_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PRailElast PBargeElast
xpose, clear varname
order _varname
replace _varname = "`g'" in 1
outsheet using "Tables/`g'ElastTable.csv", comma replace
}

*** 2. Models with truck, rail and barge as main modes
** Non-temp controlled
* Loop over commodities
#delimit ;
foreach g in "Coal" "Grain" "Gravel" 
"Transportation Equipment, not elsewhere classified" "Waste and Scrap" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 | mode == 8
* Reshape data for conditional logit implementation of mixed logit
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt8 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmBarge*Pfuel*tonmiles/10^6 if mode2 == 8
* Create some variables for error and predictions
qui: gen PBAUrand = 0
qui: gen PTruckrand_E = 0
qui: gen PRailrand_E = 0
qui: gen PBargerand_E = 0
qui: gen PAirrand_E = 0
qui: save "Data/MixedLogitTemp.dta", replace

*** Estimate Mixed Logit and Save Parameter estimates
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles Miss) base(6)
qui: matrix B = e(b)

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAUrand = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAUrand = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PBAUrand = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 8;
#delimit cr

*** Truck Elasticity Simulation
* Truck
#delimit ;
qui: replace PTruckrand_E = exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PTruckrand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PTruckrand_E = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 8;
#delimit cr

*** Rail Elasticity Simulation
* Truck
#delimit ;
qui: replace PRailrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PRailrand_E = exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PRailrand_E = exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp($gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 8;
#delimit cr
* Find most probably mode
qui: by id: egen tempRail_E = max(PRailrand_E)
qui: gen maxprobRail_E = tempRail_E == PRailrand_E
* Set mode to rail
qui: replace maxprobRail_E = 0 if PRailrand_E == . & mode2 == 4
qui: replace maxprobRail_E = 1 if PRailrand_E == . & mode2 == 6
qui: replace maxprobRail_E = 0 if PRailrand_E == . & mode2 == 8
qui: drop tempRail_E

*** Barge Elasticity Simulation
* Truck
#delimit ;
qui: replace PBargerand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp((1.01)*$gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBargerand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp((1.01)*$gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 6;
#delimit cr
* Barge
#delimit ;
qui: replace PBargerand_E = exp((1.01)*$gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + Miss*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) 
+ exp((1.01)*$gptmBarge*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + Miss*B[1,6] + B[1,7]   )) if mode2 == 8;
#delimit cr

* Initialize variables
gen tonmilesBAU = .
gen tonmilesTruck_E = .
gen tonmilesRail_E = .
gen tonmilesBarge_E = .
gen tonmilesAir_E = .

* BAU
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 4
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 6
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 8
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 11
* Truck Elasticity Scenario
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 4
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 6
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 8
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 11
* Rail Elasticity Scenario
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 4
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 6
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 8
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 11
* Barge Elasticity Scenario
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 4
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 6
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 8
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 11
* Air Elasticity Scenario
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 4
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 6
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 8
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 11

qui: collapse (first) sctg_txt (sum) tonmilesBAU tonmilesTruck_E tonmilesRail_E tonmilesBarge_E tonmilesAir_E [iweight = wgt_factor], by(mode2) 

* Unit conversion
qui: replace tonmilesBAU = tonmilesBAU/10^9
qui: replace tonmilesTruck_E = tonmilesTruck_E/10^9
qui: replace tonmilesRail_E = tonmilesRail_E/10^9
qui: replace tonmilesBarge_E = tonmilesBarge_E/10^9
qui: replace tonmilesAir_E = tonmilesAir_E/10^9

qui: save "Tables/`g'Elast.dta", replace

}
* Create elasticity table for each commodity
#delimit ;
foreach g in "Coal" "Grain" "Gravel" 
"Transportation Equipment, not elsewhere classified" "Waste and Scrap" { ;
#delimit cr
use "Tables/`g'Elast.dta", clear
collapse (mean) tonmiles* , by(sctg_txt mode2)
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PRailElast = ((tonmilesRail_E-tonmilesBAU)/tonmilesBAU)/.01
gen PBargeElast = ((tonmilesBarge_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PRailElast PBargeElast
xpose, clear varname
order _varname
replace _varname = "`g'" in 1
outsheet using "Tables/`g'ElastTable.csv", comma replace
}

*** 3. Models with truck and rail as main modes
** Temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Alcohol" "Animal Feed" "Milled Grain" "Metallic Ores"
"Other Chemical Products" "Other Prepared Foodstuffs"  { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
* Create some variables for error and predictions
qui: gen PBAUrand = 0
qui: gen PTruckrand_E = 0
qui: gen PRailrand_E = 0
qui: gen PBargerand_E = 0
qui: gen PAirrand_E = 0
qui: save "Data/MixedLogitTemp.dta", replace

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles TempCont) base(6)
qui: matrix B = e(b)

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAUrand = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAUrand = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 6;
#delimit cr
* Find most probably mode
qui: by id: egen tempBAU = max(PBAUrand)
qui: gen maxprobBAU = tempBAU == PBAUrand
* Set mode to rail
qui: replace maxprobBAU = 0 if PBAUrand == . & mode2 == 4
qui: replace maxprobBAU = 1 if PBAUrand == . & mode2 == 6
qui: drop tempBAU

** Truck Elasticity Scenario
* Truck
#delimit ;
qui: replace PTruckrand_E = exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PTruckrand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 6;
#delimit cr
* Find most probably mode
qui: by id: egen tempTruck_E = max(PTruckrand_E)
qui: gen maxprobTruck_E = tempTruck_E == PTruckrand_E
* Set mode to rail
qui: replace maxprobTruck_E = 0 if PTruckrand_E == . & mode2 == 4
qui: replace maxprobTruck_E = 1 if PTruckrand_E == . & mode2 == 6
qui: drop tempTruck_E

** Rail Elasticity Scenario
* Truck
#delimit ;
qui: replace PRailrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PRailrand_E = exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 6;
#delimit cr

* Initialize variables
gen tonmilesBAU = .
gen tonmilesTruck_E = .
gen tonmilesRail_E = .
gen tonmilesBarge_E = .
gen tonmilesAir_E = .

* BAU
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 4
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 6
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 8
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 11
* Truck Elasticity Scenario
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 4
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 6
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 8
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 11
* Rail Elasticity Scenario
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 4
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 6
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 8
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 11
* Barge Elasticity Scenario
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 4
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 6
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 8
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 11
* Air Elasticity Scenario
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 4
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 6
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 8
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 11

qui: collapse (first) sctg_txt (sum) tonmilesBAU tonmilesTruck_E tonmilesRail_E tonmilesBarge_E tonmilesAir_E [iweight = wgt_factor], by(mode2) 

* Unit conversion
qui: replace tonmilesBAU = tonmilesBAU/10^9
qui: replace tonmilesTruck_E = tonmilesTruck_E/10^9
qui: replace tonmilesRail_E = tonmilesRail_E/10^9
qui: replace tonmilesBarge_E = tonmilesBarge_E/10^9
qui: replace tonmilesAir_E = tonmilesAir_E/10^9

qui: save "Tables/`g'Elast.dta", replace

}
* Create elasticity table for each commodity
#delimit ;
foreach g in "Alcohol" "Animal Feed" "Milled Grain" "Metallic Ores"
"Other Chemical Products" "Other Prepared Foodstuffs"  { ;
#delimit cr
use "Tables/`g'Elast.dta", clear
collapse (mean) tonmiles* , by(sctg_txt mode2)
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PRailElast = ((tonmilesRail_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PRailElast
xpose, clear varname
order _varname
replace _varname = "`g'" in 1
outsheet using "Tables/`g'ElastTable.csv", comma replace
}

*** 4. Models with truck and rail as main modes
** Non temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Articles of Base Metal" "Logs and Other Wood in the Rough" "Non-Metallic Mineral Products"
"Paper" "Plastics and Rubber" "Primary Base Metal" "Pulp, Newsprint, Paper, and Paperboard" 
"Sand" "Textiles" "Vehicles" "Wood Products" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
* Create some variables for error and predictions
qui: gen PBAUrand = 0
qui: gen PTruckrand_E = 0
qui: gen PRailrand_E = 0
qui: gen PBargerand_E = 0
qui: gen PAirrand_E = 0
qui: save "Data/MixedLogitTemp.dta", replace

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(6)
qui: matrix B = e(b)

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAUrand = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAUrand = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 6;
#delimit cr
* Find most probably mode
qui: by id: egen tempBAU = max(PBAUrand)
qui: gen maxprobBAU = tempBAU == PBAUrand
* Set mode to rail
qui: replace maxprobBAU = 0 if PBAUrand == . & mode2 == 4
qui: replace maxprobBAU = 1 if PBAUrand == . & mode2 == 6
qui: drop tempBAU

** Truck Elasticity Scenario
* Truck
#delimit ;
qui: replace PTruckrand_E = exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PTruckrand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 6;
#delimit cr
* Find most probably mode
qui: by id: egen tempTruck_E = max(PTruckrand_E)
qui: gen maxprobTruck_E = tempTruck_E == PTruckrand_E
* Set mode to rail
qui: replace maxprobTruck_E = 0 if PTruckrand_E == . & mode2 == 4
qui: replace maxprobTruck_E = 1 if PTruckrand_E == . & mode2 == 6
qui: drop tempTruck_E

** Rail Elasticity Scenario
* Truck
#delimit ;
qui: replace PRailrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PRailrand_E = exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )) if mode2 == 6;
#delimit cr

* Initialize variables
gen tonmilesBAU = .
gen tonmilesTruck_E = .
gen tonmilesRail_E = .
gen tonmilesBarge_E = .
gen tonmilesAir_E = .

* BAU
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 4
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 6
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 8
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 11
* Truck Elasticity Scenario
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 4
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 6
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 8
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 11
* Rail Elasticity Scenario
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 4
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 6
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 8
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 11
* Barge Elasticity Scenario
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 4
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 6
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 8
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 11
* Air Elasticity Scenario
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 4
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 6
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 8
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 11

qui: collapse (first) sctg_txt (sum) tonmilesBAU tonmilesTruck_E tonmilesRail_E tonmilesBarge_E tonmilesAir_E [iweight = wgt_factor], by(mode2) 

* Unit conversion
qui: replace tonmilesBAU = tonmilesBAU/10^9
qui: replace tonmilesTruck_E = tonmilesTruck_E/10^9
qui: replace tonmilesRail_E = tonmilesRail_E/10^9
qui: replace tonmilesBarge_E = tonmilesBarge_E/10^9
qui: replace tonmilesAir_E = tonmilesAir_E/10^9

qui: save "Tables/`g'Elast.dta", replace

}
* Create elasticity table for each commodity
#delimit ;
foreach g in "Articles of Base Metal" "Logs and Other Wood in the Rough" "Non-Metallic Mineral Products"
"Paper" "Plastics and Rubber" "Primary Base Metal" "Pulp, Newsprint, Paper, and Paperboard" 
"Sand" "Textiles" "Vehicles" "Wood Products" { ;
#delimit cr
use "Tables/`g'Elast.dta", clear
collapse (mean) tonmiles* , by(sctg_txt mode2)
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PRailElast = ((tonmilesRail_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PRailElast
xpose, clear varname
order _varname
replace _varname = "`g'" in 1
outsheet using "Tables/`g'ElastTable.csv", comma replace
}

*** 5. Models with Air and Truck
** Temp controlled
* Loop over commodities
#delimit ;
foreach g in  "Animals" "Pharmaceuticals" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 11
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
* Create some variables for error and predictions
qui: gen PBAUrand = 0
qui: gen PTruckrand_E = 0
qui: gen PRailrand_E = 0
qui: gen PBargerand_E = 0
qui: gen PAirrand_E = 0
qui: save "Data/MixedLogitTemp.dta", replace

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles TempCont) base(11)
qui: matrix B = e(b)

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAUrand = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PBAUrand = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempBAU = max(PBAUrand)
qui: gen maxprobBAU = tempBAU == PBAUrand
qui: drop tempBAU

** Truck Elasticity Scenario
* Truck
#delimit ;
qui: replace PTruckrand_E = exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PTruckrand_E = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempTruck_E = max(PTruckrand_E)
qui: gen maxprobTruck_E = tempTruck_E == PTruckrand_E
qui: drop tempTruck_E

** Air Elasticity Scenario
* Truck
#delimit ;
qui: replace PAirrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PAirrand_E = exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 11;
#delimit cr

* Initialize variables
gen tonmilesBAU = .
gen tonmilesTruck_E = .
gen tonmilesRail_E = .
gen tonmilesBarge_E = .
gen tonmilesAir_E = .

* BAU
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 4
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 6
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 8
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 11
* Truck Elasticity Scenario
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 4
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 6
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 8
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 11
* Rail Elasticity Scenario
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 4
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 6
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 8
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 11
* Barge Elasticity Scenario
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 4
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 6
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 8
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 11
* Air Elasticity Scenario
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 4
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 6
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 8
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 11

qui: collapse (first) sctg_txt (sum) tonmilesBAU tonmilesTruck_E tonmilesRail_E tonmilesBarge_E tonmilesAir_E [iweight = wgt_factor], by(mode2) 

* Unit conversion
qui: replace tonmilesBAU = tonmilesBAU/10^9
qui: replace tonmilesTruck_E = tonmilesTruck_E/10^9
qui: replace tonmilesRail_E = tonmilesRail_E/10^9
qui: replace tonmilesBarge_E = tonmilesBarge_E/10^9
qui: replace tonmilesAir_E = tonmilesAir_E/10^9

qui: save "Tables/`g'Elast.dta", replace

}
* Create elasticity table for each commodity
#delimit ;
foreach g in  "Animals" "Pharmaceuticals" { ;
#delimit cr
use "Tables/`g'Elast.dta", clear
collapse (mean) tonmiles* , by(sctg_txt mode2)
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PAirElast = ((tonmilesAir_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PAirElast
xpose, clear varname
order _varname
replace _varname = "`g'" in 1
outsheet using "Tables/`g'ElastTable.csv", comma replace
}

*** 6. Models with air and truck 
** Not temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Precision Instruments" "Printed Products" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 11
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
* Create some variables for error and predictions
qui: gen PBAUrand = 0
qui: gen PTruckrand_E = 0
qui: gen PRailrand_E = 0
qui: gen PBargerand_E = 0
qui: gen PAirrand_E = 0
qui: save "Data/MixedLogitTemp.dta", replace

*** Estimate Mixed Logit and Save Parameter estimates 
qui: asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(11)
qui: matrix B = e(b)

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAUrand = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PBAUrand = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempBAU = max(PBAUrand)
qui: gen maxprobBAU = tempBAU == PBAUrand
qui: drop tempBAU

** Truck Elastcitity Scenario
* Truck
#delimit ;
qui: replace PTruckrand_E = exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PTruckrand_E = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempTruck_E = max(PTruckrand_E)
qui: gen maxprobTruck_E = tempTruck_E == PTruckrand_E
qui: drop tempTruck_E

** Air Elastcitity Scenario
* Truck
#delimit ;
qui: replace PAirrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 4;
#delimit cr
* Air
#delimit ;
qui: replace PAirrand_E = exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] )) if mode2 == 11;
#delimit cr

* Initialize variables
gen tonmilesBAU = .
gen tonmilesTruck_E = .
gen tonmilesRail_E = .
gen tonmilesBarge_E = .
gen tonmilesAir_E = .

* BAU
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 4
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 6
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 8
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 11
* Truck Elasticity Scenario
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 4
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 6
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 8
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 11
* Rail Elasticity Scenario
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 4
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 6
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 8
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 11
* Barge Elasticity Scenario
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 4
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 6
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 8
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 11
* Air Elasticity Scenario
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 4
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 6
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 8
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 11

qui: collapse (first) sctg_txt (sum) tonmilesBAU tonmilesTruck_E tonmilesRail_E tonmilesBarge_E tonmilesAir_E [iweight = wgt_factor], by(mode2) 

* Unit conversion
qui: replace tonmilesBAU = tonmilesBAU/10^9
qui: replace tonmilesTruck_E = tonmilesTruck_E/10^9
qui: replace tonmilesRail_E = tonmilesRail_E/10^9
qui: replace tonmilesBarge_E = tonmilesBarge_E/10^9
qui: replace tonmilesAir_E = tonmilesAir_E/10^9

qui: save "Tables/`g'Elast.dta", replace

}
* Create elasticity table for each commodity
#delimit ;
foreach g in "Precision Instruments" "Printed Products" { ;
#delimit cr
use "Tables/`g'Elast.dta", clear
collapse (mean) tonmiles* , by(sctg_txt mode2)
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PAirElast = ((tonmilesAir_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PAirElast
xpose, clear varname
order _varname
replace _varname = "`g'" in 1
outsheet using "Tables/`g'ElastTable.csv", comma replace
}

*** 7. Models with truck, rail and air as main modes
** Temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Mixed Freight" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 | mode == 11
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
* Create some variables for error and predictions
qui: gen PBAUrand = 0
qui: gen PTruckrand_E = 0
qui: gen PRailrand_E = 0
qui: gen PBargerand_E = 0
qui: gen PAirrand_E = 0
qui: save "Data/MixedLogitTemp.dta", replace

*** Estimate Mixed Logit and Save Parameter estimates 
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles TempCont) base(6)
qui: matrix B = e(b)

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAUrand = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PBAUrand = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PBAUrand = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempBAU = max(PBAUrand)
qui: gen maxprobBAU = tempBAU == PBAUrand
* Set mode to rail
qui: replace maxprobBAU = 0 if PBAUrand == . & mode2 == 4
qui: replace maxprobBAU = 1 if PBAUrand == . & mode2 == 6
qui: replace maxprobBAU = 0 if PBAUrand == . & mode2 == 11
qui: drop tempBAU

** Truck Elasticity Scenario
* Truck  
#delimit ;
qui: replace PTruckrand_E = exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PTruckrand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PTruckrand_E = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 11;
#delimit cr

** Rail Elasticity Scenario
* Truck  
#delimit ;
qui: replace PRailrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PRailrand_E = exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PRailrand_E = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempRail_E = max(PRailrand_E)
qui: gen maxprobRail_E = tempRail_E == PRailrand_E
* Set mode to rail
qui: replace maxprobRail_E = 0 if PRailrand_E == . & mode2 == 4
qui: replace maxprobRail_E = 1 if PRailrand_E == . & mode2 == 6
qui: replace maxprobRail_E = 0 if PRailrand_E == . & mode2 == 11
qui: drop tempRail_E

** Air Elasticity Scenario
* Truck  
#delimit ;
qui: replace PAirrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 4;
#delimit cr
* Rail
#delimit ;
qui: replace PAirrand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PAirrand_E = exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + TempCont*B[1,3] + B[1,4] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,5] + TempCont*B[1,6] + B[1,7] ) ) if mode2 == 11;
#delimit cr

* Initialize variables
gen tonmilesBAU = .
gen tonmilesTruck_E = .
gen tonmilesRail_E = .
gen tonmilesBarge_E = .
gen tonmilesAir_E = .

* BAU
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 4
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 6
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 8
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 11
* Truck Elasticity Scenario
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 4
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 6
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 8
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 11
* Rail Elasticity Scenario
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 4
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 6
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 8
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 11
* Barge Elasticity Scenario
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 4
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 6
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 8
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 11
* Air Elasticity Scenario
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 4
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 6
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 8
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 11

qui: collapse (first) sctg_txt (sum) tonmilesBAU tonmilesTruck_E tonmilesRail_E tonmilesBarge_E tonmilesAir_E [iweight = wgt_factor], by(mode2) 

* Unit conversion
qui: replace tonmilesBAU = tonmilesBAU/10^9
qui: replace tonmilesTruck_E = tonmilesTruck_E/10^9
qui: replace tonmilesRail_E = tonmilesRail_E/10^9
qui: replace tonmilesBarge_E = tonmilesBarge_E/10^9
qui: replace tonmilesAir_E = tonmilesAir_E/10^9

qui: save "Tables/`g'Elast.dta", replace

}
* Create elasticity table for each commodity
#delimit ;
foreach g in "Mixed Freight" { ;
#delimit cr
use "Tables/`g'Elast.dta", clear
collapse (mean) tonmiles* , by(sctg_txt mode2)
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PRailElast = ((tonmilesRail_E-tonmilesBAU)/tonmilesBAU)/.01
gen PAirElast = ((tonmilesAir_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PRailElast PAirElast
xpose, clear varname
order _varname
replace _varname = "`g'" in 1
outsheet using "Tables/`g'ElastTable.csv", comma replace
}

*** 8. Models with truck, rail and air as main modes
** Not temperature controlled
* Loop over commodities
#delimit ;
foreach g in "Machinery" "Miscellaneous Manufactured Products" { ;
#delimit cr

** Prepare data for mixed logit implementation
qui: use "Data/CFS2012.dta", clear
qui: if $useboth == 1{
append using "Data/CFS2017.dta",
}
qui: keep if sctg_txt == "`g'"
display "`g'"
qui: keep if mode == 4 | mode == 6 | mode == 11
** Reshape data to long form
* Make id variable
qui: gen id = _n
* Create a variable for reshape
qui: gen rt4 = routeid
qui: gen rt6 = routeid
qui: gen rt11 = routeid
* Reshape to multiple observations per case
qui: reshape long rt, i(id) j(mode2)
* Create interaction variables for mixed logit implementation
qui: gen choice = mode==mode2
qui: gen ValMiles = value*miles/10^6
* Create variable implementing differences in fuel intensity across modes
qui: gen PdTM = .
qui: replace PdTM = $gptmTruck*Pfuel*tonmiles/10^6 if mode2 == 4
qui: replace PdTM = $gptmRail*Pfuel*tonmiles/10^6 if mode2 == 6
qui: replace PdTM = $gptmAir*Pfuel*tonmiles/10^6 if mode2 == 11
* Create some variables for error and predictions
qui: gen PBAUrand = 0
qui: gen PTruckrand_E = 0
qui: gen PRailrand_E = 0
qui: gen PBargerand_E = 0
qui: gen PAirrand_E = 0
qui: save "Data/MixedLogitTemp.dta", replace

*** Estimate Mixed Logit and Save Parameter estimates 
asclogit choice PdTM [iweight = wgt_factor], case(id) alternatives(mode2) casevars(ValMiles) base(6)
qui: matrix B = e(b)

** Estimate BAU mode probabilities
* Truck
#delimit ;
qui: replace PBAUrand = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] ) ) if mode2 == 4;
#delimiti cr
* Rail
#delimit ;
qui: replace PBAUrand = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] ) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PBAUrand = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] )  /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] )  ) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempBAU = max(PBAUrand)
qui: gen maxprobBAU = tempBAU == PBAUrand
* Set mode to rail
qui: replace maxprobBAU = 0 if PBAUrand == . & mode2 == 4
qui: replace maxprobBAU = 1 if PBAUrand == . & mode2 == 6
qui: replace maxprobBAU = 0 if PBAUrand == . & mode2 == 11
qui: drop tempBAU

** Truck Elasticity Scenario  
* Truck
#delimit ;
qui: replace PTruckrand_E = exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] ) ) if mode2 == 4;
#delimiti cr
* Rail
#delimit ;
qui: replace PTruckrand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] ) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PTruckrand_E = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] )  /
(exp((1.01)*$gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] )  ) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempTruck_E = max(PTruckrand_E)
qui: gen maxprobTruck_E = tempTruck_E == PTruckrand_E
* Set mode to rail
qui: replace maxprobTruck_E = 0 if PTruckrand_E == . & mode2 == 4
qui: replace maxprobTruck_E = 1 if PTruckrand_E == . & mode2 == 6
qui: replace maxprobTruck_E = 0 if PTruckrand_E == . & mode2 == 11
qui: drop tempTruck_E

** Rail Elasticity Scenario  
* Truck
#delimit ;
qui: replace PRailrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] ) ) if mode2 == 4;
#delimiti cr
* Rail
#delimit ;
qui: replace PRailrand_E = exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] ) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PRailrand_E = exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] )  /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp((1.01)*$gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp($gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] )  ) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempRail_E = max(PRailrand_E)
qui: gen maxprobRail_E = tempRail_E == PRailrand_E
   
* Set mode to rail
qui: replace maxprobRail_E = 0 if PRailrand_E == . & mode2 == 4
qui: replace maxprobRail_E = 1 if PRailrand_E == . & mode2 == 6
qui: replace maxprobRail_E = 0 if PRailrand_E == . & mode2 == 11
qui: drop tempRail_E

** Air Elasticity Scenario  
* Truck
#delimit ;
qui: replace PAirrand_E = exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )/
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] ) ) if mode2 == 4;
#delimiti cr
* Rail
#delimit ;
qui: replace PAirrand_E = exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  ) /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] ) ) if mode2 == 6;
#delimit cr
* Air
#delimit ;
qui: replace PAirrand_E = exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] )  /
(exp($gptmTruck*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,2] + B[1,3] )
+ exp($gptmRail*Pfuel*(tonmiles/10^6)*B[1,1]  )
+ exp((1.01)*$gptmAir*Pfuel*(tonmiles/10^6)*B[1,1] + ValMiles*B[1,4] + B[1,5] )  ) if mode2 == 11;
#delimit cr
* Find most probably mode
qui: by id: egen tempAir_E = max(PAirrand_E)
qui: gen maxprobAir_E = tempAir_E == PAirrand_E
* Set mode to rail
qui: replace maxprobAir_E = 0 if PAirrand_E == . & mode2 == 4
qui: replace maxprobAir_E = 1 if PAirrand_E == . & mode2 == 6
qui: replace maxprobAir_E = 0 if PAirrand_E == . & mode2 == 11
qui: drop tempAir_E

* Initialize variables
gen tonmilesBAU = .
gen tonmilesTruck_E = .
gen tonmilesRail_E = .
gen tonmilesBarge_E = .
gen tonmilesAir_E = .

* BAU
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 4
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 6
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 8
replace tonmilesBAU = PBAUrand*tonmiles if mode2 == 11
* Truck Elasticity Scenario
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 4
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 6
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 8
replace tonmilesTruck_E = PTruckrand_E*tonmiles if mode2 == 11
* Rail Elasticity Scenario
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 4
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 6
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 8
replace tonmilesRail_E = PRailrand_E*tonmiles if mode2 == 11
* Barge Elasticity Scenario
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 4
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 6
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 8
replace tonmilesBarge_E = PBargerand_E*tonmiles if mode2 == 11
* Air Elasticity Scenario
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 4
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 6
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 8
replace tonmilesAir_E = PAirrand_E*tonmiles if mode2 == 11

qui: collapse (first) sctg_txt (sum) tonmilesBAU tonmilesTruck_E tonmilesRail_E tonmilesBarge_E tonmilesAir_E [iweight = wgt_factor], by(mode2) 

* Unit conversion
qui: replace tonmilesBAU = tonmilesBAU/10^9
qui: replace tonmilesTruck_E = tonmilesTruck_E/10^9
qui: replace tonmilesRail_E = tonmilesRail_E/10^9
qui: replace tonmilesBarge_E = tonmilesBarge_E/10^9
qui: replace tonmilesAir_E = tonmilesAir_E/10^9

qui: save "Tables/`g'Elast.dta", replace

}
* Create elasticity table for each commodity
#delimit ;
foreach g in "Machinery" "Miscellaneous Manufactured Products" { ;
#delimit cr
use "Tables/`g'Elast.dta", clear
collapse (mean) tonmiles* , by(sctg_txt mode2)
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PRailElast = ((tonmilesRail_E-tonmilesBAU)/tonmilesBAU)/.01
gen PAirElast = ((tonmilesAir_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PRailElast PAirElast
xpose, clear varname
order _varname
replace _varname = "`g'" in 1
outsheet using "Tables/`g'ElastTable.csv", comma replace
}

* Horrizonatal bar chart for truck-rail cross-elasticity
* Calculate mean elasticities
#delimit ;
*foreach g in "Articles of Base Metal" { ;
foreach g in "Agricultural Products" "Basic Chemicals"  "Fertilizers" "Other Coal and Petroleum" 
"Coal" "Grain" "Gravel" "Transportation Equipment, not elsewhere classified" "Waste and Scrap" 
"Alcohol" "Animal Feed" "Milled Grain" "Metallic Ores"
"Other Chemical Products" "Other Prepared Foodstuffs"  
"Articles of Base Metal" "Logs and Other Wood in the Rough" "Non-Metallic Mineral Products"
"Paper" "Plastics and Rubber" "Primary Base Metal" "Pulp, Newsprint, Paper, and Paperboard" 
"Sand" "Textiles" "Vehicles" "Wood Products" "Mixed Freight" 
"Machinery" "Miscellaneous Manufactured Products"{ ;
#delimit cr
* Open data file
use "Tables/`g'Elast.dta", clear
* Create elasticity table for each commodity
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PRailElast = ((tonmilesRail_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PRailElast
save "Tables/`g'Elast_summary.dta", replace
}

* Combine files
drop _all
#delimit ;
foreach g in "Agricultural Products" "Basic Chemicals"  "Fertilizers" "Other Coal and Petroleum" 
"Coal" "Grain" "Gravel" "Transportation Equipment, not elsewhere classified" "Waste and Scrap" 
"Alcohol" "Animal Feed" "Milled Grain" "Metallic Ores"
"Other Chemical Products" "Other Prepared Foodstuffs"  
"Articles of Base Metal" "Logs and Other Wood in the Rough" "Non-Metallic Mineral Products"
"Paper" "Plastics and Rubber" "Primary Base Metal" "Pulp, Newsprint, Paper, and Paperboard" 
"Sand" "Textiles" "Vehicles" "Wood Products" "Mixed Freight" 
"Machinery" "Miscellaneous Manufactured Products"{ ;
# delimit cr
append using "Tables/`g'Elast_summary.dta",
}
keep sctg_txt mode2 PTruckElast
gsort - mode2 + PTruckElast 
* Clean up transportation label
replace sctg_txt = "Transportation Equipment" if sctg_txt == "Transportation Equipment, not elsewhere classified"
* Keep only cross-price elasticities
keep if mode2 == 6
gsort - PTruckElast
gen order = _n
#delimit;
graph hbar (mean) PTruckElast, over(sctg_txt,sort(order))
title("Rail to Truck Cross-Price Elasticity") 
ytitle("Elasticity")
bar(1,color(purple%75));
#delimit cr
graph export "Figures/ElastRailtoTruck.pdf", replace

* Horrizonatal bar chart for truck-Air cross-elasticity
* Calculate mean elasticities
#delimit ;
foreach g in "Animals" "Pharmaceuticals" 
"Precision Instruments" "Printed Products" { ;
#delimit cr

* Open data file
use "Tables/`g'Elast.dta", clear
* Create elasticity table for each commodity
gen PTruckElast = ((tonmilesTruck_E-tonmilesBAU)/tonmilesBAU)/.01
gen PAirElast = ((tonmilesAir_E-tonmilesBAU)/tonmilesBAU)/.01
keep sctg_txt mode2 PTruckElast PAirElast
save "Tables/`g'Elast_summary.dta", replace
}

* Combine files
drop _all
#delimit ;
foreach g in "Animals" "Pharmaceuticals" 
"Precision Instruments" "Printed Products"{ ;
# delimit cr
append using "Tables/`g'Elast_summary.dta",
}
keep sctg_txt mode2 PTruckElast
gsort - mode2 + PTruckElast 
* Keep only cross-price elasticities
keep if mode2 == 11
gsort - PTruckElast
gen order = _n
#delimit;
graph hbar (mean) PTruckElast, over(sctg_txt,sort(order))
title("Air to Truck Cross-Price Elasticity") 
ytitle("Elasticity")
bar(1,color(gold%50));
#delimit cr
graph export "Figures/ElastAirtoTruck.pdf", replace
