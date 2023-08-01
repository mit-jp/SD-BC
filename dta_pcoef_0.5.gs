function dcoef(args)
* this script first interpolate the variables on 2x2.5 grid, then calculates
* the difference of HFD coefficients difference between 71th~80th year and 
* 1st 10-year of 1pctCO2 run
* usage example dta_pcoef_2x2.5.gs tas(variable1) pr(variable 2) dpcxy
* does not include the model HadGEM3-GC31-MM

'reinit'
var1 = subwrd(args,1)
var2 = subwrd(args,2)
var3 = subwrd(args,3)

models = 'ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR AWI-ESM-1-1-LR BCC-CSM2-MR BCC-ESM1 CAMS-CSM1-0 CanESM5 CanESM5-CanOE CAS-ESM2-0 CESM2 CESM2-WACCM CESM2-WACCM-FV2 CIESM CMCC-CM2-HR4 CMCC-CM2-SR5 CMCC-ESM2 CNRM-CM6-1 CNRM-CM6-1-HR CNRM-ESM2-1 E3SM-1-0 EC-Earth3 EC-Earth3-AerChem EC-Earth3-CC EC-Earth3-Veg FGOALS-f3-L FGOALS-g3 FIO-ESM-2-0 GISS-E2-1-G GISS-E2-1-H GISS-E2-2-G HadGEM3-GC31-LL IITM-ESM INM-CM4-8 INM-CM5-0 IPSL-CM5A2-INCA IPSL-CM6A-LR KACE-1-0-G MCM-UA-1-0 MIROC6 MIROC-ES2L MPI-ESM-1-2-HAM MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorCPM1 NorESM2-LM NorESM2-MM SAM0-UNICON TaiESM1 UKESM1-0-LL'

'open ../../dumgrid_0.5.ctl'
'q file 1'
grdinfo = sublin(result,5)
xmax = subwrd(grdinfo,3)
ymax = subwrd(grdinfo,6)
'set x 1 'xmax
'set y 1 'ymax
'set z 1'
'set gxout fwrite'

ll = 1
while(ll <= 52)
  mdl = subwrd(models,ll)
  'open ../../../'var1'/'var1'_Amon_'mdl'.ctl'
  'open ../../'var2'_Amon_'mdl'.ctl'
  'set fwrite -sq 'mdl'_1pctCO2_'var3'_0360x0720.gr'
  'q file 2'
  grdinfo = sublin(result,5)
  tmx = subwrd(grdinfo,12)
  say mdl' 'tmx
  'set t 1 'tmx
* interpolate the variables to 2x2.5 grid
* v1 - tas at 2x2.5 of the entire time series
* v2 - pr at 2x2.5 of the entire time series (unit: mm/day)
*zv2 - zonal pr at 2x1 grid of the entire time series
  'v1 = lterp('var1'.2,grid.1(t=1),aave)'
  'v2 = lterp('var2'.3*86400,grid.1(t=1),aave)'
  'set x 1'
  'zv2 = ave(v2,x=1,x='xmax')'
  mm = 1
  while(mm <= 12)
* use the first 10 years of 1pctCO2 run for averaging period
* use the 71th ~ 80 years of 1pctCO2 run for averaging period 
* as 70 years is the point that co2 gets doubled.
* cznlv2 - climatology of a month for zonal pr over the first 10-year
* cgrdv2 - climatology of a month for gridded pr at 2x2.5 over the first 10-year
* cgrdv1 - climatology of a month for gridded tas at 2x2.5 over the first 10-year
* ccxy - grid to zonal coefficient for pr at 2x2.5 over the first 10-year
* znlv2 - climatology of a month for zonal pr over the 71st-80th year
* grdv2 - climatology of a month for gridded pr at 2x2.5 over the 71st-80th year
* grdv1 - climatology of a month for gridded tas at 2x2.5 over the 71st-80th year
* cxy - grid to zonal coefficient for pr at 2x2.5 over the 71st-80th year
* denominator of dcxy is the global mean tas difference between 71st-80th and the first 10 years
    t2 = 840 + mm
    tmx2 = 840 + 120
    'set t 'mm
    'cznlv2= ave(zv2,t+0,t=120,12)'
    'set x 1 'xmax
    'cgrdv2= ave(v2,t+0,t=120,12)'
    'cgrdv1= ave(v1,t+0,t=120,12)'
    'ccxy = cgrdv2/cznlv2'
    'set t 't2
    'set x 1'
    'znlv2 = ave(zv2,t+0,t='tmx2',12)'
    'set x 1 'xmax
    'grdv2 = ave(v2,t+0,t='tmx2',12)'
    'grdv1 = ave(v1,t+0,t='tmx2',12)'
    'cxy = grdv2/znlv2'
    'dcxy =(cxy-ccxy(t='mm'))/aave(grdv1-cgrdv1(t='mm'),x=1,x='xmax',y=1,y='ymax')'
    'd dcxy'
    mm = mm + 1
  endwhile
  'disable fwrite'
  'close 3'
  'close 2'
  ll = ll + 1
endwhile

'quit'


