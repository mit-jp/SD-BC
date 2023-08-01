function mdpcxy(args)

'reinit'
var = subwrd(args,1)

* these models are used to generate the dataset for pnnl
models = 'ACCESS-ESM1-5 AWI-ESM-1-1-LR BCC-CSM2-MR CanESM5 CMCC-ESM2 CNRM-ESM2-1 EC-Earth3-Veg FGOALS-g3 FIO-ESM-2-0 GISS-E2-2-G HadGEM3-GC31-MM INM-CM5-0 IPSL-CM6A-LR MIROC-ES2L MPI-ESM1-2-HR MRI-ESM2-0 SAM0-UNICON UKESM1-0-LL' 

ll =1 
while(ll <= 18)
  mdl = subwrd(models,ll)
  say mdl
  'open 'mdl'_1pctCO2_'var'_0360x0720.ctl'
  ll = ll + 1
endwhile

'q file 1'
grdinfo = sublin(result,5)
xmax = subwrd(grdinfo,3)
ymax = subwrd(grdinfo,6)
'set x 1 'xmax
'set y 1 'ymax
'set z 1'
'set gxout fwrite'
'set fwrite -sq Multimodel_1pctCO2_'var'_0360x0720.gr'

tt = 1
while(tt <= 12)
  'set t 'tt
  'd ('var'.1+'var'.2+'var'.3+'var'.4+'var'.5+'var'.6+'var'.7+'var'.8+'var'.9+'var'.10+'var'.11+'var'.12+'var'.13+'var'.14+'var'.15+'var'.16+'var'.17+'var'.18)/18'
  tt = tt + 1
endwhile 

'disable fwrite'

ll = 18
while(ll >= 1)
  'close 'll
  ll = ll - 1
endwhile

'quit'


