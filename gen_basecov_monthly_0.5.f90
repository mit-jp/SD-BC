PROGRAM genmon 

!** this program applies HFD coefficients(climatology cxy from MERRA2 plus
!** changes of cxy from CMIP6) directly to MESM zonal output to generate
!** monhtly meteorological forcing for Tufts from 2011 to 2150

USE NETCDF
IMPLICIT NONE

INTEGER, PARAMETER :: nx = 720, ny =360, nzy = 46
INTEGER, PARAMETER :: nyr =  80, yr0 = 2020, nm = 12, nyr1 = 100
INTEGER, PARAMETER :: nmdl=18, nens = 50
INTEGER, PARAMETER :: nyrt =   89!** the number of years for doing running mean of global temperature
INTEGER, PARAMETER :: nyrc = 10  !** the number of years for calculating moving average of gt 
INTEGER, PARAMETER :: nyrz = 20  !** the number of years for calculating moving average of zonal values 
REAL,    PARAMETER :: missing = -9.99e+08, rndns = 5.0
CHARACTER(LEN=15), PARAMETER :: sce = 'BASECOV       ' 
!CHARACTER(LEN=15), PARAMETER :: sce = 'PFCOV          ' 
!CHARACTER(LEN=15), PARAMETER :: sce = 'PARIS_2C       ' 
!CHARACTER(LEN=15), PARAMETER :: sce = 'PARIS_1p5C_2030' 

CHARACTER(LEN=3)   :: ce3
CHARACTER(LEN= 4)  :: che3,cfe3,dumtbl(33),dumc4m1,dumc4m2
CHARACTER(LEN=15)  :: cmodel(nmdl)
CHARACTER(LEN=250) :: fin0,fin1,fin2,fin3,fin4,fin5,fin6,fin7,fin8,fin9,fdir,cdum
CHARACTER(LEN=250) :: fout1,fout2,fout3,fout4,fout5,fout6,fout7,fout8,fout9
INTEGER :: dumd1,dumd2,dumyr1,dumyr2,inqtab(33),j1qt(33),inqmap(57)
INTEGER :: iret,enum(2),tstep(nyr*nm),istart(3),icount(3)
INTEGER :: oid1,oid2,oid3,oid4,oid5,oid6,oid7,oid8,oid9
INTEGER :: londim,latdim,time_dim,one_dim(1),three_dims(3)
INTEGER :: lonid,latid,timeid,prepid,tairid,tmaxid,tminid,windid,swdownid,lwdownid,rhid,shid
INTEGER :: trprepid,trtairid,trtmaxid,trtminid,trwindid,trswdownid,trlwdownid,trshid,trrhid
INTEGER :: i,iy,im,icnt,icnt0,iix,iiy,irx,iry,jindex,imdl,ie
REAL    :: gbudg(49,89,4),dumexp,dumtau,qmaps(49,57),qtable(49,12,33)
REAl    :: gt(nyrt,nm),dgt(nyr,nm),sumt,baset,tmp1,tmp2,rnum
REAL    :: summp,summt,summtx,summtn,summsw,summlw,summws,summsh,summrh
REAL    :: lon(nx),lat(ny),tmp(nx,ny)
REAL    :: pcxy(nx,ny,nm),tascxy(nx,ny,nm),tmxcxy(nx,ny,nm),tmncxy(nx,ny,nm),swcxy(nx,ny,nm)
REAL    :: lwcxy(nx,ny,nm),wscxy(nx,ny,nm),shcxy(nx,ny,nm),hurcxy(nx,ny,nm),psrf(nx,ny,nm)
REAL    :: dpcxy(nx,ny,nm,nmdl),dtascxy(nx,ny,nm,nmdl),dtmxcxy(nx,ny,nm,nmdl),dtmncxy(nx,ny,nm,nmdl)
REAL    :: dswcxy(nx,ny,nm,nmdl),dlwcxy(nx,ny,nm,nmdl),dwindcxy(nx,ny,nm,nmdl),dshcxy(nx,ny,nm,nmdl),dhurcxy(nx,ny,nm,nmdl)
REAL    :: tpxy(nx,ny),ttxy(nx,ny),ttmxxy(nx,ny),ttmnxy(nx,ny),tswxy(nx,ny),tlwxy(nx,ny),twsxy(nx,ny),tshxy(nx,ny),thurxy(nx,ny)
!** this is the total fields
REAL    :: prep(nx,ny),tair(nx,ny),tmax(nx,ny),tmin(nx,ny),sw(nx,ny),lw(nx,ny),wind(nx,ny),rh(nx,ny),huss(nx,ny)
!** this is the delta fields (IGSM trend)
!REAL    :: dlt_prep(nx,ny),dlt_tair(nx,ny),dlt_tmax(nx,ny),dlt_tmin(nx,ny),dlt_sw(nx,ny),dlt_lw(nx,ny),dlt_wind(nx,ny),dlt_huss(nx,ny)
REAL    :: dlt_prep(nx,ny,nyr,nm),dlt_tair(nx,ny,nyr,nm),dlt_tmax(nx,ny,nyr,nm),dlt_tmin(nx,ny,nyr,nm)
REAL    :: dlt_sw(nx,ny,nyr,nm),dlt_lw(nx,ny,nyr,nm),dlt_wind(nx,ny,nyr,nm),dlt_huss(nx,ny,nyr,nm),dlt_rh(nx,ny,nyr,nm)
!** future IGSM zonal values
REAL    :: zp(nzy,nyr,nm),zt(nzy,nyr,nm),ztmx(nzy,nyr,nm),ztmn(nzy,nyr,nm),zsw(nzy,nyr,nm)
REAL    :: zlw(nzy,nyr,nm),zws(nzy,nyr,nm),zsh(nzy,nyr,nm),zrh(nzy,nyr,nm)
REAL    :: zp0(nzy,nyr1,nm),zt0(nzy,nyr1,nm),ztmx0(nzy,nyr1,nm),ztmn0(nzy,nyr1,nm),zsw0(nzy,nyr1,nm)
REAL    :: zlw0(nzy,nyr1,nm),zws0(nzy,nyr1,nm),zsh0(nzy,nyr1,nm),zrh0(nzy,nyr1,nm)
!** GSWP3 1931-2010 climatology zonal values
!REAL    :: ozp(ny,nm),ozt(ny,nm),oztmx(ny,nm),oztmn(ny,nm),ozsw(ny,nm),ozlw(ny,nm),ozws(ny,nm),ozsh(ny,nm)
!** GSWP3 HFD-derived climatology field
!REAL    :: op(nx,ny,nm),ot(nx,ny,nm),otmx(nx,ny,nm),otmn(nx,ny,nm),osw(nx,ny,nm),olw(nx,ny,nm),ows(nx,ny,nm),osh(nx,ny,nm)
!** IGSM 2001-2020 climatology zonal values
REAL    :: mzp(nzy,nm),mzt(nzy,nm),mztmx(nzy,nm),mztmn(nzy,nm),mzsw(nzy,nm),mzlw(nzy,nm),mzws(nzy,nm),mzsh(nzy,nm),mzrh(nzy,nm)
!** GSWP3-based 80-yr baseline climate for delta method
REAL    :: p0(nx,ny,nyr*nm),tas0(nx,ny,nyr*nm),tmx0(nx,ny,nyr*nm),tmn0(nx,ny,nyr*nm)
REAL    :: sw0(nx,ny,nyr*nm),lw0(nx,ny,nyr*nm),ws0(nx,ny,nyr*nm),sh0(nx,ny,nyr*nm),hur0(nx,ny,nyr*nm),tmp0(nx,ny)
logical :: flag

!** historical and BASECOV
enum = (/1500,1600/)
!** historical and PFCOV
!enum = (/1500,1650/)
!** historical and PARIS_2C
!enum = (/1500,1700/)
!** historical and PARIS_1p5C
!enum = (/1500,1750/)

cmodel = (/'ACCESS-ESM1-5  ','AWI-ESM-1-1-LR ','BCC-CSM2-MR    ','CanESM5        ','CMCC-ESM2      ','CNRM-ESM2-1    ',&
           'EC-Earth3-Veg  ','FGOALS-g3      ','FIO-ESM-2-0    ','GISS-E2-2-G    ','HadGEM3-GC31-MM','INM-CM5-0      ',&
           'IPSL-CM6A-LR   ','MIROC-ES2L     ','MPI-ESM1-2-HR  ','MRI-ESM2-0     ','SAM0-UNICON    ','UKESM1-0-LL    '/)

flag = .true.

DO iix = 1, nx
  lon(iix) = -179.75+(iix-1)*0.5
ENDDO
DO iiy = 1, ny
  lat(iiy) = -89.75+ (iiy-1)*0.5
ENDDO

DO iy = 1, nyr*nm
  tstep(iy) = iy-1
ENDDO

!** read GSWP3-based baseline climate from 1931-2010
!** the orginal file (baseclim_0.5_1931_2010_v2.dat) starts from 0.25 degree instead of -179.75 as in 2x2.5 degree
fin0 = '/net/fs05/d1/xgao/cesm2/Analyses/GSWP3/forcing/Monthly/script/baseclim_0.5_1931_2010_v2_west.dat'
open(11,file=fin0,status='old',form='unformatted')
DO i = 1, nyr*nm
  read(11) p0(:,:,i)
  read(11) sw0(:,:,i)
  read(11) lw0(:,:,i)
  read(11) sh0(:,:,i)
  read(11) tas0(:,:,i)
  read(11) ws0(:,:,i)
  read(11) tmn0(:,:,i)
  read(11) tmx0(:,:,i)
  read(11) tmp0
  read(11) hur0(:,:,i)
ENDDO
close(11)

!=================================================================================
!
!  Read climatology of zonal to grid coefficients from GSWP3
!==============================================================================

open(11,file='/home/xgao/DOE/data/GSWP3/GSWP3_9vars_Cxy_9110_0.5.gr',status='old',form='unformatted')
DO i = 1, nm
  read(11) pcxy(:,:,i)
  read(11) swcxy(:,:,i)
  read(11) lwcxy(:,:,i)
  read(11) shcxy(:,:,i)
  read(11) tascxy(:,:,i)
  read(11) wscxy(:,:,i)
  read(11) tmncxy(:,:,i)
  read(11) tmxcxy(:,:,i)
  read(11) hurcxy(:,:,i)
ENDDO
close(11)

DO imdl = 1, nmdl
!=================================================================================
!
!  Read trends of forcing variables coefficients from particular
!  AR4 GCM (dcxy/dT)
!
!=================================================================================
!** The default maximum width of a line is 132,otherwise, gives error

fin1 = '/home/xgao/HFD/CMIP6/1pctCO2/pr/zonal/res0.5/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_dpcxy_0360x0720.gr'
fin2 = '/home/xgao/HFD/CMIP6/1pctCO2/tas/zonal/res0.5/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_dtcxy_0360x0720.gr' 
fin3 = '/home/xgao/HFD/CMIP6/1pctCO2/rsds/zonal/res0.5/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_dswcxy_0360x0720.gr'
fin4 = '/home/xgao/HFD/CMIP6/1pctCO2/huss/zonal/res0.5/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_dshcxy_0360x0720.gr'
fin0 = '/home/xgao/HFD/CMIP6/1pctCO2/sfcWind/global/res0.5/'
fin5 = fin0(1:len_trim(fin0))//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_dwindcxy_0360x0720_global.gr'
fin6 = '/home/xgao/HFD/CMIP6/1pctCO2/tasmax/zonal/res0.5/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_dtmaxcxy_0360x0720.gr'
fin7 = '/home/xgao/HFD/CMIP6/1pctCO2/tasmin/zonal/res0.5/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_dtmincxy_0360x0720.gr'
fin8 = '/home/xgao/HFD/CMIP6/1pctCO2/rlds/zonal/res0.5/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_dlwcxy_0360x0720.gr'
fin9 = '/home/xgao/HFD/CMIP6/1pctCO2/hurs/zonal/res0.5/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))//'_1pctCO2_drhcxy_0360x0720.gr'
open(11,file=fin1,status='old',form='unformatted')
open(12,file=fin2,status='old',form='unformatted')
open(13,file=fin3,status='old',form='unformatted')
open(14,file=fin4,status='old',form='unformatted')
open(15,file=fin5,status='old',form='unformatted')
open(16,file=fin6,status='old',form='unformatted')
open(17,file=fin7,status='old',form='unformatted')
open(18,file=fin8,status='old',form='unformatted')
open(19,file=fin9,status='old',form='unformatted')

DO i = 1, nm
  read(11) dpcxy(:,:,i,imdl)
  read(12) dtascxy(:,:,i,imdl)
  read(13) dswcxy(:,:,i,imdl)
  read(14) dshcxy(:,:,i,imdl)
  read(15) dwindcxy(:,:,i,imdl)
  read(16) dtmxcxy(:,:,i,imdl)
  read(17) dtmncxy(:,:,i,imdl)
  read(18) dlwcxy(:,:,i,imdl)
  read(19) dhurcxy(:,:,i,imdl)
ENDDO
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
close(17)
close(18)
close(19)

!*** deal with missing data to accomodate the 4 combinations of cases
!*** cxy is missing, no matter dcxy is missing or not, the result is 1.0*IGSM zonal value
!*** Cxy is not mssing, if dcxy is missing, the result is Cxy*IGSM zonal value;
!*** dcxy is not missing, use the full equation - (cxy+dt*(dcxy/dt))*IGSM zonal value
where(dpcxy(:,:,:,imdl) == missing .or. pcxy == missing) dpcxy(:,:,:,imdl) = 0.
where(dtascxy(:,:,:,imdl) == missing .or. tascxy == missing) dtascxy(:,:,:,imdl) = 0.
where(dtmxcxy(:,:,:,imdl) == missing .or. tmxcxy == missing) dtmxcxy(:,:,:,imdl) = 0.
where(dtmncxy(:,:,:,imdl) == missing .or. tmncxy == missing) dtmncxy(:,:,:,imdl) = 0.
where(dswcxy(:,:,:,imdl) == missing .or. swcxy == missing) dswcxy(:,:,:,imdl) = 0.
where(dwindcxy(:,:,:,imdl) == missing .or. wscxy == missing) dwindcxy(:,:,:,imdl) = 0.
where(dshcxy(:,:,:,imdl) == missing .or. shcxy == missing) dshcxy(:,:,:,imdl) = 0.
where(dlwcxy(:,:,:,imdl) == missing .or. lwcxy == missing) dlwcxy(:,:,:,imdl) = 0.
where(dhurcxy(:,:,:,imdl) == missing .or. hurcxy == missing) dhurcxy(:,:,:,imdl) = 0.

!** the following where statements have to be put after the previous where statements
where(pcxy == missing) pcxy = 1.
where(tascxy == missing) tascxy = 1.
where(tmxcxy == missing) tmxcxy = 1.
where(tmncxy == missing) tmncxy = 1.
where(swcxy == missing) swcxy = 1.
where(wscxy == missing) wscxy = 1.
where(shcxy == missing) shcxy = 1.
where(lwcxy == missing) lwcxy = 1.
where(hurcxy == missing) hurcxy = 1.

ENDDO !** end loop of reading HFD coefficients for all the models

DO ie = 1, nens
  write(ce3,'(I3.3)') ie
  write(che3,'(I4)') enum(1)+ie
  write(cfe3,'(I4)') enum(2)+ie
!** (2012-2100) to calculate the  80-year (2021-2100) moving temperature change relative to 2012-2021 
  fin1='/net/fs10/d1/sokolov/MESM2.2/20C/ENSM_ALL50_XIANG/'//che3//'.21/plot'//che3//'.21'
!  if(ie == 1) write(*,*) fin1
  im = 1
  icnt = 0
  icnt0 = 0
  gt = missing
  dgt = missing
!** 100-year future zonal values
  zp0 = missing
  zt0 = missing
  ztmx0 = missing
  ztmn0 = missing
  zsw0 = missing
  zlw0 = missing
  zws0 = missing
  zsh0 = missing
  zrh0 = missing
!** 80-yr moving-average future zonal values
  zp = missing
  zt = missing
  ztmx = missing
  ztmn = missing
  zsw = missing
  zlw = missing
  zws = missing
  zsh = missing
  zrh = missing
!** 2001-2020 zonal values
  mzp = missing
  mzt = missing
  mztmx = missing
  mztmn = missing
  mzsw = missing
  mzlw = missing
  mzws = missing
  mzsh = missing
  mzrh = missing
  open(11,file=fin1,status='old',form='unformatted')
  read(11) dumexp,dumtau,dumtbl
  do while(flag)
      read(11,end=99) dumexp,dumd1,dumc4m1,dumyr1,dumd2,dumc4m2,dumyr2,gbudg,&
                      qmaps,qtable,inqtab,j1qt,inqmap
!      write(*,*) dumyr2,dumc4m2,gbudg(49,35,1)/10
      IF(dumyr2 >= 2001 .and. dumyr2 <= 2005) then
!         IF(ie == 1) write(*,*) dumyr2,im
         zp0(:,dumyr2-2001+1,im) = gbudg(1:46,50,1)            
         zt0(:,dumyr2-2001+1,im) = gbudg(1:46,35,1)/10. + 273.15
         ztmx0(:,dumyr2-2001+1,im) = gbudg(1:46,31,1)/10. + 273.15
         ztmn0(:,dumyr2-2001+1,im) = gbudg(1:46,32,1)/10. + 273.15
         zsw0(:,dumyr2-2001+1,im) = gbudg(1:46,5,1)
         zlw0(:,dumyr2-2001+1,im) = gbudg(1:46,10,1)
         zws0(:,dumyr2-2001+1,im) = gbudg(1:46,39,1)
!** specific humidity
         zsh0(:,dumyr2-2001+1,im) = qtable(1:46,1,3)*0.00001
!** relative humidity (0 - 100)
         zrh0(:,dumyr2-2001+1,im) = gbudg(1:46,33,1)
         im = im + 1
         icnt = icnt + 1
      ELSE
      ENDIF
      if(im > 12) im = 1
  enddo
99 continue
!  IF(ie == 1) write(*,*) 'icnt= ',icnt !** icnt should be  5*12 =  60
!  IF(ie == 1) write(*,*) 'im = ', im  !** should be 1
  close(11)

  fin2 ='/net/fs10/d1/sokolov/MESM2.2/21C/OUTLOOK2020_ENMS/'//sce(1:len_trim(sce))//'_50_XIANG/'//cfe3//'.21/plot'//cfe3//'.21'
!  IF(ie == 1) write(*,*) fin1,fin2
  open(11,file=fin2,status='old',form='unformatted')
  read(11) dumexp,dumtau,dumtbl
  do while(flag)
     read(11,end=101) dumexp,dumd1,dumc4m1,dumyr1,dumd2,dumc4m2,dumyr2,gbudg,&
                      qmaps,qtable,inqtab,j1qt,inqmap
!     write(*,*) dumyr2,dumc4m2,gbudg(49,35,1)/10
     IF(dumyr2 >= 2006 .and. dumyr2 <= yr0+nyr) THEN
!       IF(ie == 1) write(*,*) dumyr2,im
       zp0(:,dumyr2-2001+1,im) = gbudg(1:46,50,1)            
       zt0(:,dumyr2-2001+1,im) = gbudg(1:46,35,1)/10. + 273.15
       ztmx0(:,dumyr2-2001+1,im) = gbudg(1:46,31,1)/10. + 273.15
       ztmn0(:,dumyr2-2001+1,im) = gbudg(1:46,32,1)/10. + 273.15
       zsw0(:,dumyr2-2001+1,im) = gbudg(1:46,5,1)
       zlw0(:,dumyr2-2001+1,im) = gbudg(1:46,10,1)
       zws0(:,dumyr2-2001+1,im) = gbudg(1:46,39,1)
!** specific humidity
       zsh0(:,dumyr2-2001+1,im) = qtable(1:46,1,3)*0.00001
!** relative humidity (0 - 100)
       zrh0(:,dumyr2-2001+1,im) = gbudg(1:46,33,1)
       icnt = icnt + 1
       IF(dumyr2 >= 2012) THEN
         gt(dumyr2-2012+1,im) = gbudg(49,35,1)/10.
         icnt0 = icnt0 + 1
       ENDIF  
       im = im + 1
     else
       go to 100      
     endif
     if(im > 12) im = 1
  enddo
100 continue  
101 continue
!  IF(ie == 1) write(*,*) 'icnt= ',icnt !** icnt should be 100*12 =  1200
!  IF(ie == 1) write(*,*) 'icnt0= ',icnt0 !** icnt should be 89*12 = 1068
!  IF(ie == 1) write(*,*) 'im = ', im  !** should be 1
  close(11)

!** Calculate the zonal IGSM climatology from 2001-2020
  mzp(:,:) = sum(zp0(:,1:20,:),2)/20.
  mzt(:,:) = sum(zt0(:,1:20,:),2)/20.
  mztmx(:,:) = sum(ztmx0(:,1:20,:),2)/20.
  mztmn(:,:) = sum(ztmn0(:,1:20,:),2)/20.
  mzsw(:,:) = sum(zsw0(:,1:20,:),2)/20.
  mzlw(:,:) = sum(zlw0(:,1:20,:),2)/20.
  mzws(:,:) = sum(zws0(:,1:20,:),2)/20.
  mzsh(:,:) = sum(zsh0(:,1:20,:),2)/20.
  mzrh(:,:) = sum(zrh0(:,1:20,:),2)/20.

!** calculate the difference between 10-yr running mean of global mean temeperature and  the 2012-2021 mean for each month
  DO im = 1, nm
    DO iy = 1, nyr
      sumt = 0.0
      icnt = 0
      DO i = iy, iy+nyrc-1
        sumt = sumt + gt(i,im)
        icnt = icnt + 1
!        IF(ie == 1 .and. im == 1) write (*,*) iy, i, icnt
      ENDDO
      IF(iy == 1) baset = sumt/real(icnt)
      dgt(iy,im) = sumt/real(icnt) - baset
!      IF(ie == 1) write(*,'(4I4,2F10.3)') im,iy,iy+nyrc-1,icnt,baset,dgt(iy,im)
    ENDDO
  ENDDO

!** calculate the 2021-2100 20-yr moving average zonal values
  DO im = 1, nm
    DO iiy = 1, nzy 
      DO iy = 1, nyr
        summp = 0.0
        summt = 0.0
        summtx = 0.0
        summtn = 0.0
        summsw = 0.0
        summlw = 0.0
        summws = 0.0
        summsh = 0.0
        summrh = 0.0
        icnt = 0
        DO i = iy+1, iy+nyrz
          summp = summp + zp0(iiy,i,im)
          summt = summt + zt0(iiy,i,im)
          summtx = summtx + ztmx0(iiy,i,im)
          summtn = summtn + ztmn0(iiy,i,im)
          summsw = summsw + zsw0(iiy,i,im)
          summlw = summlw + zlw0(iiy,i,im)
          summws = summws + zws0(iiy,i,im)
          summsh = summsh + zsh0(iiy,i,im)
          summrh = summrh + zrh0(iiy,i,im)
          icnt = icnt + 1
!        IF(ie == 1 .and. im == 1 .and. iiy == 20) write (*,*) iy, i, icnt
        ENDDO
        zp(iiy,iy,im) = summp/real(icnt)
        zt(iiy,iy,im) = summt/real(icnt)
        ztmx(iiy,iy,im) = summtx/real(icnt)
        ztmn(iiy,iy,im) = summtn/real(icnt)
        zsw(iiy,iy,im) = summsw/real(icnt)
        zlw(iiy,iy,im) = summlw/real(icnt)
        zws(iiy,iy,im) = summws/real(icnt)
        zsh(iiy,iy,im) = summsh/real(icnt)
        zrh(iiy,iy,im) = summrh/real(icnt)
      ENDDO  
    ENDDO
  ENDDO  

!** start the model loop for downscaling
  DO imdl = 1, nmdl
!** define the monthly netcdf output for 2011-2150 
    fdir = '/net/fs04/d2/xgao/pnnl/d0.5/'//sce(1:len_trim(sce))//'/'//cmodel(imdl)(1:len_trim(cmodel(imdl)))
    fout1 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_Prep_e'//ce3//'_monthly.nc'
    fout2 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_Tair_e'//ce3//'_monthly.nc'
    fout3 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_Tmax_e'//ce3//'_monthly.nc'
    fout4 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_Tmin_e'//ce3//'_monthly.nc'
    fout5 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_WIND_e'//ce3//'_monthly.nc'
    fout6 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_FSDS_e'//ce3//'_monthly.nc'
    fout7 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_FLDS_e'//ce3//'_monthly.nc'
    fout8 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_Huss_e'//ce3//'_monthly.nc'
    fout9 = fdir(1:len_trim(fdir))//'/2021_2100_0.5_Hurs_e'//ce3//'_monthly.nc'
!    IF(ie == 1) write(*,*) fout
!** enter define mode
!** precipitation
    iret = NF90_CREATE(fout1,NF90_CLOBBER, oid1)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid1,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid1,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid1, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid1,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid1,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,latid,'axis','Y')
    CALL CHECK_err(iret)    

    iret = NF90_DEF_VAR(oid1,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,timeid,'axis','T')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid1,'PRECTmmd',NF90_FLOAT,three_dims,prepid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,prepid,'long_name','Precipitation rate')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,prepid,'standard_name','precipitation_flux')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,prepid,'units','mm/day')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,prepid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,prepid,'mode','time-dependent')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid1,'PRECTmmd_trend',NF90_FLOAT,three_dims,trprepid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,trprepid,'long_name','IGSM precipitation rate trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,trprepid,'units','mm/day')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,trprepid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid1,trprepid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid1, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively' 
    iret = NF90_PUT_ATT(oid1, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid1, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid1, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid1)
    CALL CHECK_err(iret)

!** air temperature
    iret = NF90_CREATE(fout2,NF90_CLOBBER, oid2)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid2,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid2,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid2, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid2,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid2,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,latid,'axis','Y')
    CALL CHECK_err(iret)    

    iret = NF90_DEF_VAR(oid2,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,timeid,'axis','T')
    CALL CHECK_err(iret)
    
    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid2,'Tair',NF90_FLOAT,three_dims,tairid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,tairid,'long_name','Near surface air temperature')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,tairid,'standard_name','air_temperature')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,tairid,'units','K')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,tairid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,tairid,'mode','time-dependent')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid2,'Tair_trend',NF90_FLOAT,three_dims,trtairid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,trtairid,'long_name','IGSM near surface air temperature trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,trtairid,'units','K')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,trtairid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid2,trtairid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid2, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively'
    iret = NF90_PUT_ATT(oid2, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid2, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid2, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid2)
    CALL CHECK_err(iret)

!** Maximum Temperature
    iret = NF90_CREATE(fout3,NF90_CLOBBER, oid3)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid3,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid3,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid3, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid3,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid3,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,latid,'axis','Y')
    CALL CHECK_err(iret)    
    
    iret = NF90_DEF_VAR(oid3,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,timeid,'axis','T')
    CALL CHECK_err(iret)
    
    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid3,'Tmax',NF90_FLOAT,three_dims,tmaxid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,tmaxid,'long_name','Monthly mean of daily maximum near surface air temperature')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,tmaxid,'standard_name','air_temperature')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,tmaxid,'units','K')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,tmaxid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,tmaxid,'mode','time-dependent')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid3,'Tmax_trend',NF90_FLOAT,three_dims,trtmaxid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,trtmaxid,'long_name','IGSM monthly mean of daily maximum near surface air temperature trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,trtmaxid,'units','K')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,trtmaxid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid3,trtmaxid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid3, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively'
    iret = NF90_PUT_ATT(oid3, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid3, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid3, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid3)
    CALL CHECK_err(iret)

!** Minimum Temperature
    iret = NF90_CREATE(fout4,NF90_CLOBBER, oid4)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid4,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid4,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid4, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid4,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid4,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,latid,'axis','Y')
    CALL CHECK_err(iret)    
    
    iret = NF90_DEF_VAR(oid4,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,timeid,'axis','T')
    CALL CHECK_err(iret)
    
    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid4,'Tmin',NF90_FLOAT,three_dims,tminid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,tminid,'long_name','Monthly mean of daily minimum near surface air temperature')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,tminid,'standard_name','air_temperature')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,tminid,'units','K')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,tminid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,tminid,'mode','time-dependent')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid4,'Tmin_trend',NF90_FLOAT,three_dims,trtminid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,trtminid,'long_name','IGSM monthly mean of daily minimum near surface air temperature trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,trtminid,'units','K')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,trtminid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid4,trtminid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid4, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively'
    iret = NF90_PUT_ATT(oid4, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid4, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid4, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid4)
    CALL CHECK_err(iret)

!** Wind Speed
    iret = NF90_CREATE(fout5,NF90_CLOBBER, oid5)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid5,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid5,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid5, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid5,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid5,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,latid,'axis','Y')
    CALL CHECK_err(iret)    

    iret = NF90_DEF_VAR(oid5,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,timeid,'axis','T')
    CALL CHECK_err(iret)
    
    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid5,'WIND',NF90_FLOAT,three_dims,windid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,windid,'long_name','Near surface wind speed')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,windid,'standard_name','wind_speed')
    CALL CHECK_err(iret)    
    iret = NF90_PUT_ATT(oid5,windid,'units','m/s')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,windid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,windid,'mode','time-dependent')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid5,'WIND_trend',NF90_FLOAT,three_dims,trwindid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,trwindid,'long_name','IGSM near surface wind speed trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,trwindid,'units','m/s')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,trwindid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid5,trwindid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid5, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively'
    iret = NF90_PUT_ATT(oid5, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid5, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid5, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid5)
    CALL CHECK_err(iret)

!** shortwave radiation
    iret = NF90_CREATE(fout6,NF90_CLOBBER, oid6)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid6,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid6,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid6, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid6,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid6,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,latid,'axis','Y')
    CALL CHECK_err(iret)    

    iret = NF90_DEF_VAR(oid6,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,timeid,'axis','T')
    CALL CHECK_err(iret)
    
    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid6,'FSDS',NF90_FLOAT,three_dims,swdownid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,swdownid,'long_name','Surface incident shortwave radiation')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,swdownid,'standard_name','surface_downwelling_shortwave_flux_in_air')
    CALL CHECK_err(iret)    
    iret = NF90_PUT_ATT(oid6,swdownid,'units','W/m2')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,swdownid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,swdownid,'mode','time-dependent')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid6,'FSDS_trend',NF90_FLOAT,three_dims,trswdownid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,trswdownid,'long_name','IGSM surface incident shortwave radiation trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,trswdownid,'units','W/m2')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,trswdownid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid6,trswdownid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid6, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively'
    iret = NF90_PUT_ATT(oid6, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid6, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid6, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid6)
    CALL CHECK_err(iret)

!** longwave radiation
    iret = NF90_CREATE(fout7,NF90_CLOBBER, oid7)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid7,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid7,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid7, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid7,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid7,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,latid,'axis','Y')
    CALL CHECK_err(iret)    
    
    iret = NF90_DEF_VAR(oid7,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,timeid,'axis','T')
    CALL CHECK_err(iret)
    
    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid7,'FLDS',NF90_FLOAT,three_dims,lwdownid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,lwdownid,'long_name','Surface incident longwave radiation')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,lwdownid,'standard_name','surface_downwelling_longwave_flux_in_air')
    CALL CHECK_err(iret)    
    iret = NF90_PUT_ATT(oid7,lwdownid,'units','W/m2')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,lwdownid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,lwdownid,'mode','time-dependent')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid7,'FLDS_trend',NF90_FLOAT,three_dims,trlwdownid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,trlwdownid,'long_name','IGSM surface incident longwave radiation trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,trlwdownid,'units','W/m2')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,trlwdownid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid7,trlwdownid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid7, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively'
    iret = NF90_PUT_ATT(oid7, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid7, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid7, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid7)
    CALL CHECK_err(iret)

!** Specific Humidity
    iret = NF90_CREATE(fout8,NF90_CLOBBER, oid8)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid8,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid8,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid8, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid8,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid8,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,latid,'axis','Y')
    CALL CHECK_err(iret)    

    iret = NF90_DEF_VAR(oid8,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,timeid,'axis','T')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid8,'Huss',NF90_FLOAT,three_dims,shid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,shid,'long_name','Near surface specific humidity')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,shid,'standard_name','surface_specific_humidity')
    CALL CHECK_err(iret)    
    iret = NF90_PUT_ATT(oid8,shid,'units','kg/kg')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,shid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,shid,'mode','time-dependent')
    CALL CHECK_err(iret)    

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid8,'Huss_trend',NF90_FLOAT,three_dims,trshid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,trshid,'long_name','IGSM near surface specific humidity trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,trshid,'units','kg/kg')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,trshid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid8,trshid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid8, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively'
    iret = NF90_PUT_ATT(oid8, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid8, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid8, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid8)
    CALL CHECK_err(iret)

!** relative humidity 
    iret = NF90_CREATE(fout9,NF90_CLOBBER, oid9)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid9,'lon',nx,londim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid9,'lat',ny,latdim)
    CALL CHECK_err(iret)
    iret = NF90_DEF_DIM(oid9, 'time', NF90_UNLIMITED, time_dim)
    CALL CHECK_err(iret)

    one_dim(1) = londim
    iret = NF90_DEF_VAR(oid9,'lon',NF90_FLOAT,one_dim,lonid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,lonid,'standard_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,lonid,'long_name','longitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,lonid,'units','degrees_east')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,lonid,'axis','X')
    CALL CHECK_err(iret)    

    one_dim(1) = latdim
    iret = NF90_DEF_VAR(oid9,'lat',NF90_FLOAT,one_dim,latid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,latid,'standard_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,latid,'long_name','latitude')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,latid,'units','degrees_north')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,latid,'axis','Y')
    CALL CHECK_err(iret)    

    iret = NF90_DEF_VAR(oid9,'time',NF90_FLOAT,time_dim,timeid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,timeid,'standard_name','time')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,timeid,'units','months since 2021-01-01 00:00:00')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,timeid,'calendar','proleptic_gregorian')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,timeid,'axis','T')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid9,'Hurs',NF90_FLOAT,three_dims,rhid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,rhid,'long_name','Near surface relative humidity')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,rhid,'standard_name','relative_humidity')
    CALL CHECK_err(iret)    
    iret = NF90_PUT_ATT(oid9,rhid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,rhid,'mode','time-dependent')
    CALL CHECK_err(iret)

    three_dims(3) = time_dim
    three_dims(2) = latdim
    three_dims(1) = londim
    iret = NF90_DEF_VAR(oid9,'Hurs_trend',NF90_FLOAT,three_dims,trrhid)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,trrhid,'long_name','IGSM near surface relative humidity trend (delta)')
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,trrhid,'_FillValue',missing)
    CALL CHECK_err(iret)
    iret = NF90_PUT_ATT(oid9,trrhid,'mode','time-dependent')
    CALL CHECK_err(iret)

    cdum = 'HFD scaled MESM 0.5 degree monthly meteorological forcing (2021-2100) for '//sce(1:len_trim(sce))
    iret = NF90_PUT_ATT(oid9, NF90_GLOBAL, 'title',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'HFD climatology and future coefficients are based on GSWP3 and CMIP6, respectively'
    iret = NF90_PUT_ATT(oid9, NF90_GLOBAL, 'remarks',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Created by MIT Joint Program on the Science and Policy of Global Change'
    iret = NF90_PUT_ATT(oid9, NF90_GLOBAL, 'sources',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)
    cdum = 'Xiang Gao:xgao304@mit.edu'
    iret = NF90_PUT_ATT(oid9, NF90_GLOBAL, 'contact',cdum(1:len_trim(cdum)))
    CALL CHECK_err(iret)

    iret = NF90_ENDDEF(oid9)
    CALL CHECK_err(iret)

!** write out time-independent variables
    iret = NF90_PUT_VAR(oid1,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid1,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid1,timeid,tstep)
    CALL CHECK_err(iret)

    iret = NF90_PUT_VAR(oid2,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid2,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid2,timeid,tstep)
    CALL CHECK_err(iret)

    iret = NF90_PUT_VAR(oid3,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid3,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid3,timeid,tstep)
    CALL CHECK_err(iret)

    iret = NF90_PUT_VAR(oid4,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid4,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid4,timeid,tstep)
    CALL CHECK_err(iret)

    iret = NF90_PUT_VAR(oid5,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid5,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid5,timeid,tstep)
    CALL CHECK_err(iret)

    iret = NF90_PUT_VAR(oid6,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid6,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid6,timeid,tstep)
    CALL CHECK_err(iret)

    iret = NF90_PUT_VAR(oid7,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid7,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid7,timeid,tstep)
    CALL CHECK_err(iret)

    iret = NF90_PUT_VAR(oid8,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid8,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid8,timeid,tstep)
    CALL CHECK_err(iret)

    iret = NF90_PUT_VAR(oid9,lonid,lon)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid9,latid,lat)
    CALL CHECK_err(iret)
    iret = NF90_PUT_VAR(oid9,timeid,tstep)
    CALL CHECK_err(iret)    

!** this is the total fields      
    prep = missing
    tair = missing
    tmax = missing
    tmin = missing
    sw = missing
    lw = missing
    wind = missing
    huss = missing 
    rh = missing
!** this is the delta fields (IGSM trend) 
    dlt_prep = missing
    dlt_tair = missing
    dlt_tmax = missing
    dlt_tmin = missing
    dlt_sw = missing
    dlt_lw = missing
    dlt_wind = missing
    dlt_huss = missing
    dlt_rh = missing
    DO iy = 1, nyr
      DO im = 1, nm
        DO iry = 1,ny
          IF(iry <= 4) THEN
            jindex = 1
          ELSE  
            IF(mod(iry-4,8) == 0) THEN
              jindex = 1 + int((iry-4)/8)
            ELSE
              jindex = 2 + int((iry-4)/8)
            ENDIF
          ENDIF  
!          IF(ie == 1 .and. imdl == 1 .and. iy == 1 .and. im == 1) write(*,*) iry, jindex
          DO irx = 1,nx
            tpxy(irx,iry) = pcxy(irx,iry,im)+dgt(iy,im)*dpcxy(irx,iry,im,imdl)
            ttxy(irx,iry) = tascxy(irx,iry,im)+dgt(iy,im)*dtascxy(irx,iry,im,imdl)
            ttmxxy(irx,iry) = tmxcxy(irx,iry,im)+dgt(iy,im)*dtmxcxy(irx,iry,im,imdl)
            ttmnxy(irx,iry) = tmncxy(irx,iry,im)+dgt(iy,im)*dtmncxy(irx,iry,im,imdl)
            tswxy(irx,iry) = swcxy(irx,iry,im)+dgt(iy,im)*dswcxy(irx,iry,im,imdl)
            tlwxy(irx,iry) = lwcxy(irx,iry,im)+dgt(iy,im)*dlwcxy(irx,iry,im,imdl)
            twsxy(irx,iry) = wscxy(irx,iry,im)+dgt(iy,im)*dwindcxy(irx,iry,im,imdl)
            tshxy(irx,iry) = shcxy(irx,iry,im)+dgt(iy,im)*dshcxy(irx,iry,im,imdl)
            thurxy(irx,iry) = hurcxy(irx,iry,im)+dgt(iy,im)*dhurcxy(irx,iry,im,imdl)
            IF(tpxy(irx,iry) < 0.001) tpxy (irx,iry) = 0.001
            IF(ttxy(irx,iry) < 0.001) ttxy (irx,iry) = 0.001
            IF(ttmxxy(irx,iry) < 0.001) ttmxxy (irx,iry) = 0.001
            IF(ttmnxy(irx,iry) < 0.001) ttmnxy (irx,iry) = 0.001
            IF(tswxy(irx,iry) < 0.001) tswxy (irx,iry) = 0.001
            IF(tlwxy(irx,iry) < 0.001) tlwxy (irx,iry) = 0.001
            IF(twsxy(irx,iry) < 0.001) twsxy (irx,iry) = 0.001
            IF(tshxy(irx,iry) < 0.001) tshxy (irx,iry) = 0.001
            IF(thurxy(irx,iry) < 0.001) thurxy (irx,iry) = 0.001
!** derived IGSM trend
            dlt_prep(irx,iry,iy,im) = tpxy(irx,iry)*zp(jindex,iy,im)-pcxy(irx,iry,im)*mzp(jindex,im)
            dlt_tair(irx,iry,iy,im) = ttxy(irx,iry)*zt(jindex,iy,im)-tascxy(irx,iry,im)*mzt(jindex,im)
            dlt_tmax(irx,iry,iy,im) = ttmxxy(irx,iry)*ztmx(jindex,iy,im)-tmxcxy(irx,iry,im)*mztmx(jindex,im)
            dlt_tmin(irx,iry,iy,im) = ttmnxy(irx,iry)*ztmn(jindex,iy,im)-tmncxy(irx,iry,im)*mztmn(jindex,im)
            dlt_sw(irx,iry,iy,im) = tswxy(irx,iry)*zsw(jindex,iy,im)-swcxy(irx,iry,im)*mzsw(jindex,im)
            dlt_lw(irx,iry,iy,im) = tlwxy(irx,iry)*zlw(jindex,iy,im)-lwcxy(irx,iry,im)*mzlw(jindex,im)
            dlt_wind(irx,iry,iy,im) = twsxy(irx,iry)*zws(jindex,iy,im)-wscxy(irx,iry,im)*mzws(jindex,im)
            dlt_huss(irx,iry,iy,im) = tshxy(irx,iry)*zsh(jindex,iy,im)-shcxy(irx,iry,im)*mzsh(jindex,im)
            dlt_rh(irx,iry,iy,im) = thurxy(irx,iry)*zrh(jindex,iy,im)-hurcxy(irx,iry,im)*mzrh(jindex,im)
          ENDDO
        ENDDO
      ENDDO  
    ENDDO
!** apply the 20-yr running average for IGSM trend
    DO im = 1, nm
      DO iry = 1, ny
        DO irx = 1, nx
          DO iy = 1, nyr
            summp  = 0.0           
            summt  = 0.0           
            summtx = 0.0           
            summtn = 0.0           
            summsw = 0.0           
            summlw = 0.0           
            summws = 0.0           
            summsh = 0.0           
            summrh = 0.0           
            IF(iy < 11) THEN
!** year 2021-2030
              DO i = 1,iy+10
!                IF(im == 1 .and. irx == 48 .and. iry == 41) write(*,*) iy,yr0+iy,yr0+i
                summp  = summp  + dlt_prep(irx,iry,i,im) 
                summt  = summt  + dlt_tair(irx,iry,i,im)
                summtx = summtx + dlt_tmax(irx,iry,i,im)
                summtn = summtn + dlt_tmin(irx,iry,i,im)
                summsw = summsw + dlt_sw(irx,iry,i,im)
                summlw = summlw + dlt_lw(irx,iry,i,im)
                summws = summws + dlt_wind(irx,iry,i,im)
                summsh = summsh + dlt_huss(irx,iry,i,im)
                summrh = summrh + dlt_rh(irx,iry,i,im)
              ENDDO
              dlt_prep(irx,iry,iy,im) = summp/(iy+10.0)
              dlt_tair(irx,iry,iy,im) = summt/(iy+10.0)
              dlt_tmax(irx,iry,iy,im) = summtx/(iy+10.0)
              dlt_tmin(irx,iry,iy,im) = summtn/(iy+10.0)
              dlt_sw(irx,iry,iy,im)   = summsw/(iy+10.0)
              dlt_lw(irx,iry,iy,im)   = summlw/(iy+10.0)
              dlt_wind(irx,iry,iy,im) = summws/(iy+10.0)
              dlt_huss(irx,iry,iy,im) = summsh/(iy+10.0)    
              dlt_rh(irx,iry,iy,im) = summrh/(iy+10.0)    
            ELSEIF(iy > 70) THEN
              DO i = iy-10,nyr
!                IF(im == 1 .and. irx == 48 .and. iry == 41) write(*,*) iy,yr0+iy,yr0+i
                summp  = summp  + dlt_prep(irx,iry,i,im)
                summt  = summt  + dlt_tair(irx,iry,i,im)
                summtx = summtx + dlt_tmax(irx,iry,i,im)
                summtn = summtn + dlt_tmin(irx,iry,i,im)
                summsw = summsw + dlt_sw(irx,iry,i,im)
                summlw = summlw + dlt_lw(irx,iry,i,im)
                summws = summws + dlt_wind(irx,iry,i,im)
                summsh = summsh + dlt_huss(irx,iry,i,im)
                summrh = summrh + dlt_rh(irx,iry,i,im)
              ENDDO
              dlt_prep(irx,iry,iy,im) = summp/(nyr-iy+11)
              dlt_tair(irx,iry,iy,im) = summt/(nyr-iy+11)
              dlt_tmax(irx,iry,iy,im) = summtx/(nyr-iy+11)
              dlt_tmin(irx,iry,iy,im) = summtn/(nyr-iy+11)
              dlt_sw(irx,iry,iy,im)   = summsw/(nyr-iy+11)
              dlt_lw(irx,iry,iy,im)   = summlw/(nyr-iy+11)
              dlt_wind(irx,iry,iy,im) = summws/(nyr-iy+11)
              dlt_huss(irx,iry,iy,im) = summsh/(nyr-iy+11)                    
              dlt_rh(irx,iry,iy,im) = summrh/(nyr-iy+11)                    
            ELSE
!** year 2031 ~ 2090
              DO i = iy-10,iy+10
!                IF(im == 1 .and. irx == 48 .and. iry == 41) write(*,*) iy,yr0+iy,yr0+i
                summp  = summp  + dlt_prep(irx,iry,i,im)
                summt  = summt  + dlt_tair(irx,iry,i,im)
                summtx = summtx + dlt_tmax(irx,iry,i,im)
                summtn = summtn + dlt_tmin(irx,iry,i,im)
                summsw = summsw + dlt_sw(irx,iry,i,im)
                summlw = summlw + dlt_lw(irx,iry,i,im)
                summws = summws + dlt_wind(irx,iry,i,im)
                summsh = summsh + dlt_huss(irx,iry,i,im)
                summrh = summrh + dlt_rh(irx,iry,i,im)
              ENDDO
              dlt_prep(irx,iry,iy,im) = summp/21.0
              dlt_tair(irx,iry,iy,im) = summt/21.0
              dlt_tmax(irx,iry,iy,im) = summtx/21.0
              dlt_tmin(irx,iry,iy,im) = summtn/21.0
              dlt_sw(irx,iry,iy,im)   = summsw/21.0
              dlt_lw(irx,iry,iy,im)   = summlw/21.0
              dlt_wind(irx,iry,iy,im) = summws/21.0
              dlt_huss(irx,iry,iy,im) = summsh/21.0
              dlt_rh(irx,iry,iy,im) = summrh/21.0
            ENDIF 
          ENDDO !** year loop
        ENDDO  !** latitude loop
      ENDDO !** longitude loop
    ENDDO !** month loop  
!** derived total field
    DO iy = 1, nyr
      DO im = 1, nm
        DO iry = 1,ny
          DO irx = 1,nx
            prep(irx,iry) = p0(irx,iry,(iy-1)*nm+im) + dlt_prep(irx,iry,iy,im)
            tair(irx,iry) = tas0(irx,iry,(iy-1)*nm+im) + dlt_tair(irx,iry,iy,im)
            tmax(irx,iry) = tmx0(irx,iry,(iy-1)*nm+im) + dlt_tmax(irx,iry,iy,im)
            tmin(irx,iry) = tmn0(irx,iry,(iy-1)*nm+im) + dlt_tmin(irx,iry,iy,im)
            sw(irx,iry)   = sw0(irx,iry,(iy-1)*nm+im) + dlt_sw(irx,iry,iy,im)
            lw(irx,iry)   = lw0(irx,iry,(iy-1)*nm+im) + dlt_lw(irx,iry,iy,im)
            wind(irx,iry) = ws0(irx,iry,(iy-1)*nm+im) + dlt_wind(irx,iry,iy,im)
            huss(irx,iry) = sh0(irx,iry,(iy-1)*nm+im) + dlt_huss(irx,iry,iy,im)
            rh(irx,iry) = hur0(irx,iry,(iy-1)*nm+im) + dlt_rh(irx,iry,iy,im)
            IF(prep(irx,iry) < 0.) prep(irx,iry) = 0.
            CALL RANDOM_Number(rnum)            
            IF(rh(irx,iry) > 100.) THEN
              rh(irx,iry) = 95. + rnum*rndns      
            ELSEIF(rh(irx,iry) < 0.) THEN
              rh(irx,iry) = 5.0 - rnum*rndns
            ELSE
            ENDIF        
!            IF(rh(irx,iry) > 95.) rh(irx,iry) = 95.
!            IF(rh(irx,iry) <  5.) rh(irx,iry) =  5.
!** convert specific humidity to relative humidity
!** es =  6.112 * exp((17.67 * temp)/(temp + 243.5))
!** e  = qair * press / (0.378 * qair + 0.622)
!** rh = 100 * e / es
!** temp = Celsius; press = mb (hpa), 1 pa = 0.01 mbar; rh = 100* e/es
!            tmp1 = huss(irx,iry)*psrf(irx,iry,im)*0.01/(0.378*huss(irx,iry)+0.622)
!            tmp2 = 6.112 * exp((17.67 * (tair(irx,iry)-273.15))/((tair(irx,iry)-273.15) + 243.5))
!            rh(irx,iry) = 100. * tmp1 / tmp2
          ENDDO  !** longitude loop
        ENDDO  !** latitude loop
!** write out the 2D output for 2021-2100
        istart(3) = (iy-1)*nm+im
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = ny
        icount(1) = nx
!** write out total fields        
        iret = NF90_PUT_VAR(oid1,prepid,prep,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid1,trprepid,dlt_prep(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)

        iret = NF90_PUT_VAR(oid2,tairid,tair,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid2,trtairid,dlt_tair(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)

        iret = NF90_PUT_VAR(oid3,tmaxid,tmax,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid3,trtmaxid,dlt_tmax(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)

        iret = NF90_PUT_VAR(oid4,tminid,tmin,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid4,trtminid,dlt_tmin(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)

        iret = NF90_PUT_VAR(oid5,windid,wind,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid5,trwindid,dlt_wind(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)

        iret = NF90_PUT_VAR(oid6,swdownid,sw,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid6,trswdownid,dlt_sw(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)

        iret = NF90_PUT_VAR(oid7,lwdownid,lw,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid7,trlwdownid,dlt_lw(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)

        iret = NF90_PUT_VAR(oid8,shid,huss,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid8,trshid,dlt_huss(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)

        iret = NF90_PUT_VAR(oid9,rhid,rh,istart,icount)
        CALL CHECK_err(iret)
        iret = NF90_PUT_VAR(oid9,trrhid,dlt_rh(:,:,iy,im),istart,icount)
        CALL CHECK_err(iret)
      ENDDO  !** month loop
    ENDDO !** year loop  
    iret = NF90_CLOSE(oid1)
    CALL CHECK_err(iret)
    iret = NF90_CLOSE(oid2)
    CALL CHECK_err(iret)
    iret = NF90_CLOSE(oid3)
    CALL CHECK_err(iret)
    iret = NF90_CLOSE(oid4)
    CALL CHECK_err(iret)
    iret = NF90_CLOSE(oid5)
    CALL CHECK_err(iret)
    iret = NF90_CLOSE(oid6)
    CALL CHECK_err(iret)
    iret = NF90_CLOSE(oid7)
    CALL CHECK_err(iret)
    iret = NF90_CLOSE(oid8)
    CALL CHECK_err(iret)
    iret = NF90_CLOSE(oid9)
    CALL CHECK_err(iret)
  ENDDO !** model loop
ENDDO !** ensemble loop
    
END PROGRAM genmon

SUBROUTINE CHECK_ERR(iret)

USE NETCDF
IMPLICIT NONE
  
INTEGER                :: iret
  
IF(iret /= NF90_NOERR) THEN
  PRINT*, NF90_STRERROR(iret)
  STOP 'check_err'
ENDIF

END SUBROUTINE CHECK_ERR

