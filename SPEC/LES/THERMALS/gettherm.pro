PRO gettherm

spawn, 'clear'
print, ''
print, '** Thermals Analysis **'
print, ' (usine à gaz) '
print, ''

full='true'
f_offset='false'
overplot_convadj='false'
plot_3d = 'false'

lctu_gcm = 8.

; *********** Best values for LES thermals with tau =0.5
betalpha = 1.3
afact = 1.8
fact_epsilon = 0.0008
detr_min = 0.0007


;betalpha = 1.
;afact = 2.4
;fact_epsilon = 0.0007
;detr_min = 0.0007

; --------

; *********** Best values for LES thermals with tau =1.
;betalpha = 1.
;afact = 1.1
;fact_epsilon = 0.00025
;detr_min = 0.0007
; --------

;datname='thermiques.dat.scale1.2'
;datname='thermiques.dat.scale1.4'
;datname='thermiques.dat.scale0.6'
datname='thermiques.dat'
;datname='thermiques.dat.scale0.8'
;datname='thermiques.dat'           ;scale =1.0, sigmao =1.0

ns = 0     ; number of points for time-smoothing of LES data : 2*ns+1 points, ns = 9 eq to 30mn (-15mn//+15mn)
nstot = float(2.*ns+1.)

GcmSubCase = ''
LayerCase=''
s_trac1 = 'qtrac1'
s_trac2 = 'qtrac2'
got_pdt = 'true'
TestCase = 'Case_A'
SubCase = '_11_shorter'
Histo = 'true'
newtest = ''
visualization_mode = 'false'
label_init:
spawn, 'clear'
print, ' Available simulations :'
print, '----'
print, ' 1/ Case_A_11 : 45x45x71,   ztop=10km,dx=100m,dz=140m,Ls=47.1°,(21.8N;205.0E),55tiu,A=0.275,Tau=0.5'
print, ' 2/ Case_A_4  : 101x101x201 ztop=15km,dx=100m,dz= 75m,Ls=47.1°,(21.8N;205.0E),55tiu,A=0.275,Tau=0.5'
print, ' 3/ Case_A_4_shorter  : Case A4 with Dtrac1 = 5 mn and Dtrac2 = 10 mn (compared to 20 and 100)'
print, ' 4/ Case_A_4_shorter_winds  : Case A4_shorter with bckgrnd wind u=10 m/s'
print, ' 5/ Case_A_4_shorter_winds_30  : Case A4_shorter with bckgrnd wind u=30 m/s'
print, ' 6/ Case_A_4_shorter_winds_tau1  : Case A4_shorter with bckgrnd wind u=10 m/s and tau=1'
print, ' 7/ Case_A_4_shorter_winds_tau2  : Case A4_shorter with bckgrnd wind u=10 m/s and tau=2'
print, ' 8/ Case_ExtremeCase  : 101x101x201 ztop=15km,dx=100m,dz=75m,Ls=0°,(0.N;0.E),50tiu,A=0.1,Tau=0.05'
print, '----'
print, ' 9/ Case_C_4_shorter_winds  : '
print, ' 10/ Case_I_4_shorter_winds  : '
print, ' 11/ Case_Z_4_shorter_winds  : '
print, ' 12/ 1D : 124 layers'
print, ' 13/ 1D : 32 layers'
print, ' 14/ 1D : low dt'
print, ''
print, ' 0/ PLUME VISUALISATION : '+visualization_mode
print, ' 999/ CLEAR thermiques.dat for considered case'
print, ' 100/ OVERPLOT CONVADJ ONLY RESULTS : '+overplot_convadj
print, ''
print, ' SIMULATION NUMBER  : '
print, ' ** '+TestCase+SubCase+LayerCase+' ** '
print, ''
print, 'Any change ? (number of new case to change, or any other key to continue)'
read, newtest
if (newtest eq '1') then begin
TestCase = 'Case_A'
SubCase = '_11'
pGround = 867.5594
goto,label_init
endif
if (newtest eq '2') then begin
TestCase = 'Case_A'
SubCase = '_4'
got_pdt = 'false'
pGround = 867.5594
goto,label_init
endif
if (newtest eq '3') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter'
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 867.5594
goto,label_init
endif
if (newtest eq '4') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter_wind'
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 867.5594
goto,label_init
endif
if (newtest eq '5') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter_wind_30'
GcmSubCase = '_wind_30'
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
pGround = 867.5594
f_offset='false'
goto,label_init
endif
if (newtest eq '6') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter_wind_tau1'
GcmSubCase = '_tau1'
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 867.5594
goto,label_init
endif
if (newtest eq '7') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter_wind_tau2'
GcmSubCase = '_tau2'
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 867.5594
goto,label_init
endif
if (newtest eq '8') then begin
TestCase = 'ExtremeCase'
SubCase = ''
GcmSubCase = ''
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 677.722
;lctu_gcm = 6.
goto,label_init
endif
if (newtest eq '9') then begin
TestCase = 'Case_C'
SubCase = '_4_shorter_wind'
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 483.
goto,label_init
endif
if (newtest eq '10') then begin
TestCase = 'Case_I'
SubCase = '_4_shorter_wind'
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 630.
goto,label_init
endif
if (newtest eq '11') then begin
TestCase = 'Case_Z'
SubCase = '_4_shorter_wind'
;s_trac1 = 'qtrac2'
;s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 266.
goto,label_init
endif
if (newtest eq '12') then begin
LayerCase=''
goto,label_init
endif
if (newtest eq '13') then begin
LayerCase='_32lev'
goto,label_init
endif
if (newtest eq '14') then begin
LayerCase='_low_dt'
goto,label_init
endif
if (newtest eq '100') then begin
if (overplot_convadj eq 'true') then overplot_convadj = 'false' else overplot_convadj = 'true'
goto,label_init
endif
if (newtest eq '0') then begin
visualization_mode = 'true'
spawn, 'clear'
print, 'The first timestep of the considered file will be shown.'
print, 'Defaut file is : file 6 (lt ~ 13h10)'
print, 'Which file do you want ? (lt ~= file_number + 7)'
loop_special = '6'
read, loop_special
print, 'Do you wish to plot histograms or manipulate volumic data ?'
print, '1/ Histogram'
print, '2/ Volumic data'
x=''
read, x
if (x eq '1') then Histo='true' else Histo='false'
goto,label_init
endif

;les_path='/san0/acolmd/SIMUS/LES_'+TestCase+SubCase
;gcm_path='/san0/acolmd/SIMUS/GCM_'+TestCase+LayerCase+'_2'
;gcm_path='/san0/acolmd/SIMUS/GCM_'+TestCase+GcmSubCase+LayerCase
;gcm_convadj_path=gcm_path+'_convadj'

;les_path='/data/acolmd/Thermiques/LES_Case_A_4_2trac_wind10_tau05/'
;les_path='/data/acolmd/Thermiques/LES_Case_A_257x257x301_wind10_tau05'
;les_path='/data/acolmd/Thermiques/LES_Case_E_257x257x301_wind10_tau005'
;les_path='/data/acolmd/Thermiques/LES_Case_A_4_2trac_wind30_tau05/'
;les_path='/data/acolmd/Thermiques/LES_Case_A_4_2trac_wind10_tau1/'
;les_path='/data/acolmd/Thermiques/LES_Case_A_4_2trac_wind10_tau2/'
;les_path='/data/acolmd/Thermiques/LES_Case_C_4_2trac_wind10_tau05/'
;les_path='/data/acolmd/Thermiques/LES_Case_I_4_2trac_wind10_tau05/'
;les_path='/data/acolmd/Thermiques/LES_Case_Z_4_2trac_wind10_tau05/'
;les_path='/data/acolmd/Thermiques/ExtremeCase/'

;les_path='/data/acolmd/Thermiques/LES_Case_A_101x101x201_tracup_wind10_tau05_gcmsoil/'
;les_path='/data/acolmd/Thermiques/LES_Case_E_101x101x201_tracup_wind10_tau005_gcmsoil/'
;les_path='/data/acolmd/Thermiques/LES_Case_C_101x101x201_tracup_wind10_tau05_gcmsoil'
;les_path='/data/acolmd/Thermiques/LES_Case_I_101x101x201_tracup_wind10_tau05_gcmsoil'
;les_path='/data/acolmd/Thermiques/LES_Case_Z_101x101x201_tracup_wind10_tau05_gcmsoil'
;les_path='/data/acolmd/Thermiques/LES_Case_A_101x101x201_tracup_wind30_tau05_gcmsoil'
les_path='/data/acolmd/Thermiques/LES_Case_A_101x101x201_tracup_wind10_tau2_gcmsoil'

gcm_path='/data/acolmd/Thermiques/THgcm/'
gcm_convadj_path=gcm_path+'_convadj'


if (newtest eq '999') then spawn, 'rm -f '+les_path+'/'+datname
print, ''
print, ' -- Loading LES data -- '

print, 'LES DATA IN : '
print, les_path
print, 'GCM DATA IN : '
print, gcm_path

p0=610. & t0=220. & r_cp=1./3.89419 & grav=3.72 & R=191.182
history_interval_s = 100.
;lctu_gcm = 8.                      ; Initial local time of gcm 1d simu
scale = 1.                         ; Scaling factor for conditional sampling
decimate = 10.			   ; Coeff for subsampling the data for sigma integral
sigmao= 0.3			   ; multiplicative coeff for the computation of Sigma0 in the CS
sigmao_ude = 0.2		   ; number of standard deviation away from mean for the selection of downdraft in UDE
NBINS=100.

openr,unit,les_path+'/'+datname,/get_lun,error=err
IF (err ne 0) THEN BEGIN

OPENR, 22, les_path+'/input_coord' & READF, 22, lonu & READF, 22, latu & READF, 22, lsu & READF, 22, lctu & CLOSE, 22
OPENR, 23, les_path+'/input_more' & READF, 23, hgtu, tsurfu & CLOSE, 23

domain='d01' & filesWRF = FindFile(les_path+'/wrfout_'+domain+'_????-??-??_??:??:??') & nf=n_elements(filesWRF)
;print, filesWRF
;domain='d01' & filesWRF = les_path+'/'+['wrfout_d01_9999-01-01_03:05:00','wrfout_d01_9999-01-01_04:06:40'] & nf=n_elements(filesWRF)
;ce fichier utilise aussi offset_localtime = 3.05

; WARNING WARNING : FOR THE CASE 4_SHORTER, THE TRAC2 HAS 10 MN LIFETIME, WE WANT TO USE IT MORE EXTENSIVELY THAN THE TRAC1 (5mn) SO
; we switch the names of trac1 and trac2 in the initialization of this routine, in the "case".

; WARNING WARNING : in this version (thermiques2), one of the tracers in the data is just a neutral tracer which will be used to compute the turbulent flux of tracer MMR.  This is a priori done by the tracer 2 (after inversion)

id=ncdf_open(filesWRF(0))
NCDF_DIMINQ, id, NCDF_DIMID(id, 'west_east'    ), dummy, nx & NCDF_DIMINQ, id, NCDF_DIMID(id, 'south_north'  ), dummy, ny
NCDF_DIMINQ, id, NCDF_DIMID(id, 'bottom_top'   ), dummy, nz & NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), dummy, nt
NCDF_CLOSE, id
id=ncdf_open(filesWRF(nf-1))  ;; for interrupted runs
NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), dummy, ntlast
NCDF_CLOSE, id
nttot = (nf-1)*nt + ntlast
wt  = fltarr(nz,nttot) & wq = fltarr(nz,nttot) & wq_updraft = fltarr(nz,nttot) & wq_downdraft = fltarr(nz,nttot) & wq_env_ude=fltarr(nz,nttot)
q_mean_up = fltarr(nz,nttot) & q_mean_down = fltarr(nz,nttot) & q_mean_env_ude = fltarr(nz,nttot) & q_mean = fltarr(nz,nttot)
tke_les = fltarr(nz,nttot) & ztke = fltarr(nz,nttot) & t = fltarr(nz,nttot)
p = fltarr(nz) & ph = fltarr(nz) & pht = fltarr(nz,nttot) & pt = fltarr(nz,nttot)
xtke = fltarr(nz,nttot) & ytke = fltarr(nz,nttot) & temp_les = fltarr(nz,nttot)
wmax = fltarr(nttot) 
alpha1 = fltarr(nz) & alpha2 = fltarr(nz)
alpha1out = fltarr(nz,nttot) & alpha2out = fltarr(nz,nttot)
zqtrac1 = dblarr(nx,ny,nz) & zqtrac2 = dblarr(nx,ny,nz)
sigmazqtrac1 = fltarr(nz) & sigmazqtrac2 = fltarr(nz)
sigmazminqtrac1 = fltarr(nz) & sigmazminqtrac2 = fltarr(nz)
fm_trac1_les = fltarr(nz,nttot) & fm_trac2_les = fltarr(nz,nttot)
anomalqtrac1 = fltarr(nx,ny,nz) & anomalqtrac2 = fltarr(nx,ny,nz)
e_trac1_les = fltarr(nz,nttot) & e_trac2_les = fltarr(nz,nttot)
dtempdztmp = fltarr(nx,ny,nz)
localtime = lctu + history_interval_s*findgen(nttot)/3700.
w_mean1 = fltarr(nz,nttot) & w_mean2 = fltarr(nz,nttot)
w_mean1_env = fltarr(nz,nttot) & w_mean1_down = fltarr(nz,nttot)
w_mean1_env_ude = fltarr(nz,nttot) & w_mean1_full = fltarr(nz,nttot)
buoyancy1_les = fltarr(nz,nttot) & buoyancy2_les = fltarr(nz,nttot)
e1_term2 = fltarr(nz,nttot) & e1_term3 = fltarr(nz,nttot) 
dtetadttmp = fltarr(nx,ny,nz)
rhomoy1 = fltarr(nz,nttot)
plumeIndex1out = make_array(nx*ny,nz,VALUE=-1.)	& envIndex1out = make_array(nx*ny,nz,VALUE=-1.)
hf1tmp = fltarr(nz,nttot) & hf1tmpenv = fltarr(nz,nttot)
tplume1moy = fltarr(nz,nttot) & tenv1moy = fltarr(nz,nttot) & tenv1moy_ude = fltarr(nz,nttot)
tmoy_full = fltarr(nz,nttot) & tdown1moy = fltarr(nz,nttot)
dteta1moydt_entr = fltarr(nz,nttot) & dteta1moydt_detr = fltarr(nz,nttot)
d1_term1 = fltarr(nz,nttot) & d1_term2 = fltarr(nz,nttot) & d1_term3=fltarr(nz,nttot)
hf1_term1 = fltarr(nz,nttot) & hf1_term2 = fltarr(nz,nttot) & hf1_term3 = fltarr(nz,nttot)
d1_term1_ude = fltarr(nz,nttot) & d1_term2_ude = fltarr(nz,nttot) & d1_term3_ude=fltarr(nz,nttot)
e1_term1_ude = fltarr(nz,nttot) & e1_term2_ude = fltarr(nz,nttot) & e1_term3_ude=fltarr(nz,nttot)
downward_flux1 = fltarr(nz,nttot) & beta1out = fltarr(nz,nttot)
hf1_ude_term1 = fltarr(nz,nttot) & hf1_ude_term2 = fltarr(nz,nttot) & hf1_ude_term3 = fltarr(nz,nttot) & hf1_ude_term4 = fltarr(nz,nttot)
hf1tmpenv_ude = fltarr(nz,nttot) & hf1tmp_down = fltarr(nz,nttot)
dTeta_phys = make_array(nz,nttot) 
exner = fltarr(nz,nttot)
uv_moy = fltarr(nz,nttot)
tsurf = fltarr(nttot)
Gamma_1 = fltarr(nz,nttot) & Gamma_2 = fltarr(nz,nttot) & Gamma_3 = fltarr(nz,nttot)
Gamma_1_tmp = fltarr(nz,nttot) & dgamma1tmp = fltarr(nz,nttot)
ptotprime = fltarr(nx,ny,nz) & anomalptot = fltarr(nx,ny,nz) & dptotprimedztmp  = fltarr(nx,ny,nz)
hfx = fltarr(nttot) & flxrad = fltarr(nttot) & flxgrd = fltarr(nttot) & lwdownz = fltarr(nttot) & swdownz = fltarr(nttot)

;wBin = fltarr(NBINS,nttot) & wBinEnv_ude = fltarr(NBINS,nttot) & wBinUp = fltarr(NBINS,nttot) & wBinDown = fltarr(NBINS,nttot)
l=0

FOR loop  = 0, nf-1 DO BEGIN
  timetime = SYSTIME(1)
  if (loop ne nf-1) then nloop2=nt else nloop2=ntlast
  if (loop ne 0) then loop2_init=0 else loop2_init=1                ;le tout premier pas de temps est l'initialisation, certains champs sont à 0 => bug
        FOR loop2 = loop2_init, nloop2-1 DO BEGIN
	pht(*,l) = TOTAL(TOTAL(getget(filesWRF(loop), 'PHTOT',  count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny) / 1000. / 3.72
        ph = TEMPORARY(ph) + pht(*,l) / (nttot-1)
        ENDFOR
        print, 'computing altitudes, file '+string(loop+1,'(I0)'), SYSTIME(1) - timetime, ' s'
ENDFOR
altitudes_LES  = 1000.*(TEMPORARY(ph)  - hgtu/1000.)  ;; altitude above ground
pht = fltarr(nz,nttot)
ph = fltarr(nz)

FOR loop  = 0, nf-1 DO BEGIN
  timetime = SYSTIME(1)
  if (loop ne nf-1) then nloop2=nt else nloop2=ntlast
  if (loop ne 0) then loop2_init=0 else loop2_init=1                ;le tout premier pas de temps est l'initialisation, certains champs sont à 0 => bug
	FOR loop2 = loop2_init, nloop2-1 DO BEGIN
	
	anomalt = 1. & anomalu = 1. & anomalv = 1. & anomalw = 1.
	; --------------------------------------------------------
	; u' = u and v' = v   (car PAS de background wind !)
	; tke = 0.5 ( <u'^2> + <v'^2> + <w'^2> ) ; u' = u ; v' = v
	; --------------------------------------------------------
	
	tprime = getget(filesWRF(loop), 'T', anomaly=anomalt, count=[0,0,0,1], offset=[0,0,0,loop2])  ;; t' = t - <t>
        t(*,l) = t0 + temporary(anomalt)
	ztke(*,l) = 0.5 * TOTAL(TOTAL(getget(filesWRF(loop), 'W', anomaly=anomalw, count=[0,0,0,1], offset=[0,0,0,loop2])^2,1),1) / float(nx) / float(ny)
	xtke(*,l) = 0.5 * TOTAL(TOTAL(getget(filesWRF(loop), 'U', anomaly=anomalu, count=[0,0,0,1], offset=[0,0,0,loop2])^2,1),1) / float(nx) / float(ny)
	ytke(*,l) = 0.5 * TOTAL(TOTAL(getget(filesWRF(loop), 'V', anomaly=anomalv, count=[0,0,0,1], offset=[0,0,0,loop2])^2,1),1) / float(nx) / float(ny)
	uv_moy(*,l) = TOTAL(TOTAL(sqrt(getget(filesWRF(loop), 'U', count=[0,0,0,1], offset=[0,0,0,loop2])^2 + getget(filesWRF(loop), 'V', count=[0,0,0,1], offset=[0,0,0,loop2])^2),1),1) / float(nx) / float(ny)
	tke_les(*,l) = xtke(*,l) + ytke(*,l) + ztke(*,l)
        wprime = getget(filesWRF(loop), 'W', anomaly=anomalw, count=[0,0,0,1], offset=[0,0,0,loop2])
        pht(*,l) = TOTAL(TOTAL(getget(filesWRF(loop), 'PHTOT',  count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny) / 1000. / 3.72
	pt(*,l)  = TOTAL(TOTAL(getget(filesWRF(loop), 'PTOT' ,  count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny)
        tsurf(l) = TOTAL(TOTAL(getget(filesWRF(loop), 'TSURF' ,  count=[0,0,1], offset=[0,0,loop2]),1),1) / float(nx) / float(ny)
        hfx(l) = TOTAL(TOTAL(getget(filesWRF(loop), 'HFX' ,  count=[0,0,1], offset=[0,0,loop2]),1),1) / float(nx) / float(ny)
        flxrad(l) = TOTAL(TOTAL(getget(filesWRF(loop), 'FLXRAD' ,  count=[0,0,1], offset=[0,0,loop2]),1),1) / float(nx) / float(ny)
        flxgrd(l) = TOTAL(TOTAL(getget(filesWRF(loop), 'FLXGRD' ,  count=[0,0,1], offset=[0,0,loop2]),1),1) / float(nx) / float(ny)
        lwdownz(l) = TOTAL(TOTAL(getget(filesWRF(loop), 'LWDOWNZ' ,  count=[0,0,1], offset=[0,0,loop2]),1),1) / float(nx) / float(ny)
        swdownz(l) = TOTAL(TOTAL(getget(filesWRF(loop), 'SWDOWNZ' ,  count=[0,0,1], offset=[0,0,loop2]),1),1) / float(nx) / float(ny)

        temp_les(*,l) = t(*,l)*(pt(*,l)/p0)^r_cp
	IF (got_pdt eq 'true') then begin
	exner(*,l) = (pt(*,l)/p0)^r_cp
	dTeta_phys(*,l) = (TOTAL(TOTAL(getget(filesWRF(loop), 'PDT', count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny))/exner(*,l)
	ENDIF
	ph = TEMPORARY(ph) + pht(*,l) / (nttot-1)
	p  = TEMPORARY(p ) + pt(*,l) / nttot
        ptotprime = getget(filesWRF(loop), 'PTOT', count=[0,0,0,1], offset=[0,0,0,loop2])
	FOR k=0, nz-1 DO BEGIN
		rhomoy1(*,l) = TOTAL(TOTAL(reform(((ptotprime(*,*,k)/(R*(t(k,l)+tprime(*,*,k))))*(p0/ptotprime(*,*,k))^r_cp)),1),1)/(float(nx)*float(ny))
		anomalptot(*,*,k) = ptotprime(*,*,k) - total(total(ptotprime(*,*,k),1),1)/ float(nx) / float(ny)    ; prime par rapport a la moyenne du niveau (plus continu)
        ENDFOR
	zqtrac1 = getget(filesWRF(loop), s_trac1, count=[0,0,0,1], offset=[0,0,0,loop2])
;        zqtrac2 = getget(filesWRF(loop), s_trac2, count=[0,0,0,1], offset=[0,0,0,loop2])
	FOR i=0,nx-1 DO BEGIN
		FOR j=0, ny-1 DO BEGIN
			dtempdztmp(i,j,*) = deriv(altitudes_LES, tprime(i,j,*) + t(*,l))
	                dptotprimedztmp(i,j,*) = deriv(altitudes_LES, anomalptot(i,j,*))
		ENDFOR
	ENDFOR

        FOR k=0,nz-1 DO BEGIN
	        anomalqtrac1(*,*,k) = zqtrac1(*,*,k) - TOTAL(TOTAL(REFORM(zqtrac1(*,*,k)),1),1)/ float(nx) / float(ny)
;                anomalqtrac2(*,*,k) = zqtrac2(*,*,k) - TOTAL(TOTAL(REFORM(zqtrac2(*,*,k)),1),1)/ float(nx) / float(ny)
                sigmazqtrac1(k) = STDDEV(REFORM(zqtrac1(*,*,k)))
	        IF (k ne 0) THEN BEGIN
         	        subsampledAltitudes = INTERPOL(altitudes_LES(0:k),findgen(k+1),findgen(decimate*k+1)/decimate)
        	       	sigmazminqtrac1(k) = (sigmao/(altitudes_LES(k)-altitudes_LES(0)))*INT_TABULATED(subsampledAltitudes,INTERPOL(sigmazqtrac1(0:k),altitudes_LES(0:k),subsampledAltitudes),/DOUBLE)
	        ENDIF ELSE BEGIN
        	        sigmazminqtrac1(k) = sigmazqtrac1(k)
	        ENDELSE
		plumeIndex1 =  WHERE((anomalqtrac1(*,*,k) GT scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)])) AND ((anomalw(k)+wprime(*,*,k)) GT 0.))
		envIndex1 = WHERE((anomalqtrac1(*,*,k) LE scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)])) OR ((anomalw(k)+wprime(*,*,k)) LE 0.))
;		plumeIndex1 =  WHERE(anomalqtrac1(*,*,k) GT scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)]))
;		envIndex1 = WHERE(anomalqtrac1(*,*,k) LE scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)]))
		IF(plumeIndex1(0) EQ -1) THEN BEGIN
		fm_trac1_les(k,l)=0.
	        e_trac1_les(k,l)=0.
		alpha1out(k,l)=0.
		buoyancy1_les(k,l)=0.
		w_mean1(k,l)=0.
		w_mean1_env(k,l)=0.
		w_mean1_down(k,l)=0.
		w_mean1_full(k,l)=0.
		w_mean1_env_ude(k,l)=0.
		e1_term2(k,l)=0.
		e1_term3(k,l)=0.
                e1_term1_ude(k,l)=0.
                e1_term2_ude(k,l)=0.
                e1_term3_ude(k,l)=0.
		hf1tmp(k,l)=0.
		hf1tmpenv(k,l)=0.
		Gamma_1_tmp(k,l)=0.
		plumeIndex1out(*,k)=-1.
		envIndex1out(*,k)=-1.
		d1_term1(k,l)=0.
		d1_term2(k,l)=0.
		d1_term3(k,l)=0.
                d1_term1_ude(k,l)=0.
                d1_term2_ude(k,l)=0.
                d1_term3_ude(k,l)=0.
		downward_flux1(k,l)=0.
		beta1out(k,l)=0.
		tmoy_full(k,l)=0.
		tdown1moy(k,l)=0.
                wq_updraft(k,l)=0.
                wq_downdraft(k,l)=0.
                q_mean_up(k,l)=0.
                q_mean_down(k,l)=0.
		Gamma_2(k,l)=0.
		Gamma_3(k,l)=0.
		ENDIF ELSE BEGIN
		FOR n=0L,n_elements(plumeIndex1)-1 DO BEGIN
 	        	plumeIndex1out(n,k)=plumeIndex1(n)
		ENDFOR
		FOR n=0L,n_elements(envIndex1)-1 DO BEGIN
                        envIndex1out(n,k)=envIndex1(n)
                ENDFOR
		alpha1(k) = n_elements(plumeIndex1) / float(nx) / float(ny)
                wprimetmp = reform(reform((anomalw(k)+wprime(*,*,k))),[nx*ny,1])
		w_mean1_full(k,l) = mean(wprimetmp)
                w_mean1(k,l) = mean(wprimetmp(plumeIndex1))
                w_mean1_env(k,l) = mean(wprimetmp(envIndex1))
		downdraft_index1 = WHERE((abs(anomalw(k)+wprime(*,*,k)) gt sigmao_ude*STDDEV(wprimetmp(envIndex1))) and (anomalw(k)+wprime(*,*,k) lt 0.))
		
		envIndex1_ude = WHERE(((abs(anomalw(k)+wprime(*,*,k)) le sigmao_ude*STDDEV(wprimetmp(envIndex1))) or (anomalw(k)+wprime(*,*,k) ge 0.)) AND ((anomalqtrac1(*,*,k) LE scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)])) OR ((anomalw(k)+wprime(*,*,k)) LE 0.)))
		IF (envIndex1_ude(0) ne -1) THEN w_mean1_env_ude(k,l) = mean(wprimetmp(envIndex1_ude)) ELSE w_mean1_env_ude(k,l) =0.
		if (downdraft_index1(0) ne -1) then begin
		w_mean1_down(k,l)=mean(wprimetmp(downdraft_index1))
		wprimetmp=0.
		beta1 = n_elements(downdraft_index1) / float(nx) / float(ny) 
		beta1out(k,l)=beta1
		downward_flux1(k,l) = beta1*rhomoy1(k,l)*w_mean1_down(k,l)
		endif else begin
		downward_flux1(k,l)=0.
		beta1out(k,l)=0.
		w_mean1_down(k,l)=0.
		tdown1moy(k,l)=0.
		endelse
		fm_trac1_les(k,l) = alpha1(k)*rhomoy1(k,l)*w_mean1(k,l)
		dtempdztmplin = reform(reform(dtempdztmp(*,*,k)),[nx*ny,1])
                alpha1out(k,l)=alpha1(k)
		tfull=reform(tprime(*,*,k)+t(k,l),[nx*ny,1])
		if (downdraft_index1(0) ne -1) then tdown1moy(k,l)=mean(tfull(downdraft_index1)) 
		tplume1moy(k,l)=mean(tfull(plumeIndex1))
		tenv1moy(k,l)=mean(tfull(envIndex1))
		if (envIndex1_ude(0) ne -1) then tenv1moy_ude(k,l) = mean(tfull(envIndex1_ude)) else tenv1moy_ude(k,l)=0.
		tmoy_full(k,l) = mean(tfull)
		buoyancy1_les(k,l)=grav*(tplume1moy(k,l)/tenv1moy(k,l)-1.)
		e_trac1_les(k,l) = fm_trac1_les(k,l)*TOTAL((1./(tenv1moy(k,l)-tplume1moy(k,l)))*(dtempdztmplin(plumeIndex1)),1)/float(n_elements(plumeIndex1))
		d1_term1(k,l) = fm_trac1_les(k,l)*TOTAL((1./(tenv1moy(k,l)-tplume1moy(k,l)))*(temporary(dtempdztmplin(envIndex1))),1)/float(n_elements(envIndex1))
		if (envIndex1_ude(0) ne -1) then begin
                e1_term1_ude(k,l) = fm_trac1_les(k,l)*TOTAL((1./(tenv1moy_ude(k,l)-tplume1moy(k,l)))*(dtempdztmplin(plumeIndex1)),1)/float(n_elements(plumeIndex1))
                d1_term1_ude(k,l) = fm_trac1_les(k,l)*TOTAL((1./(tenv1moy_ude(k,l)-tplume1moy(k,l)))*(temporary(dtempdztmplin(envIndex1_ude))),1)/float(n_elements(envIndex1_ude))
	     	endif else begin
		e1_term1_ude(k,l) = 0.
		d1_term1_ude(k,l) = 0.
		endelse
		wtmp=reform(wprime(*,*,k)+anomalw(k),[nx*ny,1])
		ttmp=reform(tprime(*,*,k)+t(k,l),[nx*ny,1])
                qtmp=reform(zqtrac1(*,*,k),[nx*ny,1])
                pttmp=reform(ptotprime(*,*,k),[nx*ny,1])
                pt_mean1=mean(pttmp(plumeIndex1))
                pttmp=0.
;                anomalptot(*,*,k) = ptotprime(*,*,k) - temporary(pt_mean1)    ;prime par rapport a la moyenne dans la plume (oscille fortement : probleme de bonne definition verticale de la plume)
                Gamma_1_tmp(k,l) = alpha1out(k,l)*rhomoy1(k,l)*TOTAL((wtmp(plumeIndex1)-w_mean1(k,l))^2,1) / float(n_elements(plumeIndex1))
		anomalptotlin = reform(anomalptot(*,*,k),[nx*ny,1])
                dptotprimedztmplin = reform(dptotprimedztmp(*,*,k),[nx*ny,1])
                Gamma_2(k,l) = -(TOTAL(dptotprimedztmplin(plumeIndex1),1)/float(n_elements(plumeIndex1)))/rhomoy1(k,l)
                Gamma_3(k,l) = -grav*(TOTAL(anomalptotlin(plumeIndex1),1)/float(n_elements(plumeIndex1)))/pt(k,l)
		hf1tmp(k,l) = TOTAL((wtmp(plumeIndex1)-w_mean1(k,l))*(ttmp(plumeIndex1)-tplume1moy(k,l)),1) / float(n_elements(plumeIndex1))
		hf1tmpenv(k,l) = TOTAL((wtmp(envIndex1)-w_mean1_env(k,l))*(ttmp(envIndex1)-tenv1moy(k,l)),1) / float(n_elements(envIndex1))
                q_mean_up(k,l)=mean(qtmp(plumeIndex1))
		wq_updraft(k,l) = TOTAL((wtmp(plumeIndex1)-w_mean1(k,l))*(qtmp(plumeIndex1)-mean(qtmp(plumeIndex1))),1) / float(n_elements(plumeIndex1))

		if (envIndex1_ude(0) ne -1) then begin
                hf1tmpenv_ude(k,l) = TOTAL((wtmp(envIndex1_ude)-w_mean1_env_ude(k,l))*(ttmp(envIndex1_ude)-tenv1moy_ude(k,l)),1) / float(n_elements(envIndex1_ude)) 
                q_mean_env_ude(k,l)=mean(qtmp(envIndex1_ude))
                wq_env_ude(k,l)= TOTAL((wtmp(envIndex1_ude)-mean(wtmp(envIndex1_ude)))*(qtmp(envIndex1_ude)-mean(qtmp(envIndex1_ude))),1) / float(n_elements(envIndex1_ude))
                endif else begin
                hf1tmpenv_ude(k,l) =0.
                wq_env_ude(k,l)=0.
                endelse
		if (downdraft_index1(0) ne -1) then begin
                hf1tmp_down(k,l) = TOTAL((wtmp(downdraft_index1)-w_mean1_down(k,l))*(ttmp(downdraft_index1)-tdown1moy(k,l)),1) / float(n_elements(downdraft_index1))
                q_mean_down(k,l)=mean(qtmp(downdraft_index1))
                wq_downdraft(k,l) = TOTAL((wtmp(downdraft_index1)-mean(wtmp(downdraft_index1)))*(qtmp(downdraft_index1)-mean(qtmp(downdraft_index1))),1) / float(n_elements(downdraft_index1))
                endif else begin
                hf1tmp_down(k,l)=0.
		wq_downdraft(k,l)=0.
                endelse
                q_mean(k,l)=mean(qtmp)
		IF((n_elements(plumeIndex1) + n_elements(envIndex1)) ne float(nx*ny)) then print, 'WARNING : INDEX PROBLEM : plume / env : ', n_elements(plumeIndex1), n_elements(envIndex1)
;		IF((n_elements(plumeIndex1) + n_elements(envIndex1_ude) + n_elements(downdraft_index1)) ne float(nx*ny)) then print, 'WARNING : INDEX PROBLEM : plume / env / downdraft : ', n_elements(plumeIndex1), n_elements(envIndex1_ude), n_elements(downdraft_index1)
		ENDELSE
	ENDFOR

;        FOR i=0, nx-1 DO BEGIN
;            FOR j=0, ny-1 DO BEGIN
;                dptotprimedztmp(i,j,*) = deriv(altitudes_LES, anomalptot(i,j,*))
;            ENDFOR
;        ENDFOR
;
;        FOR k=0, nz-1 DO BEGIN
;		IF(plumeIndex1out(0,k) eq -1) THEN BEGIN
;			Gamma_2(k,l)=0.
;                ENDIF ELSE BEGIN
;			wpi=where(plumeIndex1out ne -1)
;			plumeIndex1=reform(plumeIndex1out(temporary(wpi),k))
;		        dptotprimedztmplin = reform(dptotprimedztmp(*,*,k),[nx*ny,1])
;		        Gamma_2(k,l) = -(TOTAL(dptotprimedztmplin(plumeIndex1),1)/float(n_elements(plumeIndex1)))/rhomoy1(k,l)
;		ENDELSE
;	ENDFOR

        dgamma1tmp(*,l) = deriv(altitudes_LES,Gamma_1_tmp(*,l))
        FOR k=0, nz-1 DO BEGIN
	        IF (alpha1out(k,l) ne 0.) THEN BEGIN
		        Gamma_1(k,l) = -(1./(alpha1out(k,l)*rhomoy1(k,l)))*dgamma1tmp(k,l)
	        ENDIF ELSE BEGIN
		        Gamma_1(k,l)=0.
	 	ENDELSE
        ENDFOR
        drhoahfdztmp = deriv(altitudes_LES,rhomoy1(*,l)*alpha1out(*,l)*hf1tmp(*,l))
        drhoahfdztmpDetr = deriv(altitudes_LES,rhomoy1(*,l)*(1.-alpha1out(*,l))*hf1tmpenv(*,l))
        drhoahfdztmpDetr_ude = deriv(altitudes_LES,rhomoy1(*,l)*(1.-alpha1out(*,l)-beta1out(*,l))*hf1tmpenv_ude(*,l))

;	wBin(*,l)=HISTOGRAM(wtmp,nbins=NBINS)
;       wBinEnv_ude(*,l)=HISTOGRAM(wtmp(envIndex1_ude),nbins=NBINS)
;	wBinUp(*,l)=HISTOGRAM(wtmp(plumeIndex1),nbins=NBINS)
;	wBinDown(*,l)=HISTOGRAM(wtmp(downdraft_index1),nbins=NBINS)
	wtmp=0.
	ttmp=0.
	
	FOR k=0,nz-1 DO BEGIN
		IF(plumeIndex1out(0,k) eq -1) THEN BEGIN
		e1_term2(k,l)=0.
		e1_term2_ude(k,l)=0.
		ENDIF ELSE BEGIN
		e1_term2(k,l)=drhoahfdztmp(k)/(tenv1moy(k,l)-tplume1moy(k,l))
		e1_term2_ude(k,l)=drhoahfdztmp(k)/(tenv1moy_ude(k,l)-tplume1moy(k,l))
		ENDELSE

		IF(envIndex1out(0,k) eq -1) THEN BEGIN
                d1_term2(k,l)=0.
		d1_term2_ude(k,l)=0.
                ENDIF ELSE BEGIN
		d1_term2(k,l)=-drhoahfdztmpDetr(k)/(tenv1moy(k,l)-tplume1moy(k,l))
                d1_term2_ude(k,l)=-drhoahfdztmpDetr_ude(k)/(tenv1moy_ude(k,l)-tplume1moy(k,l))
		ENDELSE
	ENDFOR


	tfull1=0.

        zqtrac1=0.
;        zqtrac2 = getget(filesWRF(loop), s_trac2, count=[0,0,0,1], offset=[0,0,0,loop2])
;        FOR k=0,nz-1 DO BEGIN
;		anomalqtrac2(*,*,k) = zqtrac2(*,*,k) - TOTAL(TOTAL(REFORM(zqtrac2(*,*,k)),1),1)/ float(nx) / float(ny)
;                sigmazqtrac2(k) = STDDEV(zqtrac2(*,*,k))
;	        IF (k ne 0) THEN BEGIN
;	                subsampledAltitudes = INTERPOL(altitudes_LES(0:k),findgen(k+1),findgen(decimate*k+1)/decimate)
;	                sigmazminqtrac2(k) = (sigmao/(altitudes_LES(k)-altitudes_LES(0)))*INT_TABULATED(subsampledAltitudes,INTERPOL(sigmazqtrac2(0:k),altitudes_LES(0:k),subsampledAltitudes),/DOUBLE)
;	        ENDIF ELSE BEGIN
;        	        sigmazminqtrac2(k) = sigmazqtrac2(k)
;	        ENDELSE
;                plumeIndex2 = WHERE((anomalqtrac2(*,*,k) GT scale*MAX([sigmazqtrac2(k),sigmazminqtrac2(k)])) AND ((wprime(*,*,k)+anomalw(k)) GT 0.))
;                envIndex2 = WHERE((anomalqtrac2(*,*,k) LE scale*MAX([sigmazqtrac2(k),sigmazminqtrac2(k)])) OR ((wprime(*,*,k)+anomalw(k)) LE 0.))
;                IF(plumeIndex2(0) EQ -1) THEN BEGIN
;		 fm_trac2_les(k,l)=0.
;		 e_trac2_les(k,l)=0.
;                alpha2out(k,l)=0.
;                buoyancy2_les(k,l)=0.
;		w_mean2(k,l)=0.
;                ENDIF ELSE BEGIN
;		alpha2(k) = n_elements(plumeIndex2) / float(nx) / float(ny)
;                wprimetmp = reform(reform((anomalw(k)+wprime(*,*,k))),[nx*ny,1])
;                w_mean2(k,l) = mean(wprimetmp(plumeIndex2))
;		wprimetmp=0.
;                fm_trac2_les(k,l) = alpha2(k)*rhomoy1(k,l)*w_mean2(k,l)
;                tprimetmp = reform(reform(-tprime(*,*,k)),[nx*ny,1])
;                dtempdztmplin = reform(reform(dtempdztmp(*,*,k)),[nx*ny,1])
;                e_trac2_les(k,l) = TOTAL((1./(temporary(tprimetmp(plumeIndex2))))*(temporary(dtempdztmplin(plumeIndex2))),1)/n_elements(plumeIndex2)
;                alpha2out(k,l)=alpha2(k)
;                tfull=reform(tprime(*,*,k)+t(k,l),[nx*ny,1])
;                tplume2moy=mean(tfull(plumeIndex2))
;                tenv2moy=mean(tfull(envIndex2))
;                buoyancy2_les(k,l)=grav*(tplume2moy/tenv2moy-1.)
;		ENDELSE
;        ENDFOR
;        zqtrac2=0.

        wt(*,l) = TOTAL(TOTAL(TEMPORARY(tprime)*wprime,1),1) / float(nx) / float(ny)
	wq(*,l) = TOTAL(TOTAL(TEMPORARY(wprime)*anomalqtrac1,1),1) / float(nx) / float(ny)

	wmax(l) = max(w_mean1(*,l))
        l=l+1
	ENDFOR
   	print, 'file '+string(loop+1,'(I0)'), SYSTIME(1) - timetime, ' s'

ENDFOR

hf1_term1 = hf1tmp*alpha1out
hf1_term2 = temporary(hf1tmpenv)*(1.-alpha1out)
hf1_term3 = alpha1out*(1.-alpha1out)*(w_mean1 - w_mean1_env)*(tplume1moy - tenv1moy)

hf1_ude_term1 = temporary(hf1tmp)*alpha1out
hf1_ude_term2 = temporary(hf1tmp_down)*beta1out
hf1_ude_term3 = temporary(hf1tmpenv_ude)*(1.-(alpha1out+beta1out))
hf1_ude_term4 = alpha1out*(w_mean1 - w_mean1_full)*(tplume1moy - tmoy_full) + beta1out*(w_mean1_down - w_mean1_full)*(tdown1moy - tmoy_full) + (1.- (alpha1out+beta1out))*(w_mean1_env_ude - w_mean1_full)*(tenv1moy_ude - tmoy_full)

FOR k=0, nz-1 DO BEGIN
;       dteta1moydt_entr(k,*) = deriv(localtime,tplume1moy(k,*))/3700. - dTeta_phys(k,*)
;       dteta1moydt_detr(k,*) = deriv(localtime,tplume1moy(k,*))/3700. + dTeta_phys(k,*)
	dteta1moydt_entr(k,*) = deriv(localtime,tplume1moy(k,*))/3700. - dTeta_phys(k,*)
	dteta1moydt_detr(k,*) =  dTeta_phys(k,*) - deriv(localtime,tplume1moy(k,*))/3700. 
ENDFOR

FOR k=0, nz-1 DO BEGIN
	FOR l=0, nttot-1 DO BEGIN
		IF (tenv1moy(k,l) ne tplume1moy(k,l)) THEN e1_term3(k,l) = rhomoy1(k,l)*alpha1out(k,l)*dteta1moydt_entr(k,l)/(tenv1moy(k,l)-tplume1moy(k,l)) ELSE e1_term3(k,l)=0.
                IF (tenv1moy_ude(k,l) ne tplume1moy(k,l)) THEN e1_term3_ude(k,l) = rhomoy1(k,l)*alpha1out(k,l)*dteta1moydt_entr(k,l)/(tenv1moy_ude(k,l)-tplume1moy(k,l)) ELSE e1_term3_ude(k,l)=0
;		IF (tenv1moy(k,l) ne tplume1moy(k,l)) THEN d1_term3(k,l) = rhomoy1(k,l)*(1.-alpha1out(k,l))*dteta1moydt_detr(k,l)/(tenv1moy(k,l)-tplume1moy(k,l)) ELSE d1_term3(k,l)=0.
		IF (tenv1moy(k,l) ne tplume1moy(k,l)) THEN d1_term3(k,l) = rhomoy1(k,l)*(1.-alpha1out(k,l))*dteta1moydt_detr(k,l)/(tenv1moy(k,l)-tplume1moy(k,l)) ELSE d1_term3(k,l)=0.
                IF (tenv1moy_ude(k,l) ne tplume1moy(k,l)) THEN d1_term3_ude(k,l) = rhomoy1(k,l)*(1.-alpha1out(k,l)-beta1out(k,l))*dteta1moydt_detr(k,l)/(tenv1moy_ude(k,l)-tplume1moy(k,l)) ELSE d1_term3_ude(k,l)=0.
	ENDFOR
ENDFOR


ht = TEMPORARY(pht) - hgtu/1000.
save, hfx, flxrad, flxgrd, lwdownz, swdownz, q_mean_up, q_mean_down, q_mean_env_ude, q_mean, Gamma_1, Gamma_2, Gamma_3, w_mean1_env, tsurf, wq, wq_updraft, wq_downdraft, wq_env_ude, d1_term1_ude, d1_term2_ude, d1_term3_ude, e1_term1_ude, e1_term2_ude, e1_term3_ude, tplume1moy, tdown1moy, w_mean1_full, tmoy_full, tenv1moy_ude, w_mean1_env_ude, uv_moy, hf1_ude_term1, hf1_ude_term2, hf1_ude_term3, hf1_ude_term4, w_mean1_down, downward_flux1, beta1out, hf1_term1, hf1_term2, hf1_term3, d1_term1, d1_term2, d1_term3, e1_term2, e1_term3, buoyancy1_les, buoyancy2_les, w_mean1, w_mean2, nx, ny, alpha1out, alpha2out, e_trac1_les, e_trac2_les, tke_les, ztke, altitudes_LES, ht, t, p, pt, localtime, xtke, ytke, wt, temp_les, wmax, fm_trac1_les, fm_trac2_les,filename=les_path+'/'+datname

nz = n_elements(altitudes_LES)

ENDIF ELSE BEGIN

print, 'OK, file is here'
restore, filename=les_path+'/'+datname
nz = n_elements(altitudes_LES)
nttot = n_elements(tmoy_full(0,*))

OPENR, 23, les_path+'/input_more' & READF, 23, hgtu, tsurfu & CLOSE, 23

ENDELSE

tenv1moy = tplume1moy/((buoyancy1_les/grav)+1.)

taverage = string((localtime(nstot)-localtime(1))*3700./60.)
print, ''
print, ' -- Loading testphys1d data -- '

file1=gcm_path+'/diagfi.nc'
file2=gcm_convadj_path+'/diagfi.nc'
file3='/san0/acolmd/SIMUS/GCM3D_TestBed/diagfi.nc'

getcdf, file=file1, charvar='q2', invar=tke_gcm
getcdf, file=file1, charvar='aps', invar=aps
getcdf, file=file1, charvar='bps', invar=bps
getcdf, file=file1, charvar='co2col', invar=co2_col
;getcdf, file=file1, charvar='arcol', invar=ar_col
getcdf, file=file1, charvar='tkecol', invar=tke_col
;getcdf, file=file1, charvar='ar', invar=ar
getcdf, file=file1, charvar='heatFlux_up', invar=heatFlux_up
getcdf, file=file1, charvar='heatFlux_down', invar=heatFlux_down
getcdf, file=file1, charvar='pplay', invar=pplay
getcdf, file=file1, charvar='pplev', invar=pplev
getcdf, file=file1, charvar='temp', invar=temp_gcm
getcdf, file=file1, charvar='zw2', invar=zw2_lev
getcdf, file=file1, charvar='fm_therm', invar=fm_therm_gcm_lev
getcdf, file=file1, charvar='entr_therm', invar=zdz_entr_therm_gcm
getcdf, file=file1, charvar='detr_therm', invar=zdz_detr_therm_gcm
getcdf, file=file1, charvar='fraca', invar=alpha_gcm_lev
getcdf, file=file1, charvar='buoyancyOut', invar=buoyancy_gcm
getcdf, file=file1, charvar='buoyancyEst', invar=buoyancy_est_gcm
getcdf, file=file1, charvar='zkh', invar=zkh
getcdf, file=file1, charvar='zh', invar=zh
getcdf, file=file1, charvar='tsurf', invar=tsurf_gcm

if (overplot_convadj eq 'true') then begin
getcdf, file=file2, charvar='temp', invar=temp_gcm_convadj
getcdf, file=file2, charvar='pplay', invar=pplay_convadj
endif
if (plot_3d eq 'true') then begin
getcdf, file=file3, charvar='temp', invar=temp_gcm_3d
getcdf, file=file3, charvar='pplay', invar=pplay_3d
getcdf, file=file3, charvar='latitude', invar=latitude_3d
getcdf, file=file3, charvar='longitude', invar=longitude_3d
nWEmx_3d = n_elements(reform(temp_gcm_3d(*,0,0,0)))
nNSmx_3d = n_elements(reform(temp_gcm_3d(0,*,0,0)))
nZmx_3d = n_elements(reform(temp_gcm_3d(0,0,*,0)))
nTmx_3d = n_elements(reform(temp_gcm_3d(0,0,0,*)))
ndays_3d = 1.
lctu_gcm_3d = 0.
history_interval_s_gcm_3d = ndays_3d*88800./float(nTmx_3d)  ; Timestep interval of gcm 1d simu in sec
localtime_lon0 = lctu_gcm_3d + history_interval_s_gcm_3d*findgen(nTmx_3d)/3700.
Xc = 205.
Yc = 21.8
plot_index_x = (Xc-longitude_3d(0))/(longitude_3d(1)-longitude_3d(0))
plot_index_y = (Yc-latitude_3d(0))/(latitude_3d(1)-latitude_3d(0))
localtime_true = localtime_lon0 -(12./180.)*Xc
endif


nTmx = n_elements(reform(temp_gcm(0,*)))
if (overplot_convadj eq 'true') then begin
nTmx_convadj = n_elements(reform(temp_gcm_convadj(0,*)))
endif else begin
nTmx_convadj =  10000.
endelse
ndays = 1.
print, ''
print, 'WARNING ----------------------- '
print, 'CONFIGURATION : '+string(ndays,format='(I0)')+' days simulation'
print, ''
history_interval_s_gcm = ndays*88800./float(nTmx)  ; Timestep interval of gcm 1d simu in sec
history_interval_s_gcm_convadj = ndays*88800./float(nTmx_convadj)
localtime_gcm = lctu_gcm + history_interval_s_gcm*findgen(nTmx)/3700.
localtime_gcm_convadj = lctu_gcm + history_interval_s_gcm_convadj*findgen(nTmx_convadj)/3700.
; **********************************
; ******** PLOTS ******************
if (f_offset eq 'true') then begin
offset_localtime = 3.108100
endif else begin
offset_localtime = 0.
endelse
localtime=localtime+offset_localtime
print, '****************************************************************************************************'
print, 'local time LES'
print, localtime
print, '****************************************************************************************************'
print, 'local time GCM'
print, localtime_gcm
print, '****************************************************************************************************'

time_offset = (ndays-1.)*24.

lt_plot_ini = 8.
lt_plotindex_les_ini = where(localtime eq lt_plot_ini)
lt_plotindex_gcm_ini = where(localtime_gcm eq (lt_plot_ini+time_offset))

lt_plot0 = 10.
lt_plotindex_les0 = where(localtime eq lt_plot0)
lt_plotindex_gcm0 = where(localtime_gcm eq (lt_plot0+time_offset))
lt_plotindex_gcm_convadj0 = where(localtime_gcm_convadj eq (lt_plot0+time_offset))

;lt_plot = 8.25
lt_plot = 12.
lt_plotindex_les = where((localtime lt lt_plot+0.01) and (localtime gt lt_plot-0.01))
lt_plotindex_gcm = where(localtime_gcm eq (lt_plot+time_offset))
lt_plotindex_gcm_convadj = where(localtime_gcm_convadj eq (lt_plot+time_offset))
print, 'lt plotindex les 12h'
print, lt_plotindex_les

lt_plot2 = 14.
lt_plotindex_les2 = where(localtime eq lt_plot2)
lt_plotindex_gcm2 = where(localtime_gcm eq (lt_plot2+time_offset))
lt_plotindex_gcm_convadj2 = where(localtime_gcm_convadj eq (lt_plot2+time_offset))

lt_plot3 = 16.
lt_plotindex_les3 = where(localtime eq lt_plot3)
lt_plotindex_gcm3 = where(localtime_gcm eq (lt_plot3+time_offset))
lt_plotindex_gcm_convadj3 = where(localtime_gcm_convadj eq (lt_plot3+time_offset))

lt_plot4 = 18.
lt_plotindex_les4 = where(localtime eq lt_plot4)
lt_plotindex_gcm4 = where(localtime_gcm eq (lt_plot4+time_offset))
lt_plotindex_gcm_convadj4 = where(localtime_gcm_convadj eq (lt_plot4+time_offset))
;--------------------------------------------------------------------------------
;---------------------------------------------------------------------------------

nTmx_les=n_elements(reform(wt(0,*)))
nZmx=n_elements(aps)          ; number of vertical levels
H_low=9650.                   ; scale height at low altitudes
H_high=15000.                 ; scale height at high altitudes
trans_window=10.              ; # of levels over which H(:) goes from H_low to H_high
lev_trans=32.+trans_window/2. ; level at which H(lev_trans)=(H_low+H_high)/2
P_ref=p0                    ; reference surface pressure used to build zsurface -610 Pa-
Hgcm = make_array(nZmx)
altitudes_GCM = make_array(nZmx)
; Build scale heights
;FOR k=0,nZmx-1 DO BEGIN
;	 Hgcm(k)=H_low+(H_high-H_low)*0.5*(1.0+tanh(6.*(k-lev_trans)/trans_window))
;ENDFOR

FOR k=0,nZmx-1 DO BEGIN
        Hgcm(k)=R*temp_gcm(k,lt_plotindex_gcm)/grav
ENDFOR
print, 'Hgcm'
print, Hgcm
; Compute altitudes_GCM
FOR k=0,nZmx-1 DO BEGIN
	altitudes_GCM(k)=-Hgcm(k)*alog(pplay(k,lt_plotindex_gcm)/pGround)
ENDFOR
Hgcm=0.

teta_gcm = temp_gcm * (p0/pplay)^r_cp
if (overplot_convadj eq 'true') then begin
teta_gcm_convadj = temp_gcm_convadj * (p0/pplay_convadj)^r_cp
endif

OPENR, 1, gcm_path+'/profile'
data=FLTARR(nZmx+1)
READF, 1, data
temp_gcm_0_ground = data(0)
temp_gcm_0 = data(1:nZmx-1)
data = 0.
CLOSE, 1

teta_gcm_0 = temp_gcm_0 * (p0/pplay)^r_cp
approx_zdz_gcm = make_array(nZmx)
approx_zdz_gcm(0)=altitudes_GCM(1)
FOR k=1, nZmx-2 DO BEGIN
       approx_zdz_gcm(k) = altitudes_GCM(k+1) - altitudes_GCM(k)
ENDFOR
approx_zdz_gcm(nZmx-1)=approx_zdz_gcm(nZmx-2)

print, 'approx zdz gcm'
print, approx_zdz_gcm

print, '****************************************************************************************************'
print, 'altitudes LES based on phtot : inter-levels'
print, altitudes_LES
print, '****************************************************************************************************'
print, 'altitudes GCM based on pplay : inter-levels'
print, altitudes_GCM
print, '****************************************************************************************************'

; Compute tracer deviation :

co2_col = co2_col/co2_col(0)
;ar_col = ar_col/ar_col(0)
tke_col = tke_col+1.

; Compute <teta> les

teta_les = temporary(t) 


; ========================================================================
; ========================================================================

IF (visualization_mode eq 'true') THEN BEGIN

print,' *****************************************-----------------------------------'
print,' ************ PLUME **********************-----------------------------------'
print,' *****************************************-----------------------------------'

; We are evaluating the first time-step element of the file number 'loop-1' :
; file 1 starts at 8h (loop =0, loop2 =0)
; file 6 starts at 13h  (roughly)  (loop =5, loops2=0)
; file 12 starts at 18h (roughly) (loop = 11,loop2 = 0)

loop=uint(loop_special)-1
;loop2=34
loop2=10
domain='d01'
filesWRF = FindFile(les_path+'/wrfout_'+domain+'_????-??-??_??:??:??')
anomalw=1.

zqtrac1 = dblarr(nx,ny,nz) & zqtrac2 = dblarr(nx,ny,nz)
sigmazqtrac1 = fltarr(nz) & sigmazqtrac2 = fltarr(nz)
sigmazminqtrac1 = fltarr(nz) & sigmazminqtrac2 = fltarr(nz)
anomalqtrac1 = fltarr(nx,ny,nz) & anomalqtrac2 = fltarr(nx,ny,nz)
zqtrac1 = getget(filesWRF(loop), 'qtrac1', count=[0,0,0,1], offset=[0,0,0,loop2])
zqtrac2 = getget(filesWRF(loop), 'qtrac2', count=[0,0,0,1], offset=[0,0,0,loop2])
wprime = getget(filesWRF(loop), 'W', anomaly = anomalw, count=[0,0,0,1], offset=[0,0,0,loop2])
supermask1 = fltarr(nx,ny,nz)
supermask2 = fltarr(nx,ny,nz)
k_out_histo = 8
k_out_hist = [1,10,25,50,70,85]
nbtest=n_elements(k_out_hist)
b=0
FOR k=0,nz-1 DO BEGIN
        mask1=fltarr(nx*ny)
        mask2=fltarr(nx*ny)
        anomalqtrac1(*,*,k) = zqtrac1(*,*,k) - TOTAL(TOTAL(REFORM(zqtrac1(*,*,k)),1),1)/ float(nx) / float(ny)
        sigmazqtrac1(k) = STDDEV(REFORM(zqtrac1(*,*,k)))
	IF (k ne 0) THEN BEGIN
	        subsampledAltitudes = INTERPOL(altitudes_LES(0:k),findgen(k+1),findgen(decimate*k+1)/decimate)
		sigmazminqtrac1(k) = (sigmao/(altitudes_LES(k)-altitudes_LES(0)))*INT_TABULATED(subsampledAltitudes,INTERPOL(sigmazqtrac1(0:k),altitudes_LES(0:k),subsampledAltitudes),/DOUBLE)
	ENDIF ELSE BEGIN
		sigmazminqtrac1(k) = sigmazqtrac1(k)
	ENDELSE
	print, sigmazqtrac1(k),sigmazminqtrac1(k)
        plumeIndex1 =  WHERE((anomalqtrac1(*,*,k) GT scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)])) AND ((anomalw(k)+wprime(*,*,k)) GT 0.))
;        anomalqtrac2(*,*,k) = zqtrac2(*,*,k) - TOTAL(TOTAL(REFORM(zqtrac2(*,*,k)),1),1)/ float(nx) / float(ny)
;        sigmazqtrac2(k) = STDDEV(REFORM(zqtrac2(*,*,k)))
        IF (k ne 0) THEN BEGIN
                subsampledAltitudes = INTERPOL(altitudes_LES(0:k),findgen(k+1),findgen(decimate*k+1)/decimate)
;                sigmazminqtrac2(k) = (sigmao/(altitudes_LES(k)-altitudes_LES(0)))*INT_TABULATED(subsampledAltitudes,INTERPOL(sigmazqtrac2(0:k),altitudes_LES(0:k),subsampledAltitudes),/DOUBLE)
	ENDIF ELSE BEGIN
;                sigmazminqtrac2(k) = sigmazqtrac2(k)
        ENDELSE
;        plumeIndex2 =  WHERE((anomalqtrac2(*,*,k) GT scale*MAX([sigmazqtrac2(k),sigmazminqtrac2(k)])) AND ((anomalw(k)+wprime(*,*,k)) GT 0.))
        IF(plumeIndex1(0) NE -1 ) THEN BEGIN
        mask1(plumeIndex1) = 1.
        ENDIF ELSE BEGIN
        mask1(*)=0.
        ENDELSE
        IF(plumeIndex2(0) NE -1 ) THEN BEGIN
        mask2(plumeIndex2) = 1.
        ENDIF ELSE BEGIN
        mask2(*)=0.
        ENDELSE
        mask1 = reform(mask1,[nx,ny])
        supermask1(*,*,k) = mask1(*,*)
        mask2 = reform(mask2,[nx,ny])
        supermask2(*,*,k) = mask2(*,*)
;	IF (k eq k_out_histo) THEN BEGIN
;		plume1_out = plumeIndex1
;		plume2_out = plumeIndex2
;	ENDIF	

	IF (k eq k_out_hist(0)) THEN BEGIN
		c1=plumeIndex1
		c2=plumeIndex2
	ENDIF
        IF (k eq k_out_hist(1)) THEN BEGIN
                d1=plumeIndex1
                d2=plumeIndex2
        ENDIF
        IF (k eq k_out_hist(2)) THEN BEGIN
                e1=plumeIndex1
                e2=plumeIndex2
        ENDIF
        IF (k eq k_out_hist(3)) THEN BEGIN
                f1=plumeIndex1
                f2=plumeIndex2
        ENDIF
        IF (k eq k_out_hist(4)) THEN BEGIN
                g1=plumeIndex1
                g2=plumeIndex2
        ENDIF
        IF (k eq k_out_hist(5)) THEN BEGIN
                h1=plumeIndex1
                h2=plumeIndex2
        ENDIF
ENDFOR

IF (Histo eq 'false') THEN BEGIN
IVOLUME, supermask1
IVOLUME, supermask2
ENDIF ELSE BEGIN

; -------------------------------------------------------------------------------
; THIS IS THE ULTRA-GORE UN-ESTHETIC UGLY-AS-HELL LOOP FOR CONCENTRATION PLOTTING
; but well, this is just because idl cannot handle arrays as well as I would like
; -------------------------------------------------------------------------------
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Distrib.ps'
PSOPEN, THICK=100, CHARSIZE=60, FILE = filename, FONT = 5, TFONT = 5,XPLOTS=2,YPLOTS=3,MARGIN=2000 

FOR n=0, nbtest-1 DO BEGIN

CASE n OF
	0:BEGIN
		plume1_out=c1 & plume2_out=c2
	END
        1:BEGIN
                plume1_out=d1 & plume2_out=d2
        END
        2:BEGIN
                plume1_out=e1 & plume2_out=e2
        END
        3:BEGIN
                plume1_out=f1 & plume2_out=f2
        END
        4:BEGIN
                plume1_out=g1 & plume2_out=g2
        END
        5:BEGIN
                plume1_out=h1 & plume2_out=h2
        END
ENDCASE
 
ToBin1 = reform(zqtrac1(*,*,k_out_hist(n)),[nx*ny,1])
;ToBin2 = reform(zqtrac2(*,*,k_out_hist(n)),[nx*ny,1])
svmin1=min(ToBin1) & svmin2=min(ToBin2)
svmax1=max(ToBin1) & svmax2=max(ToBin2)
NBINS=100
ds1=(svmax1-svmin1+1.)/(NBINS-1) & ds2=(svmax2-svmin2+1.)/(NBINS-1)
Xaxis1 = svmin1+((svmax1-svmin1)/(NBINS-1))*indgen(NBINS)
Xaxis2 = svmin2+((svmax2-svmin2)/(NBINS-1))*indgen(NBINS)
Bin1=HISTOGRAM(ToBin1,nbins=NBINS)
Bin2=INTERPOL(HISTOGRAM(ToBin2,nbins=NBINS),Xaxis2,Xaxis1,/SPLINE)

what_I_plot = [[Bin1],[Bin2]]
labels=['LES tracer 1 conc. distrib.','LES tracer 2 conc. distrib.']

title_user = TestCase+SubCase+' LES tracer 1&2 concentration distribution at '+string(altitudes_LES(k_out_hist(n)))+'m AGL'
IF (n lt 3) THEN BEGIN
	POS, XPOS=1, YPOS=uint(n+1)
ENDIF ELSE BEGIN
	POS, XPOS=2, YPOS=uint(n+1)-3
ENDELSE
MAP
CS, SCALE=28
GSET, XMIN=0, XMAX=20, YMIN=0, YMAX=300, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=Xaxis1, Y=what_I_plot, /LEGEND, LEGPOS=9, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 4, XTITLE='Trac concentration (kg/kg)', YSTEP=50, YTITLE='Counts',NDECS=1

ToBin1 = reform(zqtrac1(*,*,k_out_hist(n)),[nx*ny,1])
ToBin2 = reform(zqtrac2(*,*,k_out_hist(n)),[nx*ny,1])
ToBin1 = ToBin1(plume1_out)
ToBin2 = ToBin2(plume2_out)
svmin1=min(ToBin1) & svmin2=min(ToBin2)
svmax1=max(ToBin1) & svmax2=max(ToBin2)
NBINS=50
ds1=(svmax1-svmin1+1.)/(NBINS-1) & ds2=(svmax2-svmin2+1.)/(NBINS-1)
Xaxis1 = svmin1+((svmax1-svmin1)/(NBINS-1))*indgen(NBINS)
Xaxis2 = svmin2+((svmax2-svmin2)/(NBINS-1))*indgen(NBINS)
Bin1=HISTOGRAM(ToBin1,nbins=NBINS)
Bin2=HISTOGRAM(ToBin2,nbins=NBINS)
oplot,  Xaxis1, Bin1, psym=4
oplot,  Xaxis2, Bin2, psym=5

ENDFOR

PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

ENDELSE

ENDIF ELSE BEGIN

; =========================================================================
; =========================================================================

print, '........ ALPHA'

alpha_interlay_gcm = make_array(nZmx)
FOR k=0, nZmx-2 DO BEGIN
        alpha_interlay_gcm(k) = (alpha_gcm_lev(k,lt_plotindex_gcm)+alpha_gcm_lev(k+1,lt_plotindex_gcm))/2.
ENDFOR
alpha_interlay_gcm(nZmx-1)=0.

smoothed_alpha1_les = make_array(nz)
smoothed_alpha2_les = make_array(nz)
smoothed_beta1_les = make_array(nz)
FOR t=-ns,ns DO BEGIN
        smoothed_alpha1_les = smoothed_alpha1_les + REFORM(alpha1out(*,lt_plotindex_les+t))
        smoothed_alpha2_les = smoothed_alpha2_les + REFORM(alpha2out(*,lt_plotindex_les+t))
        smoothed_beta1_les = smoothed_beta1_les + REFORM(beta1out(*,lt_plotindex_les+t))
ENDFOR

smoothed_alpha1_les = smoothed_alpha1_les/nstot
smoothed_alpha2_les = smoothed_alpha2_les/nstot
smoothed_beta1_les = smoothed_beta1_les/nstot

; =========================================================================

; *** Temperatures ***

print, '........ TEMPERATURES'

xmin = 190
xmax = 250
if (TestCase eq 'Case_Z') then begin
xmin = 170
xmax = 250
endif


what_I_plot = [[reform(temp_gcm(*,lt_plotindex_gcm_ini))],[reform(temp_gcm(*,lt_plotindex_gcm0))],[reform(temp_gcm(*,lt_plotindex_gcm))],[reform(temp_gcm(*,lt_plotindex_gcm2))],[reform(temp_gcm(*,lt_plotindex_gcm3))],[reform(temp_gcm(*,lt_plotindex_gcm4))]]
labels=['TH temp 1d, lt='+string(lt_plot_ini),'TH temp 1d, lt='+string(lt_plot0),'TH temp 1d, lt='+string(lt_plot),'TH temp 1d, lt='+string(lt_plot2),'TH temp 1d, lt='+string(lt_plot3),'TH temp 1d, lt='+string(lt_plot4)]
title_user = TestCase+SubCase+LayerCase+' Temperatures Comparison'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Temperature.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=xmin, XMAX=xmax, YMIN=0, YMAX=12, TITLE=title_user
cols=INDGEN(6)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 5, XTITLE='Temperature (K)', YSTEP=1, YTITLE='Altitude (km)',NDECS=1

oplot, temp_les(*,lt_plotindex_les_ini), altitudes_LES/1000., psym=4
oplot, temp_les(*,lt_plotindex_les0), altitudes_LES/1000., psym=4
oplot, temp_les(*,lt_plotindex_les), altitudes_LES/1000., psym=4
oplot, temp_les(*,lt_plotindex_les2), altitudes_LES/1000., psym=4
oplot, temp_les(*,lt_plotindex_les3), altitudes_LES/1000., psym=4
oplot, temp_les(*,lt_plotindex_les4), altitudes_LES/1000., psym=4

if(overplot_convadj EQ 'true') then begin
oplot, temp_gcm_convadj(*,lt_plotindex_gcm_convadj0), altitudes_GCM/1000., thick=0.1,color=8,linestyle=3
oplot, temp_gcm_convadj(*,lt_plotindex_gcm_convadj), altitudes_GCM/1000., thick=0.1,color=8,linestyle=3
oplot, temp_gcm_convadj(*,lt_plotindex_gcm_convadj2), altitudes_GCM/1000., thick=0.1,color=8,linestyle=3
oplot, temp_gcm_convadj(*,lt_plotindex_gcm_convadj3), altitudes_GCM/1000., thick=0.1,color=8,linestyle=3
oplot, temp_gcm_convadj(*,lt_plotindex_gcm_convadj4), altitudes_GCM/1000., thick=0.1,color=8,linestyle=3
endif


PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; *** Pressions ***

print, '........ PRESSURES'

;what_I_plot = make_array(nZmx)
;labels=['TH P 1d, lt='+string(lt_plot)]
;title_user = TestCase+SubCase+' Pressure Comparisons'
;filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Pressure.ps'
;PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
;CS, SCALE=28
;GSET, XMIN=600, XMAX=900, YMIN=0, YMAX=0.5, TITLE=title_user
;cols=INDGEN(1)+2
;GPLOT, X=pplay(*,lt_plotindex_gcm), Y=-alog(pplay(*,lt_plotindex_gcm)/pGround), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
;AXES, XSTEP = 25, XTITLE='Log P', YSTEP=0.1, YTITLE='Altitude (km)',NDECS=1
;
;oplot, pt(*,lt_plotindex_les),-alog(pt(*,lt_plotindex_les)/pGround), psym=4
;
;PSCLOSE, /NOVIEW
;
;spawn, 'ps2png '+filename

what_I_plot = make_array(nZmx)
labels=['TH P 1d, lt='+string(lt_plot)]
title_user = TestCase+SubCase+LayerCase+' Pressure Comparisons'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Pressure.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=400, XMAX=900, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=pplay(*,lt_plotindex_gcm), Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 25, XTITLE='Pression (Pa)', YSTEP=1, YTITLE='Altitude (km)',NDECS=1

oplot, pt(*,lt_plotindex_les),altitudes_LES/1000., psym=4

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


; *** Temperatures potentielles ***
if(full eq 'true') then begin

xmin = 185
xmax = 240
if (TestCase eq 'Case_C') then begin
xmin = 200
xmax = 255
endif
if (TestCase eq 'Case_I') then begin
xmin = 195
xmax = 250
endif
if (TestCase eq 'Case_Z') then begin
xmin = 230
xmax = 290
endif
if (TestCase eq 'ExtremeCase') then begin
xmin = 195
xmax = 260
endif

print, '........ POTENTIAL TEMPERATURES'
what_I_plot = [[reform(teta_gcm(*,lt_plotindex_gcm_ini))],[reform(teta_gcm(*,lt_plotindex_gcm0))],[reform(teta_gcm(*,lt_plotindex_gcm))],[reform(teta_gcm(*,lt_plotindex_gcm2))],[reform(teta_gcm(*,lt_plotindex_gcm3))],[reform(teta_gcm(*,lt_plotindex_gcm4))]]
labels=['TH teta 1d, lt='+string(lt_plot_ini),'TH teta 1d, lt='+string(lt_plot0),'TH teta 1d, lt='+string(lt_plot),'TH teta 1d, lt='+string(lt_plot2),'TH teta 1d, lt='+string(lt_plot3),'TH teta 1d, lt='+string(lt_plot4)]
title_user = TestCase+SubCase+LayerCase+' Teta comparisons (recomputed from T and P)'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Teta.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=xmin, XMAX=xmax, YMIN=0, YMAX=1.4, TITLE=title_user
cols=INDGEN(6)+2
GPLOT, X=what_I_plot, Y=-alog(pplay(*,lt_plotindex_gcm)/pGround), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 5, XTITLE='Potential temperature (K)', YSTEP=0.2, YTITLE='-Log(P/P0) ',NDECS=1

oplot, teta_les(*,lt_plotindex_les_ini), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les0), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les2), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les3), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les4), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
if(overplot_convadj EQ 'true') then begin
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj0), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj0)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj2), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj2)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj3), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj3)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj4), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj4)/pGround), thick=0.1,color=8,linestyle=3
endif
oplot, teta_gcm_0(*), -alog(pplay(*,lt_plotindex_gcm)/pGround)

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

xmin = 230
xmax = 245
what_I_plot = [[reform(teta_gcm(*,lt_plotindex_gcm)),tsurf_gcm(lt_plotindex_gcm)],[reform(teta_gcm(*,lt_plotindex_gcm2)),tsurf_gcm(lt_plotindex_gcm2)]]
labels=['TH teta 1d, lt='+string(lt_plot),'TH teta 1d, lt='+string(lt_plot2)]
title_user = TestCase+SubCase+LayerCase+' Teta comparisons (recomputed from T and P)'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Teta_zoom.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=xmin, XMAX=xmax, YMIN=0, YMAX=0.05, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=what_I_plot, Y=[0.,-alog(pplay(*,lt_plotindex_gcm)/pGround)], /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 5, XTITLE='Potential temperature (K)', YSTEP=0.005, YTITLE='-Log(P/P0) ',NDECS=1

;oplot, teta_les(*,lt_plotindex_les_ini), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
;oplot, teta_les(*,lt_plotindex_les0), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les2), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
;oplot, teta_les(*,lt_plotindex_les3), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
;oplot, teta_les(*,lt_plotindex_les4), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
if(overplot_convadj EQ 'true') then begin
;oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj0), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj0)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj2), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj2)/pGround), thick=0.1,color=8,linestyle=3
;oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj3), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj3)/pGround), thick=0.1,color=8,linestyle=3
;oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj4), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj4)/pGround), thick=0.1,color=8,linestyle=3
endif
;oplot, teta_gcm_0(*), -alog(pplay(*,lt_plotindex_gcm)/pGround)

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


if (plot_3d eq 'true') then begin
what_I_plot = [[reform(teta_gcm(*,lt_plotindex_gcm_ini))],[reform(teta_gcm(*,lt_plotindex_gcm0))],[reform(teta_gcm(*,lt_plotindex_gcm))],[reform(teta_gcm(*,lt_plotindex_gcm2))],[reform(teta_gcm(*,lt_plotindex_gcm3))],[reform(teta_gcm(*,lt_plotindex_gcm4))]]
labels=['TH teta 1d, lt='+string(lt_plot_ini),'TH teta 1d, lt='+string(lt_plot0),'TH teta 1d, lt='+string(lt_plot),'TH teta 1d, lt='+string(lt_plot2),'TH teta 1d, lt='+string(lt_plot3),'TH teta 1d, lt='+string(lt_plot4)]
title_user = TestCase+SubCase+LayerCase+' Teta comparisons (recomputed from T and P)'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Teta.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=xmin, XMAX=xmax, YMIN=0, YMAX=1.4, TITLE=title_user
cols=INDGEN(6)+2
GPLOT, X=what_I_plot, Y=-alog(pplay_3d(*,lt_plotindex_gcm)/pGround), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 5, XTITLE='Potential temperature (K)', YSTEP=0.2, YTITLE='-Log(P/P0) ',NDECS=1

oplot, teta_les(*,lt_plotindex_les_ini), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les0), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les2), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les3), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
oplot, teta_les(*,lt_plotindex_les4), -alog(pt(*,lt_plotindex_les)/pGround), psym=4

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename
endif

endif else begin

print, '........ POTENTIAL TEMPERATURES'

what_I_plot = reform(teta_gcm(*,lt_plotindex_gcm))
labels=['TH teta 1d, lt='+string(lt_plot)]
title_user = TestCase+SubCase+LayerCase+' Teta comparisons (recomputed from T and P)'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Teta.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=210, XMAX=240, YMIN=0, YMAX=2, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=-alog(pplay(*,lt_plotindex_gcm)/pGround), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 5, XTITLE='Potential temperature (K)', YSTEP=0.2, YTITLE='-Log(P/P0) ',NDECS=1

oplot, teta_les(*,lt_plotindex_les),-alog(pt(*,lt_plotindex_les)/pGround), psym=4
PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename
endelse

print, '........ SURFACE TEMPERATURES'
what_I_plot = tsurf_gcm
labels=['TH 1d tsurf']
title_user = TestCase+SubCase+LayerCase+' Surface temperatures (recomputed from T and P)'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_tsurf.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=200, YMAX=300, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime_gcm, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=5., YTITLE='Surface temperature',NDECS=1

oplot, localtime, tsurf, psym=4
PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


; ---------------------- *** Vitesses verticales *** -------------------------------
; ------------ Verification de l'approx terrestre wmax = vmoy dans la couche instable

print, '........ CHECKING wmax = vmoy in unstable layer'

what_I_plot = uv_moy(*,lt_plotindex_les)
labels=['LES uv_moy']
title_user = TestCase+SubCase+LayerCase+' LES mean UV comp to max W in plume trac1'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_UV.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0, XMAX=10, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, YSTEP = 1, YTITLE='Altitude (km)', XSTEP=1, XTITLE='Mean horizontal velocity inside domain (m/s)',NDECS=1

oplot, make_array(nz,value=wmax(lt_plotindex_les)), altitudes_LES/1000., psym=4

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; ------------ Profil de vitesse

print, '........ VERTICAL VELOCITY'

what_I_plot = make_array(nZmx)
FOR k=0, nZmx-2 DO BEGIN
        what_I_plot(k) = 0.5*(zw2_lev(k,lt_plotindex_gcm) + zw2_lev(k+1,lt_plotindex_gcm))
ENDFOR
what_I_plot(nZmx-1) = 0.

smoothed_w_mean1_les = make_array(nz)
smoothed_w_mean2_les = make_array(nz)
smoothed_w_mean1_down_les = make_array(nz)
FOR t=-ns,ns DO BEGIN
        smoothed_w_mean1_les = smoothed_w_mean1_les + REFORM(w_mean1(*,lt_plotindex_les+t))
        smoothed_w_mean2_les = smoothed_w_mean2_les + REFORM(w_mean2(*,lt_plotindex_les+t))
	smoothed_w_mean1_down_les = smoothed_w_mean1_down_les + REFORM(w_mean1_down(*,lt_plotindex_les+t))
ENDFOR

smoothed_w_mean1_les = smoothed_w_mean1_les/nstot
smoothed_w_mean2_les = smoothed_w_mean2_les/nstot
smoothed_w_mean1_down_les = smoothed_w_mean1_down_les/nstot

ratio = make_array(nz)
FOR k=0, nz-1 DO BEGIN
	IF(smoothed_w_mean1_les(k) ne 0.) then ratio(k) = smoothed_w_mean1_down_les(k)/smoothed_w_mean1_les(k) else ratio(k)=0.
ENDFOR

labels=['TH 1d w, lt='+string(lt_plot)]
title_user = TestCase+SubCase+LayerCase+' Vertical velocity comparisons (inside thermal)'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Wprofile.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-6, XMAX=8, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, YSTEP = 2, YTITLE='Altitude (km)', XSTEP=1, XTITLE='Vertical velocity inside thermal (m/s)',NDECS=1

oplot, smoothed_w_mean1_les, altitudes_LES/1000., psym=4
oplot, smoothed_w_mean2_les, altitudes_LES/1000., psym=5
oplot, smoothed_w_mean1_down_les, altitudes_LES/1000., psym=6

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


; *** Static stability ***

print, '........ STATIC STABILITY'

dteta_dz_gcm = deriv(altitudes_GCM,reform(teta_gcm(*,lt_plotindex_gcm)))
dteta_dz_les = deriv(altitudes_LES,reform(teta_les(*,lt_plotindex_les)))

what_I_plot = dteta_dz_gcm
labels=['TH static stability 1d']
title_user = TestCase+SubCase+LayerCase+' Static stability comparison'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_dTetadz.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.002, XMAX=0.006, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.0005 , XTITLE='Static stability (K.m-1)', YSTEP=1, YTITLE='Altitude (km)',NDECS=4

oplot, dteta_dz_les, altitudes_LES/1000., psym=4

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

print,' -----------------------------------------------------------------------------------------------------------------------'
print,' ***  LES diagnostics of the PLUME *** MUAHAHAHAHAHA'
print,' V2 with E and D computed for a UDE plume'
print,' -----------------------------------------------------------------------------------------------------------------------'

print, '........ EXTRACTING DATA'

; --- Reinterpolation of F

fm_therm_gcm_interlay = make_array(nZmx,nTmx)

FOR k=0, nZmx-2 DO BEGIN
        fm_therm_gcm_interlay(k,*) = (fm_therm_gcm_lev(k,*) + fm_therm_gcm_lev(k+1,*))/2.
ENDFOR
fm_therm_gcm_interlay(nZmx-1,*)=0.

; --- Calculation of gcm df/dz using entrainment and detrainments and NOT F

df_dz_gcm = deriv(altitudes_GCM,reform(fm_therm_gcm_interlay(*,lt_plotindex_gcm)))
 
; --- Smoothing of the mass flux on a user-defined window

smoothed_fm_trac1_les = make_array(nz)
smoothed_fm_trac2_les = make_array(nz)
smoothed_downward_fm_trac1_les = make_array(nz)
FOR t=-ns,ns DO BEGIN
	smoothed_fm_trac1_les = smoothed_fm_trac1_les + REFORM(fm_trac1_les(*,lt_plotindex_les+t))
	smoothed_fm_trac2_les = smoothed_fm_trac2_les + REFORM(fm_trac2_les(*,lt_plotindex_les+t))
	smoothed_downward_fm_trac1_les = smoothed_downward_fm_trac1_les + REFORM(downward_flux1(*,lt_plotindex_les+t))
ENDFOR

smoothed_fm_trac1_les = smoothed_fm_trac1_les/nstot
smoothed_fm_trac2_les = smoothed_fm_trac2_les/nstot
smoothed_downward_fm_trac1_les = smoothed_downward_fm_trac1_les/nstot

; --- Calculation of the entrainement rate according to Rio(2010)
; done in the heavy part at the begining (reeaaaally heavy)

; --- Smoothing of the entrainment on a ~30min window

; term 1


smoothed_e_term1_trac1_les = make_array(nz)
smoothed_e_term1_trac2_les = make_array(nz)
FOR t=-ns,ns DO BEGIN
        smoothed_e_term1_trac1_les = smoothed_e_term1_trac1_les + REFORM(e_trac1_les(*,lt_plotindex_les+t))
        smoothed_e_term1_trac2_les = smoothed_e_term1_trac2_les + REFORM(e_trac2_les(*,lt_plotindex_les+t))
ENDFOR

smoothed_e_term1_trac1_les = smoothed_e_term1_trac1_les/nstot
smoothed_e_term1_trac2_les = smoothed_e_term1_trac2_les/nstot

smoothed_e_rate_term1_trac1_les = make_array(nz)
smoothed_e_rate_trac2_les = smoothed_e_term1_trac2_les

; it already is an entrainment rate ! KIND OF : it is E/Mc, and Mc is not F !! NOW it is Mc/deltaTeta * dchi/dz
FOR k=0, nz-1 DO BEGIN
	IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_e_rate_term1_trac1_les(k) = smoothed_e_term1_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_e_rate_term1_trac1_les(k) =0.
;        IF (smoothed_fm_trac2_les(k) ne 0.) THEN smoothed_e_rate_trac2_les(k) = smoothed_e_trac2_les(k)/smoothed_fm_trac2_les(k) ELSE smoothed_e_rate_trac2_les(k) =0.
ENDFOR

; term 2  & 3


smoothed_e_term2_trac1_les = make_array(nz)
smoothed_e_term3_trac1_les = make_array(nz)
smoothed_e_term1_ude_trac1_les = make_array(nz)
smoothed_e_term2_ude_trac1_les = make_array(nz)
smoothed_e_term3_ude_trac1_les = make_array(nz)

FOR t=-ns,ns DO BEGIN
        smoothed_e_term2_trac1_les = smoothed_e_term2_trac1_les + REFORM(e1_term2(*,lt_plotindex_les+t))
        smoothed_e_term2_trac1_les = smoothed_e_term2_trac1_les + REFORM(e1_term2(*,lt_plotindex_les+t))
	smoothed_e_term1_ude_trac1_les = smoothed_e_term1_ude_trac1_les + REFORM(e1_term1_ude(*,lt_plotindex_les+t))
        smoothed_e_term2_ude_trac1_les = smoothed_e_term2_ude_trac1_les + REFORM(e1_term2_ude(*,lt_plotindex_les+t))
        smoothed_e_term3_ude_trac1_les = smoothed_e_term3_ude_trac1_les + REFORM(e1_term3_ude(*,lt_plotindex_les+t))
ENDFOR

smoothed_e_term2_trac1_les = smoothed_e_term2_trac1_les/nstot
smoothed_e_term3_trac1_les = smoothed_e_term3_trac1_les/nstot
smoothed_e_term1_ude_trac1_les = smoothed_e_term1_ude_trac1_les/nstot
smoothed_e_term2_ude_trac1_les = smoothed_e_term2_ude_trac1_les/nstot
smoothed_e_term3_ude_trac1_les = smoothed_e_term3_ude_trac1_les/nstot

smoothed_e_rate_term2_trac1_les = make_array(nz)
smoothed_e_rate_term3_trac1_les = make_array(nz)
smoothed_e_rate_term1_ude_trac1_les = make_array(nz)
smoothed_e_rate_term2_ude_trac1_les = make_array(nz)
smoothed_e_rate_term3_ude_trac1_les = make_array(nz)

FOR k=0, nz-1 DO BEGIN
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_e_rate_term2_trac1_les(k) = smoothed_e_term2_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_e_rate_term2_trac1_les(k) =0.
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_e_rate_term3_trac1_les(k) = smoothed_e_term3_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_e_rate_term3_trac1_les(k) =0.
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_e_rate_term1_ude_trac1_les(k) = smoothed_e_term1_ude_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_e_rate_term1_ude_trac1_les(k) =0.
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_e_rate_term2_ude_trac1_les(k) = smoothed_e_term2_ude_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_e_rate_term2_ude_trac1_les(k) =0.
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_e_rate_term3_ude_trac1_les(k) = smoothed_e_term3_ude_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_e_rate_term3_ude_trac1_les(k) =0.

ENDFOR

; --- Summing...

smoothed_e_rate_trac1_les = smoothed_e_rate_term1_trac1_les + smoothed_e_rate_term2_trac1_les + smoothed_e_rate_term3_trac1_les
smoothed_e_rate_ude_trac1_les = smoothed_e_rate_term1_ude_trac1_les + smoothed_e_rate_term2_ude_trac1_les + smoothed_e_rate_term3_ude_trac1_les

;print, 'ommiting term3'
;smoothed_e_rate_trac1_les = smoothed_e_rate_term1_trac1_les + smoothed_e_rate_term2_trac1_les

; --- Smoothing of the detrainment rate

smoothed_d_term1_trac1_les = make_array(nz)
smoothed_d_term2_trac1_les = make_array(nz)
smoothed_d_term3_trac1_les = make_array(nz)
smoothed_d_term1_ude_trac1_les = make_array(nz)
smoothed_d_term2_ude_trac1_les = make_array(nz)
smoothed_d_term3_ude_trac1_les = make_array(nz)

FOR t=-ns,ns DO BEGIN
        smoothed_d_term1_trac1_les = smoothed_d_term1_trac1_les + REFORM(d1_term1(*,lt_plotindex_les+t))
        smoothed_d_term2_trac1_les = smoothed_d_term2_trac1_les + REFORM(d1_term2(*,lt_plotindex_les+t))
        smoothed_d_term3_trac1_les = smoothed_d_term3_trac1_les + REFORM(d1_term3(*,lt_plotindex_les+t))
        smoothed_d_term1_ude_trac1_les = smoothed_d_term1_ude_trac1_les + REFORM(d1_term1_ude(*,lt_plotindex_les+t))
        smoothed_d_term2_ude_trac1_les = smoothed_d_term2_ude_trac1_les + REFORM(d1_term2_ude(*,lt_plotindex_les+t))
        smoothed_d_term3_ude_trac1_les = smoothed_d_term3_ude_trac1_les + REFORM(d1_term3_ude(*,lt_plotindex_les+t))
ENDFOR

smoothed_d_term1_trac1_les = smoothed_d_term1_trac1_les/nstot
smoothed_d_term2_trac1_les = smoothed_d_term2_trac1_les/nstot
smoothed_d_term3_trac1_les = smoothed_d_term3_trac1_les/nstot
smoothed_d_term1_ude_trac1_les = smoothed_d_term1_ude_trac1_les/nstot
smoothed_d_term2_ude_trac1_les = smoothed_d_term2_ude_trac1_les/nstot
smoothed_d_term3_ude_trac1_les = smoothed_d_term3_ude_trac1_les/nstot

smoothed_d_rate_term1_trac1_les = make_array(nz)
smoothed_d_rate_term2_trac1_les = make_array(nz)
smoothed_d_rate_term3_trac1_les = make_array(nz)
smoothed_d_rate_term1_ude_trac1_les = make_array(nz)
smoothed_d_rate_term2_ude_trac1_les = make_array(nz)
smoothed_d_rate_term3_ude_trac1_les = make_array(nz)

FOR k=0, nz-1 DO BEGIN
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_d_rate_term1_trac1_les(k) = smoothed_d_term1_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_d_rate_term1_trac1_les(k) =0.
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_d_rate_term2_trac1_les(k) = smoothed_d_term2_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_d_rate_term2_trac1_les(k) =0.
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_d_rate_term3_trac1_les(k) = smoothed_d_term3_trac1_les(k)/smoothed_fm_trac1_les(k) ELSE smoothed_d_rate_term3_trac1_les(k) =0.
        IF (smoothed_fm_trac1_les(k) ne 0.) THEN BEGIN
	smoothed_d_rate_term1_ude_trac1_les(k) = smoothed_d_term1_ude_trac1_les(k)/smoothed_fm_trac1_les(k)
	smoothed_d_rate_term2_ude_trac1_les(k) = smoothed_d_term2_ude_trac1_les(k)/smoothed_fm_trac1_les(k)
	smoothed_d_rate_term3_ude_trac1_les(k) = smoothed_d_term3_ude_trac1_les(k)/smoothed_fm_trac1_les(k)
	ENDIF ELSE BEGIN
	smoothed_d_rate_term1_ude_trac1_les(k)=0.
	smoothed_d_rate_term2_ude_trac1_les(k)=0.
	smoothed_d_rate_term3_ude_trac1_les(k)=0.
	ENDELSE
ENDFOR

; --- Summing...

smoothed_d_rate_trac1_les = smoothed_d_rate_term1_trac1_les+smoothed_d_rate_term2_trac1_les+smoothed_d_rate_term3_trac1_les
smoothed_d_rate_ude_trac1_les = smoothed_d_rate_term1_ude_trac1_les+smoothed_d_rate_term2_ude_trac1_les+smoothed_d_rate_term3_ude_trac1_les
;print, 'ommiting term3'
;smoothed_d_rate_trac1_les = smoothed_d_rate_term1_trac1_les+smoothed_d_rate_term2_trac1_les

; --- PLOTTING : BUOYANCY TERM

smoothed_buoyancy_trac1_les = make_array(nz)
smoothed_buoyancy_ude_trac1_les = make_array(nz)
smoothed_buoyancy_trac2_les = make_array(nz)
smoothed_buoyancy_downdraft1_les_ude = make_array(nz)

FOR t=-ns,ns DO BEGIN
        smoothed_buoyancy_trac1_les = smoothed_buoyancy_trac1_les + REFORM(buoyancy1_les(*,lt_plotindex_les+t))
	smoothed_buoyancy_ude_trac1_les = smoothed_buoyancy_ude_trac1_les + REFORM(grav*(tplume1moy(*,lt_plotindex_les+t)/tenv1moy_ude(*,lt_plotindex_les+t)-1.))
        smoothed_buoyancy_trac2_les = smoothed_buoyancy_trac2_les + REFORM(buoyancy2_les(*,lt_plotindex_les+t))
	smoothed_buoyancy_downdraft1_les_ude = smoothed_buoyancy_downdraft1_les_ude + REFORM(grav*(tdown1moy(*,lt_plotindex_les+t)/tenv1moy_ude(*,lt_plotindex_les+t)-1.))
ENDFOR

smoothed_buoyancy_trac1_les = smoothed_buoyancy_trac1_les/nstot
smoothed_buoyancy_ude_trac1_les = smoothed_buoyancy_ude_trac1_les/nstot
smoothed_buoyancy_trac2_les = smoothed_buoyancy_trac2_les/nstot
smoothed_buoyancy_downdraft1_les_ude = smoothed_buoyancy_downdraft1_les_ude/nstot

print, '........ BUOYANCY'

what_I_plot = [[buoyancy_gcm(*,lt_plotindex_gcm)],[buoyancy_est_gcm(*,lt_plotindex_gcm)]]
labels=['TH buoyancy term','TH estimated buoyancy in plume']
title_user = TestCase+SubCase+LayerCase+' UDE plume buoyancy'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_B.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.06, XMAX=0.06, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.01 , XTITLE='N.m-1', YSTEP=1, YTITLE=' Altitude (km)',NDECS=3

;oplot, smoothed_buoyancy_trac1_les, altitudes_LES/1000., psym=4
;oplot, smoothed_buoyancy_trac2_les, altitudes_LES/1000., psym=5
print, smoothed_buoyancy_ude_trac1_les
oplot, smoothed_buoyancy_ude_trac1_les, altitudes_LES/1000., psym=4
oplot, smoothed_buoyancy_downdraft1_les_ude, altitudes_LES/1000., psym=7

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


; --- PLOTTING : MASS FLUX 

print, '........ MASS FLUX'

f_gcm = fm_therm_gcm_interlay(*,lt_plotindex_gcm)
what_I_plot = f_gcm
labels=['TH mass flux']
title_user = TestCase+SubCase+LayerCase+' mass flux comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_f.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.008, XMAX=0.008, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.002 , XTITLE='Kg.m-2.s-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=4

oplot, smoothed_fm_trac1_les, altitudes_LES/1000., psym=4
oplot, smoothed_fm_trac2_les, altitudes_LES/1000., psym=5
oplot, smoothed_downward_fm_trac1_les, altitudes_LES/1000., psym=6

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


; --- PLOTTING : MASS FLUX DERIVATIVE

print, '........ MASS FLUX DERIVATIVE'

df_dz_les1 = deriv(altitudes_LES,reform(smoothed_fm_trac1_les))
df_dz_les2 = deriv(altitudes_LES,reform(smoothed_fm_trac2_les))

what_I_plot = df_dz_gcm
labels=['TH mass flux vertical derivative']
title_user = TestCase+SubCase+LayerCase+' mass flux vertical derivative comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_dfdz.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.00002, XMAX=0.00002, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.000005 , XTITLE='Kg.m-3.s-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=6
oplot, df_dz_les1, altitudes_LES/1000., psym=4
oplot, df_dz_les2, altitudes_LES/1000., psym=5

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; --- PLOTTING : ENTRAINMENT RATE e = E/f 

print, '........ ENTRAINMENT RATE'

e_gcm = make_array(nZmx)

FOR k=0, nZmx-1 DO BEGIN
	IF (fm_therm_gcm_interlay(k,lt_plotindex_gcm) ne 0.) THEN BEGIN
		e_gcm(k) = zdz_entr_therm_gcm(k,lt_plotindex_gcm)/(approx_zdz_gcm(k)*fm_therm_gcm_interlay(k,lt_plotindex_gcm))
	ENDIF ELSE BEGIN
		e_gcm(k) = 0.
	ENDELSE
ENDFOR


what_I_plot = e_gcm
labels=['TH entrainment rate']
title_user = TestCase+SubCase+LayerCase+' UDE entrainment rate comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_e.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.003, XMAX=0.003, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.0006 , XTITLE='entrainment rate m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=4

oplot, smoothed_e_rate_ude_trac1_les, altitudes_LES/1000., psym=4
;oplot, smoothed_e_rate_trac1_les, altitudes_LES/1000., psym=4
;oplot, smoothed_e_rate_trac2_les, altitudes_LES/1000., psym=5

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; --- PLOTTING : EXTENDED ENTRAINMENT RATE 

print, '........ EXTENDED ENTRAINMENT RATE'

;what_I_plot = [[smoothed_e_rate_term1_trac1_les],[smoothed_e_rate_term2_trac1_les],[smoothed_e_rate_term3_trac1_les],[smoothed_e_rate_trac1_les]]
what_I_plot = [[smoothed_e_rate_term1_ude_trac1_les],[smoothed_e_rate_term2_ude_trac1_les],[smoothed_e_rate_term3_ude_trac1_les],[smoothed_e_rate_ude_trac1_les]]
labels=['LES base entrainment rate','LES term2 e rate','LES term3 e rate','LES total e rate']
title_user = TestCase+SubCase+LayerCase+' UDE entrainment rate comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_e_terms.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.01, XMAX=0.01, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(4)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.005 , XTITLE='entrainment rate m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=3

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


; --- PLOTTING : EXTENDED DETRAINMENT RATE 

print, '........ EXTENDED DETRAINMENT RATE'

what_I_plot = [[smoothed_d_rate_term1_ude_trac1_les],[smoothed_d_rate_term2_ude_trac1_les],[smoothed_d_rate_term3_ude_trac1_les],[smoothed_d_rate_ude_trac1_les]]
labels=['LES term 1 detrainment rate','LES term2 d rate','LES term3 d rate','LES Total d rate']
title_user = TestCase+SubCase+LayerCase+' UDE detrainment rate comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_d_terms.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.01, XMAX=0.01, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(4)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.002 , XTITLE='detrainment rate m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=4

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; --- PLOTTING : DETRAINMENT RATE d = D/f

print, '........ DETRAINMENT RATE'

;smoothed_d_rate_trac2_les = make_array(nz)
;
;FOR k=0, nz-1 DO BEGIN
;        IF (smoothed_fm_trac1_les(k) ne 0.) THEN smoothed_d_rate_trac1_les(k) = smoothed_e_rate_trac1_les(k) - df_dz_les1(k)/smoothed_fm_trac1_les(k) ELSE smoothed_d_rate_trac1_les(k) =0.
;        IF (smoothed_fm_trac2_les(k) ne 0.) THEN smoothed_d_rate_trac2_les(k) = smoothed_e_rate_trac2_les(k) - df_dz_les2(k)/smoothed_fm_trac2_les(k) ELSE smoothed_d_rate_trac2_les(k) =0.
;ENDFOR
;
d_gcm = make_array(nZmx)
FOR k=0, nZmx-1 DO BEGIN
        IF (fm_therm_gcm_interlay(k,lt_plotindex_gcm) ne 0.) THEN BEGIN
                d_gcm(k) = zdz_detr_therm_gcm(k,lt_plotindex_gcm)/(approx_zdz_gcm(k)*fm_therm_gcm_interlay(k,lt_plotindex_gcm))
        ENDIF ELSE BEGIN
                d_gcm(k) = 0.
        ENDELSE
ENDFOR

what_I_plot = d_gcm
labels=['TH detrainment rate']
title_user = TestCase+SubCase+LayerCase+' UDE detrainment rate comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_d.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0., XMAX=0.03, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.005 , XTITLE='m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=2

oplot, smoothed_d_rate_ude_trac1_les, altitudes_LES/1000., psym=4
;oplot, smoothed_d_rate_trac1_les, altitudes_LES/1000., psym=4
;oplot, smoothed_d_rate_trac2_les, altitudes_LES/1000., psym=5

PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

; --- PLOTTING : FRACTION COVERAGE

print, '........ EXTENDED ALPHA'

what_I_plot = alpha_interlay_gcm
labels=['TH alpha']
title_user = TestCase+SubCase+LayerCase+' fraction coverage comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_alpha.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0., XMAX=1., YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.1 , XTITLE='m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=1

oplot, smoothed_alpha1_les, altitudes_LES/1000., psym=4
oplot, smoothed_alpha2_les, altitudes_LES/1000., psym=5
oplot, smoothed_beta1_les, altitudes_LES/1000., psym=6

PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

; --- PLOTTING : THEORETICAL ENTRAINMENT RATE FROM LES DATA 

print, '........ PARAMETRIZED RATES'

approx_zdz_les = make_array(nz)

approx_zdz_les(0)=altitudes_LES(1)
FOR k=1, nz-2 DO BEGIN
	approx_zdz_les(k) = altitudes_LES(k+1) - altitudes_LES(k)
ENDFOR
approx_zdz_les(nz-1)=approx_zdz_les(nz-2)


theoretical_e_trac1_les = make_array(nz)
theoretical_e_trac2_les = make_array(nz)


FOR k=0, nz-1 DO BEGIN
        theoretical_e_trac1_les(k) = MAX([0.,(betalpha/(1.+betalpha))*((afact*smoothed_buoyancy_trac1_les(k)/((smoothed_w_mean1_les(k))^2.)) - fact_epsilon)])
        theoretical_e_trac2_les(k) = MAX([0.,(betalpha/(1.+betalpha))*((afact*smoothed_buoyancy_trac2_les(k)/((smoothed_w_mean2_les(k))^2.)) - fact_epsilon)])
ENDFOR


what_I_plot = [[theoretical_e_trac1_les],[theoretical_e_trac2_les]]
labels=['LES TH theo e rate trac1','LES TH theo e rate trac2']
title_user = TestCase+SubCase+LayerCase+' comp. theor. entr. rate comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_e_theoretical.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.015, XMAX=0.03, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.003 , XTITLE='m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=3

oplot, smoothed_e_rate_trac1_les, altitudes_LES/1000., psym=4
oplot, smoothed_e_rate_trac2_les, altitudes_LES/1000., psym=5

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; --- PLOTTING : THEORETICAL DETRAINMENT RATE FROM LES DATA
; ZDZ STUFF REMOVED

print, '........ PARAMETRIZED DETRAINMENT'

theoretical_d_trac1_les = make_array(nz)
theoretical_d_trac2_les = make_array(nz)

FOR k=0, nz-1 DO BEGIN
	theoretical_d_trac1_les(k) = MAX([detr_min,-afact*(betalpha/(1.+betalpha))*(smoothed_buoyancy_trac1_les(k)/((smoothed_w_mean1_les(k))^2.))])
	theoretical_d_trac2_les(k) = MAX([detr_min,-afact*(betalpha/(1.+betalpha))*(smoothed_buoyancy_trac2_les(k)/((smoothed_w_mean2_les(k))^2.))])
ENDFOR

what_I_plot = [[theoretical_d_trac1_les],[theoretical_d_trac2_les]]
labels=['LES TH theo d rate trac1','LES TH theo d rate trac2']
title_user = TestCase+SubCase+LayerCase+' comp. theor. detr. rate comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_d_theoretical.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.1, XMAX=0.1, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.01 , XTITLE='m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=2

oplot, smoothed_d_rate_trac1_les, altitudes_LES/1000., psym=4
;oplot, smoothed_d_rate_trac2_les, altitudes_LES/1000., psym=5

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; --- PLOTTING : THEORETICAL E-D  RATE FROM LES DATA

print, '........ PARAMETRIZED MASS FLUX DERIVATIVE'

theoretical_dfdz_f_trac1_les = make_array(nz)
theoretical_dfdz_f_trac2_les = make_array(nz)

theoretical_dfdz_f_trac1_les = theoretical_e_trac1_les - theoretical_d_trac1_les
theoretical_dfdz_f_trac2_les = theoretical_e_trac2_les - theoretical_d_trac2_les

df_dz_f_les1 = make_array(nz)
df_dz_f_les2 = make_array(nz)

FOR k=0, nz-1 DO BEGIN
	IF (smoothed_fm_trac1_les(k) ne 0.) THEN df_dz_f_les1(k) = df_dz_les1(k)/smoothed_fm_trac1_les(k) ELSE df_dz_f_les1(k)=0.
        IF (smoothed_fm_trac2_les(k) ne 0.) THEN df_dz_f_les2(k) = df_dz_les2(k)/smoothed_fm_trac2_les(k) ELSE df_dz_f_les2(k)=0.
ENDFOR

what_I_plot = [[theoretical_dfdz_f_trac1_les],[theoretical_dfdz_f_trac2_les]]
labels=['LES TH theo 1/f df/dz trac1','LES TH theo 1/f df/dz trac2']
title_user = TestCase+SubCase+LayerCase+' comp. theor. entr. - detr. rate comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_dfdzf_theoretical.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.02, XMAX=0.02, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.01 , XTITLE='entr - detr (rates) m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=4

oplot, df_dz_f_les1, altitudes_LES/1000., psym=4
oplot, df_dz_f_les2, altitudes_LES/1000., psym=5

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; --- PLOTTING : e versus B/w2

print, '........ BUOYANCY AND VERTICAL VELOCITY ENTRAINMENT RATE DEPENDENCY'

B_w2_trac1 = make_array(nz)
B_w2_trac2 = make_array(nz)
;dwdz_trac1 = deriv(altitudes_LES,smoothed_w_mean1_les)
;dwdz_trac2 = deriv(altitudes_LES,smoothed_w_mean2_les)
;full_dwdz_trac1 = make_array(nz,nttot)
;full_dadz_trac1 = make_array(nz,nttot)
;FOR l=0, nttot -1 DO BEGIN
;	full_dwdz_trac1(*,l) = deriv(altitudes_LES,w_mean1(*,l))
;	full_dadz_trac1(*,l) = deriv(altitudes_LES,alpha1out(*,l))
;ENDFOR
;alpha = 0.

FOR k=0, nz-1 DO BEGIN
	IF (smoothed_e_rate_trac1_les(k) ne 0.) THEN B_w2_trac1(k) = smoothed_buoyancy_ude_trac1_les(k)/(smoothed_w_mean1_les(k))^2 ELSE B_w2_trac1(k)=-0.
;	IF (smoothed_e_rate_trac2_les(k) ne 0.) THEN B_w2_trac2(k) = smoothed_buoyancy_trac2_les(k)/(smoothed_w_mean2_les(k))^2 ELSE B_w2_trac2(k)=-0.
ENDFOR

;FOR zzz=0.,30 DO BEGIN

;alpha = zzz/10.
;
;FOR k=0, nz-1 DO BEGIN
;        IF (smoothed_e_rate_trac1_les(k) ne 0. and smoothed_w_mean1_les(k) ne 0.) THEN B_w2_trac1(k) = 0.5*(smoothed_buoyancy_trac1_les(k)/(smoothed_w_mean1_les(k))^2 - alpha*(1./smoothed_w_mean1_les(k))*dwdz_trac1(k)) ELSE B_w2_trac1(k)=0.
;        IF (smoothed_e_rate_trac2_les(k) ne 0. and smoothed_w_mean2_les(k) ne 0.) THEN B_w2_trac2(k) = 0.5*(smoothed_buoyancy_trac2_les(k)/(smoothed_w_mean2_les(k))^2 - alpha*(1./smoothed_w_mean2_les(k))*dwdz_trac2(k)) ELSE B_w2_trac2(k)=0.
;ENDFOR

;print, smoothed_buoyancy_trac1_les(*)/(smoothed_w_mean1_les(*))^2
;print, (1./smoothed_w_mean1_les(*))*dwdz_trac1(*)

full_e1 = make_array(nz,nttot)
full_bw2 = make_array(nz,nttot)
FOR k=0, nz-1 DO BEGIN
FOR l=0, nttot-1 DO BEGIN
	if(fm_trac1_les(k,l) ne 0.) then full_e1(k,l)=(e1_term1_ude(k,l)+e1_term2_ude(k,l)+e1_term3_ude(k,l))/fm_trac1_les(k,l) else full_e1(k,l)=0.
        if(w_mean1(k,l) ne 0.) then full_bw2(k,l)=grav*(tplume1moy(k,l)/tenv1moy_ude(k,l) -1.)/(w_mean1(k,l)^2) else full_bw2(k,l)=-0.
;        if(w_mean1(k,l) ne 0.) then full_bw2(k,l)=0.5*(alpha*buoyancy1_les(k,l)/(w_mean1(k,l)^2) - full_dwdz_trac1(k,l)/w_mean1(k,l)) else full_bw2(k,l)=0.
ENDFOR
ENDFOR

what_I_plot = smoothed_e_rate_ude_trac1_les
labels=['e_rate trac1']
title_user = TestCase+SubCase+LayerCase+' LES UDE entrainment rate dep with B/w2, average over '+taverage+' mn,' 
;filename = TestCase+SubCase+LayerCase+string(alpha,format='(F3.1)')+'Gcm_Les_Comp_e_Bw2.ps'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_e_Bw2.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0., XMAX=0.4, YMIN=0., YMAX=0.4, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=B_w2_trac1, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30, SYM=5, /NOLINES
AXES, XSTEP = 0.05 , XTITLE='LES Entrainment rate m-1', YSTEP=0.05, YTITLE='Parametrized entrainement rate m-1',NDECS=4

;oplot, smoothed_e_rate_trac2_les, B_w2_trac2, psym=5
FOR l=0, nttot-1 DO BEGIN
	oplot, full_e1(*,l),full_bw2(*,l),thick=0.05,psym=1
ENDFOR
;mean_full_e1 = make_array(nz) & mean_full_bw2 = make_array(nz)
;FOR k=0, nz-1 DO BEGIN
;	mean_full_e1(k) = MEAN(reform(full_e1(k,*)))
;	mean_full_bw2(k) = MEAN(reform(full_bw2(k,*)))
;ENDFOR
;oplot, mean_full_e1, mean_full_bw2, thick=0.3, psym = 2,color=5
;oplot, theoretical_e_trac1_les, B_w2_trac1,psym=2,thick=0.8,color=7
oplot,(B_w2_trac1)/2.2222 + 0.0005,B_w2_trac1,thick=0.3,color=7
;oplot, 0.0118*(B_w2_trac1/0.043)^(1./1.65),B_w2_trac1,thick=0.3,color=7
oplot, 0.0113*(B_w2_trac1/0.043)^(1./1.65),B_w2_trac1,thick=0.3,color=7
oplot, 0.0105*(B_w2_trac1/0.048)^(1./1.7),B_w2_trac1,thick=0.3,color=6
;oplot, alog((B_w2_trac1 - 0.000942361)/0.0444855) - 3.85453, B_w2_trac1, thick=0.3,color=7

;print, alog((B_w2_trac1)/0.0444855) - 3.85453

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

;ENDFOR

what_I_plot = full_bw2(*,lt_plotindex_les)
labels=['B/w2']
title_user = TestCase+SubCase+LayerCase+' LES UDE B/w2, average over '+taverage+' mn,'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Bw2.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.01, XMAX=0.01, YMIN=0., YMAX=6., TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.002 , XTITLE='B/w2 term in LES UDE', YSTEP=0.5, YTITLE='Altitude (km)',NDECS=4

;FOR l=0, nttot-1 DO BEGIN
;        oplot, full_bw2(*,l),altitudes_LES/1000.,thick=0.05,psym=1
;ENDFOR

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


print, '........ BUOYANCY AND VERTICAL VELOCITY DETRAINMENT RATE DEPENDENCY'

full_d1 = make_array(nz,nttot)
;full_dSiebesma = make_array(nz,nttot)


FOR k=0, nz-1 DO BEGIN
FOR l=0, nttot-1 DO BEGIN
        if(fm_trac1_les(k,l) ne 0.) then full_d1(k,l)=(d1_term1_ude(k,l)+d1_term2_ude(k,l)+d1_term3_ude(k,l))/fm_trac1_les(k,l) else full_d1(k,l)=-0.
;        if(w_mean1(k,l) ne 0.) then full_dSiebesma(k,l)=0.5*buoyancy1_les(k,l)/(w_mean1(k,l)^2) -1.5*full_dwdz_trac1(k,l)/w_mean1(k,l) - full_dadz_trac1(k,l)/alpha1out(k,l) else full_dSiebesma(k,l)=-0.
ENDFOR
ENDFOR

what_I_plot = smoothed_d_rate_ude_trac1_les
labels=['d_rate trac1']
title_user = TestCase+SubCase+LayerCase+' LES UDE detrainment rate dep with B/w2, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_d_Bw2.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0, XMAX=0.4, YMIN=0., YMAX=0.4, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=full_bw2(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30, SYM=5, /NOLINES
AXES, XSTEP = 0.05 , XTITLE='Detrainment rate m-1', YSTEP=0.05, YTITLE='B/w²',NDECS=4

FOR l=0, nttot-1 DO BEGIN
        oplot, full_d1(*,l),full_bw2(*,l),thick=0.05,psym=1
ENDFOR
;oplot, theoretical_d_trac1_les, full_bw2(*,lt_plotindex_les),psym=2,thick=0.8,color=7
;oplot,B_w2_trac1/2.7 + 0.0002,B_w2_trac1,thick=0.3,color=7
oplot,B_w2_trac1/2.222 + 0.0002,B_w2_trac1,thick=0.3,color=7
oplot, 0.0105*(B_w2_trac1/0.048)^(1./1.7),B_w2_trac1,thick=0.3,color=7
oplot, 0.0118*(B_w2_trac1/0.043)^(1./1.65),B_w2_trac1,thick=0.3,color=6

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


; --- PLOTTING : 0.5*dwu2/dz
print, '........ ///////////// starting local thermal model ///////////'

print, ' -> alimentation'
a_star = make_array(nz, value=0.)
f_star = make_array(nz, value=0.)
f_star(0) = 1.
teta_est = make_array(nz, value=0.)
teta_p = make_array(nz, value=0.)
zw2 = make_array(nz,value=0.)
w_est = make_array(nz,value=0.)
entr_star = make_array(nz,value=0.)
detr_star = make_array(nz,value=0.)
zbuoy_est = make_array(nz,value=0.)
zbuoy = make_array(nz,value=0.)
a_star_tot = 0.

zw2(1)= 0.4*0.3811552*2*grav*(teta_les(0,lt_plotindex_les)/teta_les(1,lt_plotindex_les) -1.)*approx_zdz_les(0)
w_est(1) = zw2(1)
FOR k=0, nz-2 DO BEGIN
        if ((teta_les(k,lt_plotindex_les) GT (teta_les(k+1,lt_plotindex_les) +0.1)) AND (teta_les(0,lt_plotindex_les) GE teta_les(k,lt_plotindex_les))) then begin
	a_star(k) = MAX([(teta_les(k,lt_plotindex_les)-teta_les(k+1,lt_plotindex_les)),0.])*sqrt(altitudes_LES(k))
	a_star_tot = a_star_tot + a_star(k)
        lalim = k+1.
        endif
ENDFOR
FOR k=0, nz-1 DO BEGIN
        a_star(k) = a_star(k)/a_star_tot
ENDFOR
print, 'alimentation :'
;print, a_star
f_star(0)=0.
f_star(1) = a_star(0)
teta_p(0) = teta_les(0,lt_plotindex_les)
teta_est(0) = teta_les(0,lt_plotindex_les)
print, ' -> plume'
FOR k=1, nz-2 DO BEGIN
	if (k LT lalim) then begin
		teta_est(k) = (f_star(k)*teta_p(k-1)+a_star(k)*0.5*(teta_les(k,lt_plotindex_les) + teta_p(k-1)))/(f_star(k) + a_star(k))
	endif else begin
		teta_est(k) = teta_p(k-1)
	endelse
	zbuoy_est(k) = grav*(teta_est(k)/teta_les(k,lt_plotindex_les) -1.)
	zw2fact=fact_epsilon*2.*approx_zdz_les(k)/(1.+betalpha)
	zdw2=afact*zbuoy_est(k)/fact_epsilon
        w_est(k+1) = MAX([0.0001,exp(-zw2fact)*(w_est(k)-zdw2)+zdw2])
	if (w_est(k+1) lt 0.) then begin
                w_est(k+1)=zw2(k)
        endif
	if (w_est(k+1) gt 0.001) then begin
		entr_star(k)=f_star(k)*approx_zdz_les(k)*(betalpha/(1.+betalpha))*MAX([0.,afact*zbuoy_est(k)/w_est(k+1)])
		detr_star(k)=f_star(k)*approx_zdz_les(k)*MAX([detr_min,-afact*(betalpha/(1.+betalpha))*zbuoy_est(k)/w_est(k+1)])
	endif
        if (k lt lalim) then begin
          a_star(k)=max([a_star(k),entr_star(k)])
          entr_star(k)=0.
        endif
	if (w_est(k+1) gt 0.001) then begin
	f_star(k+1)=f_star(k)+a_star(k)+entr_star(k)-detr_star(k)
	if (k lt lalim) then begin
	teta_p(k)=(f_star(k)*teta_p(k-1)+(a_star(k)+entr_star(k))*0.5*(teta_p(k-1) + teta_les(k,lt_plotindex_les)))/(f_star(k+1)+detr_star(k))
	endif else begin
	teta_p(k)=(f_star(k)*teta_p(k-1)+(a_star(k)+entr_star(k))*teta_les(k,lt_plotindex_les))/(f_star(k+1)+detr_star(k))
	endelse
	zbuoy(k) = grav*(teta_p(k)/teta_les(k,lt_plotindex_les) -1.)
        zw2fact=fact_epsilon*2.*approx_zdz_les(k)/(1.+betalpha)
        zdw2=afact*zbuoy(k)/fact_epsilon
        zw2(k+1) = MAX([0.0001,exp(-zw2fact)*(zw2(k)-zdw2)+zdw2])
	endif
ENDFOR
print, ' -> done'

print, '........ CHECKING VERTICAL VELOCITY FORMULATION'
what_I_plot = make_array(nz,value=0.)
what_I_overplot = make_array(nz,value=0.)
FOR k=0, nz-2 DO BEGIN
        what_I_plot(k) = 0.5*(sqrt(zw2(k)) + sqrt(zw2(k+1)))
ENDFOR
FOR k=0, nZmx-2 DO BEGIN
        what_I_overplot(k) = 0.5*(zw2_lev(k,lt_plotindex_gcm) + zw2_lev(k+1,lt_plotindex_gcm))
ENDFOR
labels=['zw2 in les calc as in TH']
title_user = TestCase+SubCase+LayerCase+' LES vertical velocity formulation check'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_w2_check.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0., XMAX=8., YMIN=0., YMAX=7, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1 , XTITLE='w2 in LES from TH calc m/s', YSTEP=1., YTITLE='Altitude (km)',NDECS=3

oplot, w_mean1(*,lt_plotindex_les), altitudes_LES/1000.
oplot, what_I_overplot, altitudes_GCM/1000.,psym =4

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

print, '........ CHECKING TETA ESTIMATIONS FORMULATION'

what_I_plot = [[teta_est],[teta_p],[tplume1moy(*,lt_plotindex_les)]]
labels=['LES estimated teta','LES teta plume calc as in TH','LES teta plume']
title_user = TestCase+SubCase+LayerCase+' LES estimated teta formulation check'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_teta_check.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=214., XMAX=220., YMIN=0., YMAX=7, TITLE=title_user
cols=INDGEN(3)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1 , XTITLE='Teta plume and Est in LES from TH calc K', YSTEP=1., YTITLE='Altitude (km)',NDECS=3

oplot, teta_gcm(*,lt_plotindex_gcm)*(buoyancy_gcm(*,lt_plotindex_gcm)/grav +1.), altitudes_GCM/1000.
oplot, teta_gcm(*,lt_plotindex_gcm)*(buoyancy_est_gcm(*,lt_plotindex_gcm)/grav +1.), altitudes_GCM/1000.

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

print, '........ CHECKING MASS FLUX FORMULATION'

what_I_plot = [[f_star/MAX(f_star)],[smoothed_fm_trac1_les(*,lt_plotindex_les)/MAX(smoothed_fm_trac1_les(*,lt_plotindex_les))]]
labels=['LES normalized f_star ','LES normalized updraft mass flux']
title_user = TestCase+SubCase+LayerCase+' LES normalized f_star formulation check'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_fm_check.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0., XMAX=1., YMIN=0., YMAX=7, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.1 , XTITLE='f*/max(f*) in LES from TH calc', YSTEP=1., YTITLE='Altitude (km)',NDECS=3

oplot,fm_therm_gcm_interlay(*,lt_plotindex_gcm)/MAX(fm_therm_gcm_interlay(*,lt_plotindex_gcm)), altitudes_GCM/1000.,psym=4

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; COMPUTING THE CONTINUITY EQUATION IN THE QUASI-BOUSSINESQ APPROX

da_dt = make_array(nz,n_elements(localtime))
smoothed_da_dt = make_array(nz)
FOR k=0, nz-1 DO BEGIN
	da_dt(k,*) = deriv(localtime,reform(alpha1out(k,*)))/3700.
ENDFOR
FOR t=-ns,ns DO BEGIN
        smoothed_da_dt = smoothed_da_dt + REFORM(da_dt(*,lt_plotindex_les+t))
ENDFOR
smoothed_da_dt = smoothed_da_dt/nstot

rho = pt/(R*temp_les)
smoothed_rho = make_array(nz)
FOR t=-ns,ns DO BEGIN
        smoothed_rho = smoothed_rho + REFORM(rho(*,lt_plotindex_les+t))
ENDFOR
smoothed_rho = smoothed_rho/nstot

continuity1 = smoothed_rho*smoothed_da_dt + df_dz_les1 - smoothed_e_rate_ude_trac1_les*smoothed_fm_trac1_les + smoothed_d_rate_ude_trac1_les*smoothed_fm_trac1_les
print, '........ CONTINUITY CHECK'

what_I_plot = [[continuity1],[smoothed_rho*smoothed_da_dt],[df_dz_les1],[-smoothed_e_rate_ude_trac1_les*smoothed_fm_trac1_les],[smoothed_d_rate_ude_trac1_les*smoothed_fm_trac1_les]]
labels=['total continuity','rho*da/dt','df/dz','-E','D']
title_user = TestCase+SubCase+LayerCase+' LES UDE continuity check, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_continuity.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.0005, XMAX=0.0005, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(5)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30 
AXES, XSTEP = 0.0001 , XTITLE='kg.m-2.s-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=4
PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; COMPUTING THE E-D TERM FROM THE CONTINUITY EQUATION

eminusd1=make_array(nz)
FOR k=0, nz-1 DO BEGIN
	IF(smoothed_fm_trac1_les(k) ne 0.) THEN eminusd1(k) = (smoothed_rho(k)*smoothed_da_dt(k) - df_dz_les1(k))/smoothed_fm_trac1_les(k) ELSE eminusd1(k)=0.
ENDFOR
what_I_plot = eminusd1
labels=['e-d']
title_user = TestCase+SubCase+LayerCase+' LES e-d, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_EminusD.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.002, XMAX=0.002, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.0005 , XTITLE='kg.m-2.s-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=4
PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; COMPUTING THE TURBULENT FLUX DECOMPOSITION IN PASSIVE ENV AND ACTIVE PLUME
; TO CHECK CONSISTENCY

print, '........ CHECKING CONSISTENCY OF UPDRAFT/ENV DECOMPOSITION'

smoothed_hf1_term1 = make_array(nz)
smoothed_hf1_term2 = make_array(nz)
smoothed_hf1_term3 = make_array(nz)
smoothed_wt = make_array(nz)
FOR t=-ns,ns DO BEGIN
	smoothed_hf1_term1 = smoothed_hf1_term1 + REFORM(hf1_term1(*,lt_plotindex_les+t))
	smoothed_hf1_term2 = smoothed_hf1_term2 + REFORM(hf1_term2(*,lt_plotindex_les+t))
	smoothed_hf1_term3 = smoothed_hf1_term3 + REFORM(hf1_term3(*,lt_plotindex_les+t))
	smoothed_wt = smoothed_wt + REFORM(wt(*,lt_plotindex_les+t))
ENDFOR
smoothed_hf1_term1 = smoothed_hf1_term1/nstot
smoothed_hf1_term2 = smoothed_hf1_term2/nstot
smoothed_hf1_term3 = smoothed_hf1_term3/nstot
smoothed_wt = smoothed_wt/nstot

what_I_plot = [[smoothed_hf1_term1],[smoothed_hf1_term2],[smoothed_hf1_term3],[smoothed_hf1_term1+smoothed_hf1_term2+smoothed_hf1_term3]]
labels=['within plume turbulence','within env. turbulence','organized turbulence','TOTAL']
title_user = TestCase+SubCase+LayerCase+' LES turbulence decomposition check, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_turbu.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-1, XMAX=1.5, YMIN=0, YMAX=6, TITLE=title_user
cols=INDGEN(4)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=9, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.5 , XTITLE='m.K/s', YSTEP=1, YTITLE='Altitude (km)',NDECS=4
oplot, smoothed_wt, altitudes_LES/1000.,psym=3
PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename


; COMPUTING THE TURBULENT FLUX DECOMPOSITION IN PASSIVE ENV ,ACTIVE PLUME and ACTIVE DOWNDRAFT
; TO CHECK CONSISTENCY

print, '........ CHECKING CONSISTENCY OF UPDRAFT/DOWNDRAFT/ENV DECOMPOSITION'

smoothed_hf1_ude_term1 = make_array(nz)
smoothed_hf1_ude_term2 = make_array(nz)
smoothed_hf1_ude_term3 = make_array(nz)
smoothed_hf1_ude_term4 = make_array(nz)
FOR t=-ns,ns DO BEGIN
        smoothed_hf1_ude_term1 = smoothed_hf1_ude_term1 + REFORM(hf1_ude_term1(*,lt_plotindex_les+t))
        smoothed_hf1_ude_term2 = smoothed_hf1_ude_term2 + REFORM(hf1_ude_term2(*,lt_plotindex_les+t))
        smoothed_hf1_ude_term3 = smoothed_hf1_ude_term3 + REFORM(hf1_ude_term3(*,lt_plotindex_les+t))
        smoothed_hf1_ude_term4 = smoothed_hf1_ude_term4 + REFORM(hf1_ude_term4(*,lt_plotindex_les+t))
ENDFOR
smoothed_hf1_ude_term1 = smoothed_hf1_ude_term1/nstot
smoothed_hf1_ude_term2 = smoothed_hf1_ude_term2/nstot
smoothed_hf1_ude_term3 = smoothed_hf1_ude_term3/nstot
smoothed_hf1_ude_term4 = smoothed_hf1_ude_term4/nstot

what_I_plot = [[smoothed_hf1_ude_term1],[smoothed_hf1_ude_term2],[smoothed_hf1_ude_term3],[smoothed_hf1_ude_term4],[smoothed_hf1_ude_term1+smoothed_hf1_ude_term2+smoothed_hf1_ude_term3+smoothed_hf1_ude_term4]]
labels=['within plume turbulence','within downdraft turbulence','within env. turbulence','organized turbulence','TOTAL']
title_user = TestCase+SubCase+LayerCase+' LES UDE turbulence decomposition check, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_turbu_ude.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-1, XMAX=1.5, YMIN=0, YMAX=6, TITLE=title_user
cols=INDGEN(5)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=9, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.25 , XTITLE='m.K/s', YSTEP=0.5, YTITLE='Altitude (km)',NDECS=3
oplot, smoothed_wt, altitudes_LES/1000.,psym=3
PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; CHECK CONSISTENCY OF We = W and THETA e = THETA approximation

print, '........ CHECKING CONSISTENCY OF env variable (w_e,theta_e) = mean variable (w_overbar,theta_overbar)'
print, 'as well as mean(w) = 0, in the UDE decomposition'

smoothed_delta_theta_ude = make_array(nz)
smoothed_delta_w_ude = make_array(nz)
smoothed_w_mean1_full = make_array(nz)

FOR t=-ns,ns DO BEGIN
        smoothed_delta_theta_ude = smoothed_delta_theta_ude + REFORM(w_mean1_env_ude(*,lt_plotindex_les+t)-w_mean1_full(*,lt_plotindex_les+t))
        smoothed_delta_w_ude = smoothed_delta_w_ude + REFORM(tenv1moy_ude(*,lt_plotindex_les+t)-tmoy_full(*,lt_plotindex_les+t))
        smoothed_w_mean1_full = smoothed_w_mean1_full + REFORM(w_mean1_full(*,lt_plotindex_les+t))
ENDFOR

smoothed_delta_theta_ude = smoothed_delta_theta_ude/nstot
smoothed_delta_w_ude = smoothed_delta_w_ude/nstot
smoothed_w_mean1_full = smoothed_w_mean1_full/nstot

what_I_plot = [[smoothed_delta_theta_ude],[smoothed_delta_w_ude],[smoothed_w_mean1_full]]
labels=['theta env_ude - theta moy','w env_ude - w moy','mean w over domain (*,*,k)']
title_user = TestCase+SubCase+LayerCase+' LES UDE env/mean approximation check, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_approx_ude.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-2, XMAX=2, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(3)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.5 , XTITLE='(m/s) and (K)', YSTEP=1, YTITLE='Altitude (km)',NDECS=4
PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; GETTING SOME INSIGHT ON PLUME'S INSIDE TEMPERATURES

print, '........ STRUCTURE POTENTIAL TEMPERATURES'
xmin = 210
xmax = 220
if (TestCase eq 'Case_Z') then begin
xmin = 260
xmax = 270
endif
ztva = teta_gcm*(buoyancy_gcm/grav +1.)
ztva_est = teta_gcm*(buoyancy_est_gcm/grav +1.)
what_I_plot = [[tplume1moy(*,lt_plotindex_les)],[tenv1moy(*,lt_plotindex_les)],[teta_les(*,lt_plotindex_les)]]
labels=['Teta updraft','Teta env ','Teta moy']
title_user = TestCase+SubCase+LayerCase+' LES Teta in the structures, no average'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_fullTeta.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=xmin, XMAX=xmax, YMIN=0, YMAX=6, TITLE=title_user
cols=INDGEN(3)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1 , XTITLE='Potential Temperature (K)', YSTEP=1, YTITLE='Altitude (km)',NDECS=4
oplot, ztva(*,lt_plotindex_gcm), altitudes_GCM/1000., thick=0.3
oplot, ztva_est(*,lt_plotindex_gcm), altitudes_GCM/1000., thick=0.3
oplot, teta_gcm(*,lt_plotindex_gcm), altitudes_GCM/1000., thick=0.3
PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename


print, '........ UDE STRUCTURE POTENTIAL TEMPERATURES'

what_I_plot = [[tplume1moy(*,lt_plotindex_les)],[tdown1moy(*,lt_plotindex_les)],[tenv1moy_ude(*,lt_plotindex_les)],[teta_les(*,lt_plotindex_les)]]
labels=['Teta updraft','Teta downdraft','Teta env (UDE)','Teta moy']
title_user = TestCase+SubCase+LayerCase+' LES UDE Teta in the structures, no average'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_fullTeta_ude.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=210, XMAX=220, YMIN=0, YMAX=6, TITLE=title_user
cols=INDGEN(4)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1 , XTITLE='Potential Temperature (K)', YSTEP=1, YTITLE='Altitude (km)',NDECS=4
PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; -----------------------------------------------------------------------------------------------------------------------
; End of PLUME diagnostics
; -----------------------------------------------------------------------------------------------------------------------

; *** TKE ***

print, '........ TKE'

what_I_plot = reform(tke_gcm(*,lt_plotindex_gcm))
labels=['TH tke 1d']
title_user = TestCase+SubCase+LayerCase+' TKE comparison'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_tke.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-1, XMAX=8, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='Turbulent kinetic energy (kg.m-3)', YSTEP=1, YTITLE='Altitude (km)',NDECS=1

oplot, tke_les(*,lt_plotindex_les), altitudes_LES/1000., psym=4

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; *** HEAT FLUX ***

print, '........ HEAT FLUX'

lay_heatFlux_up = make_array(nZmx,nTmx)
lay_heatFlux_down = make_array(nZmx,nTmx)

FOR k=1, nZmx-1 DO BEGIN
	lay_heatFlux_up(k,*) = 0.5*(heatFlux_up(k,*) + heatFlux_up(k-1,*))
	lay_heatFlux_down(k,*) = 0.5*(heatFlux_down(k,*) + heatFlux_down(k-1,*))
ENDFOR
lay_heatFlux_up(0,*)=0.5*(heatFlux_up(0,*))
lay_heatFlux_down(0,*)=0.5*(heatFlux_down(0,*))
zkh_gcm_int = make_array(nZmx)

FOR k=0, nZmx-2 DO BEGIN
	zkh_gcm_int(k) = 0.5*(zkh(k,lt_plotindex_gcm) + zkh(k+1,lt_plotindex_gcm))
ENDFOR
MY_gcm = -zkh_gcm_int*(deriv(altitudes_GCM,reform(zh(*,lt_plotindex_gcm)))); - 0.025*max(lay_heatFlux_up(*,lt_plotindex_gcm)))


;what_I_plot = [[lay_heatFlux_up(*,lt_plotindex_gcm)],[MY_gcm],[lay_heatFlux_up(*,lt_plotindex_gcm)+MY_gcm]]
what_I_plot = [[lay_heatFlux_up(*,lt_plotindex_gcm)],[MY_gcm],[lay_heatFlux_down(*,lt_plotindex_gcm)],[lay_heatFlux_down(*,lt_plotindex_gcm)+lay_heatFlux_up(*,lt_plotindex_gcm)+MY_gcm]]

;labels=['TH updraft heat flux','Mellor and Yamada gcm heat flux','Total']
labels=['TH updraft heat flux','Mellor and Yamada gcm heat flux','TH downdraft heat flux','Total']

title_user = TestCase+SubCase+LayerCase+' TH vertical heat flux'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_WT.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-2, XMAX=3, YMIN=0, YMAX=10, TITLE=title_user
;cols=INDGEN(3)+2
cols=INDGEN(4)+2
GPLOT, X=what_I_plot, Y=altitudes_GCM/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='vertical turbulent heat flux', YSTEP=1, YTITLE='Altitude (km)',NDECS=1

;hf1_ude_term4 = alpha1out*(w_mean1 - w_mean1_full)*(tplume1moy - tmoy_full) + beta1out*(w_mean1_down - w_mean1_full)*(tdown1moy - tmoy_full) + (1.- (alpha1out+beta1out))*(w_mean1_env_ude - w_mean1_full)*(tenv1moy_ude - tmoy_full)

oplot, alpha1out(*,lt_plotindex_les)*(w_mean1(*,lt_plotindex_les) - w_mean1_full(*,lt_plotindex_les))*(tplume1moy(*,lt_plotindex_les) - tmoy_full(*,lt_plotindex_les)), altitudes_LES/1000., color=2

oplot, beta1out(*,lt_plotindex_les)*(w_mean1_down(*,lt_plotindex_les) - w_mean1_full(*,lt_plotindex_les))*(tdown1moy(*,lt_plotindex_les) - tmoy_full(*,lt_plotindex_les)), altitudes_LES/1000.,color=6

oplot, smoothed_hf1_ude_term4, altitudes_LES/1000., color = 5

;oplot, smoothed_hf1_term1, altitudes_LES/1000.,thick=0.1,LINESTYLE = 5
;oplot, smoothed_hf1_term2, altitudes_LES/1000.,color=2,thick=0.1,LINESTYLE = 5
;oplot, smoothed_hf1_term3, altitudes_LES/1000.,color=8,thick=0.1,LINESTYLE = 5

oplot, wt(*,lt_plotindex_les), altitudes_LES/1000. 

;oplot, smoothed_hf1_term3(*,lt_plotindex_les), altitudes_LES/1000., psym=7

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

; *** TRACERS ***

print, '........ TRACER PLOTS DEACTIVATED'

; trying stuff 

buoyancy_downdraft = grav*(tdown1moy/tenv1moy_ude-1.)
lmix = make_array(nttot,value=-1.)
altitudes_rel_LES = make_array(nz,nttot)
FOR l=0, nttot-1 DO BEGIN
;	kmax = where(w_mean1(*,l) eq max(w_mean1(*,l)))
;	if (kmax(0) ne -1) then lmix(l) = altitudes_LES(kmax(0)) else lmix(l) = -1.
;FOR k=nz-2, 1,-1 DO BEGIN
;	if ((buoyancy_downdraft(k,l) gt 0.) and (buoyancy_downdraft(k-1,l) lt 0.)) then lmix(l) = 0.5*(altitudes_LES(k)+altitudes_LES(k+1))
;ENDFOR
FOR k=nz-2, 1,-1 DO BEGIN
	if (tdown1moy(k,l) eq 0.) then lmix(l) = altitudes_LES(k)
ENDFOR
ENDFOR

FOR l=0, nttot-1 DO BEGIN
FOR k=0, nz-1 DO BEGIN
	altitudes_rel_LES(k,l) = altitudes_LES(k)/lmix(l)
ENDFOR
ENDFOR

print, '........ Teta down / Teta up in UDE'

stuff2=tdown1moy/tenv1moy_ude

what_I_plot = stuff2(*,lt_plotindex_les)
labels=['Teta d/Teta env 12h']
title_user = TestCase+SubCase+LayerCase+' TH trying stuff'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_stuff2.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0.992, XMAX=1.004, YMIN=0, YMAX=1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_rel_LES(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.002, XTITLE='Teta d/ Teta env', YSTEP=0.2, YTITLE='Altitude/zi ',NDECS=3

FOR i=0,nttot-1 DO BEGIN
	if(lmix(i) ne -1) then oplot, stuff2(*,i), altitudes_rel_LES(*,i), thick=0.1
ENDFOR
oplot, (altitudes_rel_LES(*,lt_plotindex_les)-0.075)/187.931 + 0.9977, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
oplot, (altitudes_rel_LES(*,lt_plotindex_les))/19.231 + 0.9938, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
oplot, (altitudes_rel_LES(*,lt_plotindex_les)-0.60)/(-1333) + 1.00025, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

print, '........ Teta down / Teta up in UDE'

stuff2=tdown1moy/tplume1moy

what_I_plot = stuff2(*,lt_plotindex_les)
labels=['Teta d/Teta u 12h']
title_user = TestCase+SubCase+LayerCase+' TH trying stuff'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_stuff2.5.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0.95, XMAX=1.1, YMIN=0, YMAX=1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_rel_LES(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.01, XTITLE='Teta d/ Teta u', YSTEP=0.2, YTITLE='Altitude/zi ',NDECS=3

FOR i=0,nttot-1 DO BEGIN
        if(lmix(i) ne -1) then oplot, stuff2(*,i), altitudes_rel_LES(*,i), thick=0.1
ENDFOR
;oplot, (altitudes_rel_LES(*,lt_plotindex_les)-0.075)/187.931 + 0.9977, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
;oplot, (altitudes_rel_LES(*,lt_plotindex_les))/19.231 + 0.9938, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
;oplot, (altitudes_rel_LES(*,lt_plotindex_les)-0.60)/(-1333) + 1.00025, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

print, '........ B down / B up in UDE'

stuff2 = (tdown1moy/tenv1moy_ude -1.)/(tplume1moy/tenv1moy_ude -1.)
what_I_plot = stuff2(*,lt_plotindex_les)
labels=['B down/B up 12h']
title_user = TestCase+SubCase+LayerCase+' TH trying stuff Bd/Bu'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_stuffBuBd.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-1, XMAX=1., YMIN=0, YMAX=1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_rel_LES(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.1, XTITLE='B down/ B up', YSTEP=0.1, YTITLE='Altitude/zi ',NDECS=1

FOR i=0,nttot-1 DO BEGIN
        if(lmix(i) ne -1) then oplot, stuff2(*,i), altitudes_rel_LES(*,i), thick=0.1
ENDFOR
;oplot, ((altitudes_rel_LES(*,lt_plotindex_les)-0.06)/0.839841)^2 - 0.3, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
oplot, ((altitudes_rel_LES(*,lt_plotindex_les)-0.06)/1.16847)^2 - 0.3, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
oplot, (altitudes_rel_LES(*,lt_plotindex_les)-0.7)/1., altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
;oplot, (altitudes_rel_LES(*,lt_plotindex_les))/0.08333333-1., altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
oplot, sqrt(altitudes_rel_LES(*,lt_plotindex_les)/0.122449)-1., altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

print, '........ F down / F up in UDE'

stuff2 = downward_flux1/fm_trac1_les
what_I_plot = stuff2(*,lt_plotindex_les)
labels=['f down/f up 12h']
title_user = TestCase+SubCase+LayerCase+' TH trying stuff f down/f up'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_stufffufd.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-8, XMAX=0.5, YMIN=0, YMAX=1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_rel_LES(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.5, XTITLE='f down/ f up', YSTEP=0.1, YTITLE='Altitude/zi ',NDECS=1

FOR i=0,nttot-1 DO BEGIN
        if(lmix(i) ne -1) then oplot, stuff2(*,i), altitudes_rel_LES(*,i), thick=0.1
ENDFOR
oplot, -alog(((altitudes_rel_LES(*,lt_plotindex_les)+0.0149259)/0.00333)), altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
;oplot, -alog(((altitudes_rel_LES(*,lt_plotindex_les)+0.02)/0.006)), altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

print, '........ dFdz down / dzFdz up in UDE'
stuff3=make_array(nz,nttot)
FOR l=0, nttot-1 DO BEGIN
	stuff3(*,l) = deriv(altitudes_LES,downward_flux1(*,l))/deriv(altitudes_LES,fm_trac1_les(*,l))
ENDFOR
what_I_plot = stuff3(*,lt_plotindex_les)
labels=['dfdz down/dfdz up 12h']
title_user = TestCase+SubCase+LayerCase+' TH trying stuff f down/f up'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_stuffdfudfd.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-30, XMAX=30, YMIN=0, YMAX=1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_rel_LES(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 5, XTITLE='dfdz down/ dfdz up', YSTEP=0.1, YTITLE='Altitude/zi ',NDECS=1

FOR i=0,nttot-1 DO BEGIN
        if(lmix(i) ne -1) then oplot, stuff3(*,i), altitudes_rel_LES(*,i), thick=0.1
ENDFOR
;oplot, -alog(((altitudes_rel_LES(*,lt_plotindex_les)+0.0149259)/0.00333)), altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
;oplot, -alog(((altitudes_rel_LES(*,lt_plotindex_les)+0.02)/0.006)), altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename
;
;what_I_plot = [[ar_col],[co2_col],[tke_col]] 
;labels=['Ar deviation','Co2 deviation','TKE deviation']
;title_user = TestCase+SubCase+' TH 1d tracer conservation'
;filename = TestCase+SubCase+'Gcm_Les_Comp_tracer.ps'
;PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
;CS, SCALE=28
;GSET, XMIN=localtime_gcm(0), XMAX=localtime_gcm(nTmx-1), YMIN=0, YMAX=2, TITLE=title_user
;cols=INDGEN(3)+2
;GPLOT, X=localtime_gcm, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
;AXES, XSTEP = 2, XTITLE='Local time (h)', YSTEP=0.1, YTITLE='Tracer integrated column mass deviation from origin (%)',NDECS=1
;
;PSCLOSE, /NOVIEW
;
;spawn, 'ps2png '+filename
;
;title_user = TestCase+SubCase+' TH argon propagation from first layer'
;PS_START, file = TestCase+SubCase+'Gcm_Les_Comp_Argon.ps'
;
;what_I_plot = transpose(ar(0:6,*))
;
;maxfield_init = 0.05
;minfield_init = 0
;pal=33
;lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
;lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
;
;section, $
;        what_I_plot, $                          ; 2D field
;        localtime_gcm, $                                ; horizontal coordinate
;        altitudes_gcm(0:6), $                                ; altitude coordinate
;        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
;        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
;;       minspace=minspace, $                    ; minimum value of space window (=0: calculate)
;;       maxspace=maxspace, $                    ; maximum value of space window (=0: calculate)
;;       overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
;;       overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
;;       overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
;;       colors=colors, $                        ; number of colors/levels (32 is default)
;        title_plot=title_user, $                ; title of the plot ('Profile' is default)
;        title_axis=['Martian hour (h)','Height above ground (m)'], $                ; title of the [x,y] axis (['Field','Altitude'] is default)
;        ct=pal, $                               ; color table (33-rainbow is default)
;;       topo=topography, $
;        format=format                           ; format of colorbar annotations ('(F6.2)' is default)
;
;PS_END, /PNG
;
;INTERVAL_VOLUME, supermask1, 0.5, 1.,verts, conn
;conn = TETRA_SURFACE(verts, conn)  
;oRain = OBJ_NEW('IDLgrPolygon', verts, POLYGONS=conn, $  
;   COLOR=[255,255,255], SHADING=1)  
;XOBJVIEW, oRain, BACKGROUND=[150,200,255]

;INTERVAL_VOLUME, supermask2, 0.5, 1.5,verts, conn
;conn = TETRA_SURFACE(verts, conn)
;oRain = OBJ_NEW('IDLgrPolygon', verts, POLYGONS=conn, $
;   COLOR=[255,255,255], SHADING=1)
;XOBJVIEW, oRain, BACKGROUND=[150,200,255]

ENDELSE

print, ''
print, '........ ALL DONE'
print, ''

END
