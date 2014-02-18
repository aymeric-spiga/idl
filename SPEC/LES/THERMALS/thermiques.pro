FUNCTION myfunct, Bw2, A
   G = A[0]*Bw2 + A[1]
   RETURN,[ [G], [Bw2] , [1.0] ]
END

FUNCTION myfunct2, Gammaw2, A
   G = A[0]*(Gammaw2)^A[1]
   RETURN,[ [G], [(Gammaw2)^A[1]] , [G*alog(Gammaw2)] ]
END

FUNCTION myfunct3, Gammaw2, A
   G = A[0]*(Gammaw2)^A[1]+A[2]
   RETURN,[ [G], [(Gammaw2)^A[1]] , [G*alog(Gammaw2)] , [1.0] ]
END

FUNCTION myfunct4, Bw2, A
   G = A*Bw2
   RETURN,[[G], [Bw2]]
END

PRO thermiques

spawn, 'clear'
print, ''
print, '** Thermals Analysis **'
print, ' (usine à gaz) '
print, ''

full='true'
f_offset='false'
overplot_convadj='false'
plot_3d = 'false'
les_special=''
lctu_gcm = 6.
got_swdownz='false'
got_tracer_flux='false'
got_hfx='false'

;got_swdownz='true'
;got_hfx='true'
;got_tracer_flux='true'

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
got_updrafts='false'

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
print, ' 15/ LES HR run (257x257x301)'
print, ' 16/ 1D : 13 layers : gcm to 15km'
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
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 867.5594
goto,label_init
endif
if (newtest eq '4') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter_wind'
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
lctu_gcm = 6.
pGround = 856.997    ;6h
;pGround = 867.5594   ;8h
;got_swdownz='true'
;got_hfx='true'
;got_tracer_flux='true'
goto,label_init
endif
if (newtest eq '5') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter_wind_30'
GcmSubCase = '_wind_30'
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
pGround = 856.997    ;6h
f_offset='false'
goto,label_init
endif
if (newtest eq '6') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter_wind_tau1'
GcmSubCase = '_tau1'
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
pGround = 856.997    ;6h
f_offset='false'
goto,label_init
endif
if (newtest eq '7') then begin
TestCase = 'Case_A'
SubCase = '_4_shorter_wind_tau2'
GcmSubCase = '_tau2'
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 856.997    ;6h
goto,label_init
endif
if (newtest eq '8') then begin
TestCase = 'ExtremeCase'
SubCase = ''
GcmSubCase = ''
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 670.018      ;6h
;pGround = 677.722    ;8h
lctu_gcm = 6.
;got_swdownz='true'
;got_hfx='true'
;got_tracer_flux='true'
goto,label_init
endif
if (newtest eq '9') then begin
TestCase = 'Case_C'
SubCase = '_4_shorter_wind'
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 480.6
goto,label_init
endif
if (newtest eq '10') then begin
TestCase = 'Case_I'
SubCase = '_4_shorter_wind'
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
pGround = 630.
goto,label_init
endif
if (newtest eq '11') then begin
TestCase = 'Case_Z'
SubCase = '_4_shorter_wind'
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
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
if (newtest eq '15') then begin
les_special='_HR'
got_swdownz='true'
got_hfx='true'
got_tracer_flux='true'
goto,label_init
endif
if (newtest eq '16') then begin
LayerCase='_13lev'
got_swdownz='true'
goto,label_init
endif
if (newtest eq '17') then begin
TestCase = 'Exomars'
SubCase = ''
GcmSubCase = ''
s_trac1 = 'qtrac2'
s_trac2 = 'qtrac1'
full = 'true'
f_offset='false'
;pGround = 717.899 ;6h
pGround = 718.
lctu_gcm = 6.
;got_swdownz='true'
;got_hfx='true'
;got_tracer_flux='true'
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

les_path='/san0/acolmd/SIMUS/LES_'+TestCase+SubCase+les_special
;gcm_path='/san0/acolmd/SIMUS/GCM_'+TestCase+LayerCase+'_2'
gcm_path='/san0/acolmd/SIMUS/GCM_'+TestCase+GcmSubCase+LayerCase
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
sigmao= 1.			   ; multiplicative coeff for the computation of Sigma0 in the CS
sigmao_ude = 0.3		   ; number of standard deviation away from mean for the selection of downdraft in UDE

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

id=ncdf_open(filesWRF(0))
NCDF_DIMINQ, id, NCDF_DIMID(id, 'west_east'    ), dummy, nx & NCDF_DIMINQ, id, NCDF_DIMID(id, 'south_north'  ), dummy, ny
NCDF_DIMINQ, id, NCDF_DIMID(id, 'bottom_top'   ), dummy, nz & NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), dummy, nt
NCDF_CLOSE, id
id=ncdf_open(filesWRF(nf-1))  ;; for interrupted runs
NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), dummy, ntlast
NCDF_CLOSE, id
nttot = (nf-1)*nt + ntlast
wt  = fltarr(nz,nttot)
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
Gamma_1 = fltarr(nz,nttot) & Gamma_2 = fltarr(nz,nttot) & Gamma_3 = fltarr(nz,nttot)
Gamma_1_tmp = fltarr(nz,nttot)
ptotprime = fltarr(nx,ny,nz) & anomalptot = fltarr(nx,ny,nz) & dptotprimedztmp  = fltarr(nx,ny,nz)
l=0
tsurf_les = fltarr(nttot)

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
	pt(*,l) = TOTAL(TOTAL(getget(filesWRF(loop), 'PTOT' ,  count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny)
        temp_les(*,l) = t(*,l)*(pt(*,l)/p0)^r_cp
	IF (got_pdt eq 'true') then begin
	exner(*,l) = (pt(*,l)/p0)^r_cp
	dTeta_phys(*,l) = (TOTAL(TOTAL(getget(filesWRF(loop), 'PDT', count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny))/exner(*,l)
	ENDIF
	ph = TEMPORARY(ph) + pht(*,l) / (nttot-1)
	p  = TEMPORARY(p ) + pt(*,l) / nttot
        tsurf_les(l)=TOTAL(TOTAL(getget(filesWRF(loop), 'TSURF', count=[0,0,1], offset=[0,0,loop2]),1),1)/ float(nx) / float(ny)

IF (got_updrafts EQ 'true') THEN BEGIN

        ptotprime(*,*,*) = getget(filesWRF(loop), 'PTOT', count=[0,0,0,1], offset=[0,0,0,loop2])
	FOR k=0, nz-1 DO BEGIN
		rhomoy1(*,l) = TOTAL(TOTAL(reform(((ptotprime(*,*,k)/(R*(t(k,l)+tprime(*,*,k))))*(p0/ptotprime(*,*,k))^r_cp)),1),1)/(float(nx)*float(ny))
		anomalptot(*,*,k) = ptotprime(*,*,k) - pt(k,l)
        ENDFOR
	zqtrac1 = getget(filesWRF(loop), s_trac1, count=[0,0,0,1], offset=[0,0,0,loop2])
	FOR i=0,nx-1 DO BEGIN
		FOR j=0, ny-1 DO BEGIN
			dtempdztmp(i,j,*) = deriv(altitudes_LES, tprime(i,j,*) + t(*,l))
                        dptotprimedztmp(i,j,*) = deriv(altitudes_LES, anomalptot(i,j,*))
		ENDFOR
	ENDFOR
	FOR k=0, nz-1 DO BEGIN
                Gamma_2(k,l) = -(TOTAL(TOTAL(dptotprimedztmp(*,*,k),1),1)/float(nx) /float(ny))/rhomoy1(k,l)
                Gamma_3(k,l) = -grav*(TOTAL(TOTAL(anomalptot(*,*,k),1),1)/float(nx) /float(ny))/pt(k,l)
	ENDFOR
        FOR k=0,nz-1 DO BEGIN
	        anomalqtrac1(*,*,k) = zqtrac1(*,*,k) - TOTAL(TOTAL(REFORM(zqtrac1(*,*,k)),1),1)/ float(nx) / float(ny)
                sigmazqtrac1(k) = STDDEV(REFORM(zqtrac1(*,*,k)))
	        IF (k ne 0) THEN BEGIN
         	        subsampledAltitudes = INTERPOL(altitudes_LES(0:k),findgen(k+1),findgen(decimate*k+1)/decimate)
        	       	sigmazminqtrac1(k) = (sigmao/(altitudes_LES(k)-altitudes_LES(0)))*INT_TABULATED(subsampledAltitudes,INTERPOL(sigmazqtrac1(0:k),altitudes_LES(0:k),subsampledAltitudes),/DOUBLE)
	        ENDIF ELSE BEGIN
        	        sigmazminqtrac1(k) = sigmazqtrac1(k)
	        ENDELSE
;		plumeIndex1 =  WHERE((anomalqtrac1(*,*,k) GT scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)])) AND ((anomalw(k)+wprime(*,*,k)) GT 0.))
;		envIndex1 = WHERE((anomalqtrac1(*,*,k) LE scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)])) OR ((anomalw(k)+wprime(*,*,k)) LE 0.))
		plumeIndex1 =  WHERE(anomalqtrac1(*,*,k) GT scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)]))
		envIndex1 = WHERE(anomalqtrac1(*,*,k) LE scale*MAX([sigmazqtrac1(k),sigmazminqtrac1(k)]))
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
		ENDIF ELSE BEGIN
		FOR n=0,n_elements(plumeIndex1)-1 DO BEGIN
 	        	plumeIndex1out(n,k)=plumeIndex1(n)
		ENDFOR
		FOR n=0,n_elements(envIndex1)-1 DO BEGIN
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
                Gamma_1_tmp (k,l) = alpha1out(k,l)*rhomoy1(k,l)*(wtmp(plumeIndex1)-w_mean1(k,l))^2
		hf1tmp(k,l) = TOTAL((wtmp(plumeIndex1)-w_mean1(k,l))*(ttmp(plumeIndex1)-tplume1moy(k,l)),1) / float(n_elements(plumeIndex1))
		hf1tmpenv(k,l) = TOTAL((wtmp(envIndex1)-w_mean1_env(k,l))*(ttmp(envIndex1)-tenv1moy(k,l)),1) / float(n_elements(envIndex1))
		if (envIndex1_ude(0) ne -1) then hf1tmpenv_ude(k,l) = TOTAL((wtmp(envIndex1_ude)-w_mean1_env_ude(k,l))*(ttmp(envIndex1_ude)-tenv1moy_ude(k,l)),1) / float(n_elements(envIndex1_ude)) else hf1tmpenv_ude(k,l) =0.
		if (downdraft_index1(0) ne -1) then hf1tmp_down(k,l) = TOTAL((wtmp(downdraft_index1)-w_mean1_down(k,l))*(ttmp(downdraft_index1)-tdown1moy(k,l)),1) / float(n_elements(downdraft_index1)) else hf1tmp_down(k,l)=0.
		IF((n_elements(plumeIndex1) + n_elements(envIndex1)) ne float(nx*ny)) then print, 'WARNING : INDEX PROBLEM : plume / env : ', n_elements(plumeIndex1), n_elements(envIndex1)
;		IF((n_elements(plumeIndex1) + n_elements(envIndex1_ude) + n_elements(downdraft_index1)) ne float(nx*ny)) then print, 'WARNING : INDEX PROBLEM : plume / env / downdraft : ', n_elements(plumeIndex1), n_elements(envIndex1_ude), n_elements(downdraft_index1)
		ENDELSE
	ENDFOR
	
        Gamma_1(*,l) = -(1./(alpha1out(*,l)*rhomoy1(*,l)))*deriv(altitudes_LES,Gamma_1_tmp(*,l))
        drhoahfdztmp = deriv(altitudes_LES,rhomoy1(*,l)*alpha1out(*,l)*hf1tmp(*,l))
        drhoahfdztmpDetr = deriv(altitudes_LES,rhomoy1(*,l)*(1.-alpha1out(*,l))*hf1tmpenv(*,l))
        drhoahfdztmpDetr_ude = deriv(altitudes_LES,rhomoy1(*,l)*(1.-alpha1out(*,l)-beta1out(*,l))*hf1tmpenv_ude(*,l))
	
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
        zqtrac2 = getget(filesWRF(loop), s_trac2, count=[0,0,0,1], offset=[0,0,0,loop2])
        FOR k=0,nz-1 DO BEGIN
		anomalqtrac2(*,*,k) = zqtrac2(*,*,k) - TOTAL(TOTAL(REFORM(zqtrac2(*,*,k)),1),1)/ float(nx) / float(ny)
                sigmazqtrac2(k) = STDDEV(zqtrac2(*,*,k))
	        IF (k ne 0) THEN BEGIN
	                subsampledAltitudes = INTERPOL(altitudes_LES(0:k),findgen(k+1),findgen(decimate*k+1)/decimate)
	                sigmazminqtrac2(k) = (sigmao/(altitudes_LES(k)-altitudes_LES(0)))*INT_TABULATED(subsampledAltitudes,INTERPOL(sigmazqtrac2(0:k),altitudes_LES(0:k),subsampledAltitudes),/DOUBLE)
	        ENDIF ELSE BEGIN
        	        sigmazminqtrac2(k) = sigmazqtrac2(k)
	        ENDELSE
                plumeIndex2 = WHERE((anomalqtrac2(*,*,k) GT scale*MAX([sigmazqtrac2(k),sigmazminqtrac2(k)])) AND ((wprime(*,*,k)+anomalw(k)) GT 0.))
                envIndex2 = WHERE((anomalqtrac2(*,*,k) LE scale*MAX([sigmazqtrac2(k),sigmazminqtrac2(k)])) OR ((wprime(*,*,k)+anomalw(k)) LE 0.))
                IF(plumeIndex2(0) EQ -1) THEN BEGIN
		 fm_trac2_les(k,l)=0.
		 e_trac2_les(k,l)=0.
                alpha2out(k,l)=0.
                buoyancy2_les(k,l)=0.
		w_mean2(k,l)=0.
                ENDIF ELSE BEGIN
		alpha2(k) = n_elements(plumeIndex2) / float(nx) / float(ny)
                wprimetmp = reform(reform((anomalw(k)+wprime(*,*,k))),[nx*ny,1])
                w_mean2(k,l) = mean(wprimetmp(plumeIndex2))
		wprimetmp=0.
                fm_trac2_les(k,l) = alpha2(k)*rhomoy1(k,l)*w_mean2(k,l)
                tprimetmp = reform(reform(-tprime(*,*,k)),[nx*ny,1])
                dtempdztmplin = reform(reform(dtempdztmp(*,*,k)),[nx*ny,1])
                e_trac2_les(k,l) = TOTAL((1./(temporary(tprimetmp(plumeIndex2))))*(temporary(dtempdztmplin(plumeIndex2))),1)/n_elements(plumeIndex2)
                alpha2out(k,l)=alpha2(k)
                tfull=reform(tprime(*,*,k)+t(k,l),[nx*ny,1])
                tplume2moy=mean(tfull(plumeIndex2))
                tenv2moy=mean(tfull(envIndex2))
                buoyancy2_les(k,l)=grav*(tplume2moy/tenv2moy-1.)
		ENDELSE
        ENDFOR
        zqtrac2=0.

ENDIF

        wt(*,l) = TOTAL(TOTAL(TEMPORARY(tprime)*TEMPORARY(wprime),1),1) / float(nx) / float(ny)
;	wmax(l) = max(w_mean1(*,l))
        l=l+1
	ENDFOR
   	print, 'file '+string(loop+1,'(I0)'), SYSTIME(1) - timetime, ' s'

ENDFOR

IF (got_updrafts EQ 'true') THEN BEGIN


hf1_term1 = hf1tmp*alpha1out
hf1_term2 = temporary(hf1tmpenv)*(1.-alpha1out)
hf1_term3 = alpha1out*(1.-alpha1out)*(w_mean1 - temporary(w_mean1_env))*(tplume1moy - tenv1moy)

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

ENDIF

ht = TEMPORARY(pht) - hgtu/1000.
save, tsurf_les, w_mean1_env, d1_term1_ude, d1_term2_ude, d1_term3_ude, e1_term1_ude, e1_term2_ude, e1_term3_ude, tplume1moy, tdown1moy, w_mean1_full, tmoy_full, tenv1moy_ude, w_mean1_env_ude, uv_moy, hf1_ude_term1, hf1_ude_term2, hf1_ude_term3, hf1_ude_term4, w_mean1_down, downward_flux1, beta1out, hf1_term1, hf1_term2, hf1_term3, d1_term1, d1_term2, d1_term3, e1_term2, e1_term3, buoyancy1_les, buoyancy2_les, w_mean1, w_mean2, nx, ny, alpha1out, alpha2out, e_trac1_les, e_trac2_les, tke_les, ztke, altitudes_LES, ht, t, p, pt, localtime, xtke, ytke, wt, temp_les, wmax, fm_trac1_les, fm_trac2_les,filename=les_path+'/'+datname

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
;getcdf, file=file1, charvar='zmax', invar=zi_gcm
getcdf, file=file1, charvar='lmax_th', invar=lmax_gcm
getcdf, file=file1, charvar='hfmax_th', invar=hfmax_th1d
getcdf, file=file1, charvar='wmax_th', invar=wmax_th1d

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
lt_plot=11.
endif else begin
offset_localtime = 0.
lt_plot=12.
endelse
localtime=localtime+offset_localtime  
;localtime_gcm=localtime_gcm+history_interval_s_gcm/3700. ; we add the offset from the fact that we output (non sense :) )
localtime_gcm=localtime_gcm
; par contre il faut prendre en compte le fait que la premiere frame du les est decalee de 1 !
localtime=localtime

print, '****************************************************************************************************'
print, 'local time LES'
print, localtime
print, '****************************************************************************************************'
print, 'local time GCM'
print, localtime_gcm
print, '****************************************************************************************************'

time_offset = (ndays-1.)*24.

print, 'TIME STEP LES : ',(localtime(1)-localtime(0))*3700.


lt_plot_ini = 6.
lt_plotindex_les_ini = where(localtime eq lt_plot_ini)
lt_plotindex_gcm_ini = where(localtime_gcm eq (lt_plot_ini+time_offset))
;lt_plotindex_gcm_ini = where(localtime_gcm eq (lt_plot_ini+time_offset+history_interval_s_gcm/3700.))

lt_plot0 = 10.
lt_plotindex_les0 = where(localtime eq lt_plot0)
lt_plotindex_gcm0 = where(localtime_gcm eq (lt_plot0+time_offset))
lt_plotindex_gcm_convadj0 = where(localtime_gcm_convadj eq (lt_plot0+time_offset))

lt_plot0a = 11.
lt_plotindex_les0a = where(localtime eq lt_plot0a)
lt_plotindex_gcm0a = where(localtime_gcm eq (lt_plot0a+time_offset))

lt_plotindex_les = where((localtime lt lt_plot+0.01) and (localtime gt lt_plot-0.01))
lt_plotindex_gcm = where(localtime_gcm eq (lt_plot+time_offset))
lt_plotindex_gcm_convadj = where(localtime_gcm_convadj eq (lt_plot+time_offset))
print, 'lt plotindex les 12h'
print, lt_plotindex_les

lt_plota = 13.
lt_plotindex_lesa = where(localtime eq lt_plota)
lt_plotindex_gcma = where(localtime_gcm eq (lt_plota+time_offset))

lt_plot2 = 14.
lt_plotindex_les2 = where(localtime eq lt_plot2)
lt_plotindex_gcm2 = where(localtime_gcm eq (lt_plot2+time_offset))
lt_plotindex_gcm_convadj2 = where(localtime_gcm_convadj eq (lt_plot2+time_offset))

lt_plot2a = 15.
lt_plotindex_les2a = where(localtime eq lt_plot2a)
lt_plotindex_gcm2a = where(localtime_gcm eq (lt_plot2a+time_offset))

lt_plot3 = 16.
lt_plotindex_les3 = where(localtime eq lt_plot3)
lt_plotindex_gcm3 = where(localtime_gcm eq (lt_plot3+time_offset))
lt_plotindex_gcm_convadj3 = where(localtime_gcm_convadj eq (lt_plot3+time_offset))

lt_plot3a = 17.
lt_plotindex_les3a = where(localtime eq lt_plot3a)
lt_plotindex_gcm3a = where(localtime_gcm eq (lt_plot3a+time_offset))

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
;tke_col = tke_col+1.

; Compute <teta> les

teta_les = temporary(t) 

rho = pt/(R*temp_les)

;print, 'bidouille'
;
;FOR l=0,nTmx -1 DO BEGIn
;print, (1300.*hfmax_th1d(l)/(TOTAL(temp_gcm(0:lmax_gcm(l),l),1)/(lmax_gcm(l)+1.)))/wmax_th1d(l)
;ENDFOR


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
        anomalqtrac2(*,*,k) = zqtrac2(*,*,k) - TOTAL(TOTAL(REFORM(zqtrac2(*,*,k)),1),1)/ float(nx) / float(ny)
        sigmazqtrac2(k) = STDDEV(REFORM(zqtrac2(*,*,k)))
        IF (k ne 0) THEN BEGIN
                subsampledAltitudes = INTERPOL(altitudes_LES(0:k),findgen(k+1),findgen(decimate*k+1)/decimate)
                sigmazminqtrac2(k) = (sigmao/(altitudes_LES(k)-altitudes_LES(0)))*INT_TABULATED(subsampledAltitudes,INTERPOL(sigmazqtrac2(0:k),altitudes_LES(0:k),subsampledAltitudes),/DOUBLE)
	ENDIF ELSE BEGIN
                sigmazminqtrac2(k) = sigmazqtrac2(k)
        ENDELSE
        plumeIndex2 =  WHERE((anomalqtrac2(*,*,k) GT scale*MAX([sigmazqtrac2(k),sigmazminqtrac2(k)])) AND ((anomalw(k)+wprime(*,*,k)) GT 0.))
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
ToBin2 = reform(zqtrac2(*,*,k_out_hist(n)),[nx*ny,1])
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
IF (got_updrafts EQ 'true') THEN BEGIN

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

ENDIF
; =========================================================================

; *** Temperatures ***

if (f_offset eq 'false') then begin
print, '........ TEMPERATURES'

xmin = 190
xmax = 250
if (TestCase eq 'Case_Z') then begin
xmin = 170
xmax = 250
endif
if (TestCase eq 'Case_C') then begin
xmin = 180
xmax = 240
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

endif

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

xmin = 210
xmax = 230
if (TestCase eq 'Case_C') then begin
xmin = 225
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
xmin = 225
xmax = 255
endif
if (TestCase eq 'Exomars') then begin
xmin = 225
xmax = 255
endif

if (f_offset eq 'false') then begin

print, '........ POTENTIAL TEMPERATURES'
;what_I_plot = [[reform(teta_gcm(*,lt_plotindex_gcm_ini))],[reform(teta_gcm(*,lt_plotindex_gcm0))],[reform(teta_gcm(*,lt_plotindex_gcm))],[reform(teta_gcm(*,lt_plotindex_gcm2))],[reform(teta_gcm(*,lt_plotindex_gcm3))],[reform(teta_gcm(*,lt_plotindex_gcm4))]]
what_I_plot = [[reform(teta_gcm(*,lt_plotindex_gcm))],[reform(teta_gcm(*,lt_plotindex_gcm2))],[reform(teta_gcm(*,lt_plotindex_gcm3))]]
labels=['TH teta 1d, lt='+string(lt_plot),'TH teta 1d, lt='+string(lt_plot2),'TH teta 1d, lt='+string(lt_plot3)]
title_user = TestCase+SubCase+LayerCase+' Teta comparisons (recomputed from T and P)'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Teta.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=xmin, XMAX=xmax, YMIN=0, YMAX=0.6, TITLE=title_user
cols=INDGEN(3)+2
GPLOT, X=what_I_plot, Y=-alog(pplay(*,lt_plotindex_gcm)/pGround), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 5, XTITLE='Potential temperature (K)', YSTEP=0.2, YTITLE='-Log(P/P0) ',NDECS=1

;oplot, teta_les(*,lt_plotindex_les_ini), -alog(pt(*,lt_plotindex_les)/pGround), psym=4,SYMSIZE=0.5, thick=0.5
;oplot, teta_les(*,lt_plotindex_les0), -alog(pt(*,lt_plotindex_les)/pGround), psym=4,SYMSIZE=0.5, thick=0.5
oplot, teta_les(*,lt_plotindex_les), -alog(pt(*,lt_plotindex_les)/pGround), psym=4,SYMSIZE=0.5, thick=0.5
oplot, teta_les(*,lt_plotindex_les2), -alog(pt(*,lt_plotindex_les)/pGround), psym=4,SYMSIZE=0.5, thick=0.5
oplot, teta_les(*,lt_plotindex_les3), -alog(pt(*,lt_plotindex_les)/pGround), psym=4,SYMSIZE=0.5, thick=0.5
;oplot, teta_les(*,lt_plotindex_les4), -alog(pt(*,lt_plotindex_les)/pGround), psym=4,SYMSIZE=0.5, thick=0.5
if(overplot_convadj EQ 'true') then begin
;oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj0), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj0)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj2), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj2)/pGround), thick=0.1,color=8,linestyle=3
oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj3), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj3)/pGround), thick=0.1,color=8,linestyle=3
;oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj4), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj4)/pGround), thick=0.1,color=8,linestyle=3
endif
;oplot, teta_gcm_0(*), -alog(pplay(*,lt_plotindex_gcm)/pGround), thick=0.5

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


xmin = 210
xmax = 225

if (TestCase eq 'Case_C') then begin
xmin = 230
xmax = 255
endif

if (TestCase eq 'ExtremeCase') then begin
xmin = 225
xmax = 255
endif

;xmin = 280
;xmax = 300
what_I_plot = [[tsurf_gcm(lt_plotindex_gcm),reform(teta_gcm(*,lt_plotindex_gcm))],[tsurf_gcm(lt_plotindex_gcm2),reform(teta_gcm(*,lt_plotindex_gcm2))]]
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


xmin = 225
xmax = 270
;xmin = 280
;xmax = 300
if (TestCase eq 'Case_C') then begin
xmin = 215
xmax = 240
endif

what_I_plot = [[tsurf_gcm(lt_plotindex_gcm),reform(temp_gcm(*,lt_plotindex_gcm))],[tsurf_gcm(lt_plotindex_gcm2),reform(temp_gcm(*,lt_plotindex_gcm2))]]
ymax=0.05
ystep=0.005
print, 'tsurf gcm :',tsurf_gcm(lt_plotindex_gcm)

if (TestCase eq 'Exomars') then begin
xmin = 235
;xmax = 255
xmax=295
what_I_plot = [[tsurf_gcm(lt_plotindex_gcm),reform(temp_gcm(*,lt_plotindex_gcm))],[10.*tsurf_gcm(lt_plotindex_gcm2),10.*reform(temp_gcm(*,lt_plotindex_gcm2))]]
;what_I_plot = 5.*what_I_plot
;ymax=0.027
ymax=0.003
ystep=0.0005
endif



labels=['TH T 1d, lt='+string(lt_plot),'TH T 1d, lt='+string(lt_plot2)]
title_user = TestCase+SubCase+LayerCase+' T comparisons'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_T_zoom.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=xmin, XMAX=xmax, YMIN=0, YMAX=ymax, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=what_I_plot, Y=[0.,-alog(pplay(*,lt_plotindex_gcm)/pGround)], /LEGEND, LEGPOS=9, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 5, XTITLE='Temperature (K)', YSTEP=ystep, YTITLE='-Log(P/P0) ',NDECS=4

;oplot, temp_les(*,lt_plotindex_les_ini), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
;oplot, temp_les(*,lt_plotindex_les0), -alog(pt(*,lt_plotindex_les)/pGround), psym=4

oplot, [tsurf_gcm(lt_plotindex_gcm),temp_les(*,lt_plotindex_les)], [0.,-alog(pt(*,lt_plotindex_les)/pGround)], psym=4
oplot, temp_les(*,lt_plotindex_les), -alog(pt(*,lt_plotindex_les)/pGround), thick=0.1
;oplot, temp_les(*,lt_plotindex_les), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
;oplot, temp_les(*,lt_plotindex_les2), -alog(pt(*,lt_plotindex_les)/pGround), psym=4

;oplot, teta_les(*,lt_plotindex_les3), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
;oplot, teta_les(*,lt_plotindex_les4), -alog(pt(*,lt_plotindex_les)/pGround), psym=4
if(overplot_convadj EQ 'true') then begin
;oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj0), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj0)/pGround), thick=0.1,color=8,linestyle=3

oplot, temp_gcm_convadj(*,lt_plotindex_gcm_convadj), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj)/pGround), thick=0.1,color=8,psym=2
oplot, temp_gcm_convadj(*,lt_plotindex_gcm_convadj), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj)/pGround), thick=0.1,color=8,linestyle=3
;oplot, temp_gcm_convadj(*,lt_plotindex_gcm_convadj2), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj2)/pGround), thick=0.1,color=8,linestyle=3


;oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj3), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj3)/pGround), thick=0.1,color=8,linestyle=3
;oplot, teta_gcm_convadj(*,lt_plotindex_gcm_convadj4), -alog(pplay_convadj(*,lt_plotindex_gcm_convadj4)/pGround), thick=0.1,color=8,linestyle=3
endif
;oplot, teta_gcm_0(*), -alog(pplay(*,lt_plotindex_gcm)/pGround)

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


endif

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


;getcdf, file=les_path+'/wrfout_d01_9999-01-01_tsurf_gcmsoil.nc', charvar='TSURF', invar=tsurf_les_tmp
;tsurf_les=make_array(nttot)
;FOR l=0,nttot-1 DO BEGIN
;        tsurf_les(l)=TOTAL(TOTAL(tsurf_les_tmp(*,*,l),1),1)/float(n_elements(reform(tsurf_les_tmp(*,0,0)))) /float(n_elements(reform(tsurf_les_tmp(0,*,0))))
;ENDFOR
;tsurf_les_tmp=0.

what_I_plot = tsurf_gcm
labels=['TH 1d tsurf']
title_user = TestCase+SubCase+LayerCase+' Surface temperatures (recomputed from T and P)'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_tsurf.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=150, YMAX=320, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime_gcm, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=5., YTITLE='Surface temperature',NDECS=1

;oplot, localtime, tsurf, psym=4,thick=0.3  ;tsurf les with les soil
;oplot, localtime, tsurf_les, psym=4,thick=0.3  ;tsurf les with gcm soil

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


if (got_swdownz eq 'true') then begin

print, '........ DOWNWARD SOLAR FLUX AT SURFACE'
print, '.. and LW FLUX AT SURFACE'

localtime=localtime-history_interval_s/3700.

;getcdf, file=les_path+'/wrfout_d01_9999-01-01_swdownz.nc', charvar='SWDOWNZ', invar=swdownz_les_tmp
getcdf, file=file1, charvar='fluxsurf_sw', invar=swdownz_gcm
;swdownz_les=make_array(nttot)
;FOR l=0,nttot-1 DO BEGIN
;	swdownz_les(l)=TOTAL(TOTAL(swdownz_les_tmp(*,*,l),1),1)/float(nx) /float(ny)
;ENDFOR
;swdownz_les_tmp=0.

swdownz_les=temporary(swdownz)
swdownz_les_Int=INTERPOL(swdownz_les,localtime,localtime_gcm)
swdownz_gcm_Int=INTERPOL(swdownz_gcm,localtime_gcm,localtime,/QUADRATIC)

what_I_plot = swdownz_gcm-swdownz_les_Int
labels=['TH1D SW flux - LES SW flux']
title_user = TestCase+SubCase+LayerCase+' Difference in Solar Fluxes at surface'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_diffswdownz.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=-20, YMAX=20, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime_gcm, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=2., YTITLE='Solar Flux (W.m-2)',NDECS=1

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

what_I_plot = swdownz_les
labels=['LES SW flux']
title_user = TestCase+SubCase+LayerCase+' Solar Fluxes at surface'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_swdownz.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=0, YMAX=500, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=50., YTITLE='Solar Flux (W.m-2)',NDECS=1

oplot, localtime_gcm, swdownz_gcm, psym=1,thick=0.2

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


;getcdf, file=les_path+'/wrfout_d01_9999-01-01_lwdownz.nc', charvar='LWDOWNZ', invar=lwdownz_les_tmp
getcdf, file=file1, charvar='fluxsurf_lw', invar=lwdownz_gcm
;lwdownz_les=make_array(nttot)
;FOR l=0,nttot-1 DO BEGIN
;        lwdownz_les(l)=TOTAL(TOTAL(lwdownz_les_tmp(*,*,l),1),1)/float(nx) /float(ny)
;ENDFOR
;lwdownz_les_tmp=0.

lwdownz_les=temporary(lwdownz)

what_I_plot = lwdownz_les
labels=['LES LW downward flux']
title_user = TestCase+SubCase+LayerCase+' LW downward fluxes at surface'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_lwdownz.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=-50, YMAX=50, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=10., YTITLE='LW downward surface flux (W.m-2)',NDECS=1

oplot, localtime_gcm, lwdownz_gcm, psym=1,thick=0.2

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename





;getcdf, file=les_path+'/wrfout_d01_9999-01-01_flxgrd.nc', charvar='FLXGRD', invar=fluxgrd_les_tmp
;getcdf, file=les_path+'/wrfout_d01_9999-01-01_flxgrd_gcmsoil.nc', charvar='FLXGRD', invar=fluxgrd_les_tmp
getcdf, file=file1, charvar='fluxgrd', invar=fluxgrd_gcm
;fluxgrd_les=make_array(nttot)
;FOR l=0,nttot-1 DO BEGIN
;        fluxgrd_les(l)=TOTAL(TOTAL(fluxgrd_les_tmp(*,*,l),1),1)/float(n_elements(reform(fluxgrd_les_tmp(*,0,0)))) /float(n_elements(reform(fluxgrd_les_tmp(0,*,0))))
;ENDFOR
;fluxgrd_les_tmp=0.

fluxgrd_les=temporary(flxgrd)

what_I_plot = fluxgrd_les
labels=['LES ground flux']
title_user = TestCase+SubCase+LayerCase+' Ground flux at surface'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_fluxgrd.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=-70, YMAX=40, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=10., YTITLE='Ground flux at surface (W.m-2)',NDECS=1

oplot, localtime_gcm, fluxgrd_gcm, psym=1,thick=0.2

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename





;getcdf, file=les_path+'/wrfout_d01_9999-01-01_flxrad.nc', charvar='FLXRAD', invar=fluxrad_les_tmp
;getcdf, file=les_path+'/wrfout_d01_9999-01-01_flxrad_gcmsoil.nc', charvar='FLXRAD', invar=fluxrad_les_tmp
getcdf, file=file1, charvar='fluxrad', invar=fluxrad_gcm
;fluxrad_les=make_array(nttot)
;FOR l=0,nttot-1 DO BEGIN
;        fluxrad_les(l)=TOTAL(TOTAL(fluxrad_les_tmp(*,*,l),1),1)/float(n_elements(reform(fluxrad_les_tmp(*,0,0)))) /float(n_elements(reform(fluxrad_les_tmp(0,*,0))))
;ENDFOR
;fluxrad_les_tmp=0.

fluxrad_les=temporary(flxrad)

what_I_plot = fluxrad_les
labels=['LES radiation flux at surface']
title_user = TestCase+SubCase+LayerCase+' Radiation flux at surface'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_fluxrad.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=-40, YMAX=80, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=10., YTITLE='Radiation flux at surface (W.m-2)',NDECS=1

oplot, localtime_gcm, fluxrad_gcm, psym=1,thick=0.2

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


localtime=localtime+history_interval_s/3700.

endif

if (got_hfx eq 'true') then begin

print, '........ SENSIBLE HEAT FLUX'

localtime=localtime-history_interval_s/3700.
;getcdf, file=les_path+'/wrfout_d01_9999-01-01_hfx.nc', charvar='HFX', invar=hfx_les_tmp
;getcdf, file=les_path+'/wrfout_d01_9999-01-01_hfx_gcmsoil.nc', charvar='HFX', invar=hfx_les_tmp

hfx_les=temporary(hfx)

getcdf, file=file1, charvar='hfx', invar=hfx_gcm

;hfx_les=make_array(nttot)
;FOR l=0,nttot-1 DO BEGIN
;        hfx_les(l)=TOTAL(TOTAL(hfx_les_tmp(*,*,l),1),1)/float(n_elements(reform(hfx_les_tmp(*,0,0)))) /float(n_elements(reform(hfx_les_tmp(0,*,0))))
;ENDFOR
;hfx_les_tmp=0.
;hfx_les_Int=INTERPOL(hfx_les,localtime,localtime_gcm)
;hfx_gcm_Int=INTERPOL(hfx_gcm,localtime_gcm,localtime,/QUADRATIC)

what_I_plot = hfx_les
labels=['LES sensible heat flux']
title_user = TestCase+SubCase+LayerCase+' Sensible heat Fluxes at surface'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_hfx.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=-50, YMAX=50, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=10., YTITLE='Sensible heat Flux (W.m-2)',NDECS=1

oplot, localtime_gcm, hfx_gcm, psym=1,thick=0.2

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename


if (1 eq 1) then begin

openw, lun, "input_hfx.def", /get_lun
for l=0, nttot-1 do printf, lun, localtime(l),hfx_les(l),fluxrad_les(l)+fluxgrd_les(l)-hfx_les(l), format='((2x,F5.2)(4x,F8.2)(4x,F8.2))'
FREE_LUN, lun
close, lun

endif



what_I_plot = fluxrad_les+fluxgrd_les-hfx_les
labels=['LES flux vdifc']
title_user = TestCase+SubCase+LayerCase+' Flux vdifc at surface'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_vdifc.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=-50, YMAX=50, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=10., YTITLE='Heat Flux (W.m-2)',NDECS=1

oplot, localtime_gcm, fluxrad_gcm+fluxgrd_gcm-hfx_gcm, psym=1,thick=0.2

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename





;what_I_plot = [[fluxrad_les-swdownz_les-lwdownz_les],[fluxgrd_les],[hfx_les],[swdownz_les],[lwdownz_les],[fluxrad_les-fluxgrd_les-hfx_les]]
what_I_plot = [[fluxrad_les],[swdownz_les],[lwdownz_les],[swdownz_les-lwdownz_les]]
;labels=['lw up','grd','hfx','sw down','lw down','total']
labels=['rad','sw down','lw down','sw +lw down']
title_user = TestCase+SubCase+LayerCase+' Fluxes at surface'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_totflux.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=-100, YMAX=100, TITLE=title_user
cols=INDGEN(4)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time (h)', YSTEP=10., YTITLE='Fluxes at surface (W.m-2)',NDECS=1

;oplot, localtime_gcm, fluxrad_gcm, psym=1,thick=0.2

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename










localtime=localtime+history_interval_s/3700.








endif

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

full_d_rate_ude = d1_term1_ude + d1_term2_ude + d1_term3_ude

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


B_w2_trac1 = make_array(nz)

FOR k=0, nz-1 DO BEGIN
        IF (smoothed_e_rate_trac1_les(k) ne 0.) THEN B_w2_trac1(k) = smoothed_buoyancy_ude_trac1_les(k)/(smoothed_w_mean1_les(k))^2 ELSE B_w2_trac1(k)=0.
;       IF (smoothed_e_rate_trac2_les(k) ne 0.) THEN B_w2_trac2(k) = smoothed_buoyancy_trac2_les(k)/(smoothed_w_mean2_les(k))^2 ELSE B_w2_trac2(k)=0.
ENDFOR

df_dz_les1 = deriv(altitudes_LES,reform(smoothed_fm_trac1_les))
df_dz_les2 = deriv(altitudes_LES,reform(smoothed_fm_trac2_les))
df_dz_param = make_array(nz)
df_dz_param2 = make_array(nz)

;dlow=0.0013   ;baseline from continuity equation
dlow=0.0003 ; svn baseline
;dlow=0.0005
dcoeff=-0.3
;dcoeff=-0.4

;aaa1=2.5
;bbb1=0.0015
;Ae1=0.045
;Be1=0.6

aaa1=1.60226
bbb1=0.0006
Ae1=0.0454
Be1=0.57

FOR k=0, nz-1 DO BEGIN
	IF (2.5*B_w2_trac1(k) GE 0.) THEN BEGIN
	        IF ((aaa1*B_w2_trac1(k)-bbb1) GE 0.) THEN BEGIN
;		df_dz_param(k)=smoothed_fm_trac1_les(k)*(Ae1*(aaa1*B_w2_trac1(k)-bbb1)^(Be1) - 0.06*(aaa1*B_w2_trac1(k))^(0.75))
		df_dz_param2(k)=smoothed_fm_trac1_les(k)*(Ae1*(aaa1*B_w2_trac1(k)-bbb1)^(Be1) - MAX([(dcoeff*B_w2_trac1(k) + dlow),0.]))
		ENDIF ELSE BEGIN
;		df_dz_param(k)=smoothed_fm_trac1_les(k)*(-0.06*(aaa1*B_w2_trac1(k))^(0.75))
                df_dz_param2(k)=smoothed_fm_trac1_les(k)*(-MAX([(dcoeff*B_w2_trac1(k) + dlow),0.]))
		ENDELSE
	ENDIF ELSE BEGIN
;        df_dz_param(k)=smoothed_fm_trac1_les(k)*(-0.06*(-aaa1*B_w2_trac1(k))^(0.75))
	df_dz_param2(k)=smoothed_fm_trac1_les(k)*(-MAX([(dcoeff*B_w2_trac1(k) + dlow),0.]))
	ENDELSE
ENDFOR

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
;oplot, df_dz_les2, altitudes_LES/1000., psym=5
oplot, df_dz_param, altitudes_LES/1000., psym=6, color=5
oplot, df_dz_param2, altitudes_LES/1000., psym=6, color=6
;print, "fm*(e-d)"
;print, smoothed_fm_trac1_les*(0.045*(2.5*B_w2_trac1-0.0015)^(0.6) - 0.06*(-2.5*B_w2_trac1)^(0.75))
;print, smoothed_fm_trac1_les
print, B_w2_trac1

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

print, '........ EXTENDED TURBULENCE'

buoyancy1_les_ude = grav*(tplume1moy/tenv1moy_ude -1.)
;Gamma_full = buoyancy1_les_ude + Gamma_1 + Gamma_2 + Gamma_3
Gamma_full = buoyancy1_les + Gamma_1 + Gamma_2 + Gamma_3
;Gamma_full = buoyancy1_les + Gamma_1

time_indices = [[lt_plotindex_les0],[lt_plotindex_les0a],[lt_plotindex_les],[lt_plotindex_lesa],[lt_plotindex_les2],[lt_plotindex_les2a],[lt_plotindex_les3]]
ctime = [[lt_plot0],[lt_plot0a],[lt_plot],[lt_plota],[lt_plot2],[lt_plot2a],[lt_plot3]]

openw, lun, "fit_ab_simple_thermiques_"+TestCase+SubCase+les_special, /get_lun
printf, lun, "    a1     ", "     b1     ","     LT"
FREE_LUN, lun
close, lun

openw, lun, "fit_ab_double_thermiques_"+TestCase+SubCase+les_special, /get_lun
printf, lun, "    a1a     ", "     b1a     ", "      a1b     ", "     b1b     ","     LT"
FREE_LUN, lun
close, lun

FOR ttt=0, n_elements(time_indices)-1 DO BEGIN

smooth_t,input=Gamma_full,nz=nz,ndt=6,t0=time_indices(ttt),output=sm_Gamma_full
smooth_t,input=buoyancy1_les,nz=nz,ndt=6,t0=time_indices(ttt),output=sm_buoyancy1_les
smooth_t,input=Gamma_1,nz=nz,ndt=6,t0=time_indices(ttt),output=sm_Gamma_1
smooth_t,input=Gamma_2,nz=nz,ndt=6,t0=time_indices(ttt),output=sm_Gamma_2
smooth_t,input=Gamma_3,nz=nz,ndt=6,t0=time_indices(ttt),output=sm_Gamma_3
smooth_t,input=w_mean1,nz=nz,ndt=6,t0=time_indices(ttt),output=sm_w_mean1
smooth_t,input=full_d_rate_ude,nz=nz,ndt=6,t0=time_indices(ttt),output=sm_full_d_rate_ude

www = where(w_mean1(*,time_indices(ttt))^2 ne 0.)
nw=n_elements(www)

Y=make_array(nw-4)
X=make_array(nw-4)
FOR zz=0, nw-5 DO BEGIN
;        Y(zz) = Gamma_full(www(zz+2),lt_plotindex_les)/w_mean1(www(zz+2),lt_plotindex_les)^2
;        X(zz) = buoyancy1_les(www(zz+2),lt_plotindex_les)/w_mean1(www(zz+2),lt_plotindex_les)^2

; Approche Rio et al (2010)
        Y(zz) = sm_Gamma_full(www(zz+2))/sm_w_mean1(www(zz+2))^2
        X(zz) = sm_buoyancy1_les(www(zz+2))/sm_w_mean1(www(zz+2))^2

; Approche Gregory et al (2001)
;        Y(zz) = sm_Gamma_full(www(zz+2))/(sm_full_d_rate_ude(www(zz+2))*sm_w_mean1(www(zz+2))^2)
;        X(zz) = sm_buoyancy1_les(www(zz+2))/(sm_full_d_rate_ude(www(zz+2))*sm_w_mean1(www(zz+2))^2)
ENDFOR

A = [2.5,-0.0015]

err_gamma=make_array(nw-4,value=0.001)
FOR zz=floor(nw/4.), nw-5 DO BEGIN
        err_gamma(zz)=0.1
ENDFOR

coefs = lmfit(X,Y,A,/DOUBLE,function_name = 'myfunct', itmax = 500, measure_error = err_gamma)

B = [2.5,-0.0015]

err_gamma=make_array(nw-4,value=0.1)
FOR zz=floor(nw/4.), floor(nw*3./4.-5) DO BEGIN
        err_gamma(zz)=0.001
ENDFOR

coefs = lmfit(X,Y,B,/DOUBLE,function_name = 'myfunct', itmax = 500, measure_error = err_gamma)

C = [1.5,-0.0010]

err_gamma=make_array(nw-4,value=0.001)
FOR zz=floor(nw*3./4.-5), nw-5 DO BEGIN
        err_gamma(zz)=0.1
ENDFOR

coefs = lmfit(X,Y,C,/DOUBLE,function_name = 'myfunct', itmax = 1000, measure_error = err_gamma)

openw, lun, 'fit_ab_simple_thermiques_'+TestCase+SubCase+les_special, /append
printf, lun, C[0],C[1],ctime(ttt), format='((2x,F6.3)(4x,F9.6)(4x,I0))'
close, lun

openw, lun, 'fit_ab_double_thermiques_'+TestCase+SubCase+les_special, /append
printf, lun, A[0],A[1],B[0],B[1],ctime(ttt), format='((2x,F6.3)(4x,F9.6)(4x,F6.3)(4x,F9.6)(4x,I0))'
close, lun


print, '~~~~~> LT: '+string(ctime(ttt))+' <~~~~~~~~~'
print, 'suggested coefs for fit, a,b in alim layer:'
print, A
print, 'suggested coefs for fit, a,b above alim layer:'
print, B
print, 'suggested coefs for uniform fit, a,b:'
print, C

;what_I_plot = [[Gamma_full(*,lt_plotindex_les)],[buoyancy1_les_ude(*,lt_plotindex_les)],[Gamma_1(*,lt_plotindex_les)],[Gamma_2(*,lt_plotindex_les)],[Gamma_3(*,lt_plotindex_les)]]

;what_I_plot = [[Gamma_full(*,lt_plotindex_les)],[buoyancy1_les(*,lt_plotindex_les)],[Gamma_1(*,lt_plotindex_les)],[Gamma_2(*,lt_plotindex_les)],[Gamma_3(*,lt_plotindex_les)]]

what_I_plot = [[sm_Gamma_full(*)],[sm_buoyancy1_les(*)],[sm_Gamma_1(*)],[sm_Gamma_2(*)],[sm_Gamma_3(*)]]
labels=['Tot','B','G1','G2','G3']
title_user = TestCase+SubCase+LayerCase+' UDE turbulence term, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Gamma'+string(ctime(ttt),format='(I0)')+'.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.12, XMAX=0.12, YMIN=0, YMAX=6, TITLE=title_user
cols=INDGEN(5)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.03 , XTITLE='Gamma term m.s-2', YSTEP=1, YTITLE='Altitude (km)',NDECS=4

oplot, -0.0015*w_mean1(*,lt_plotindex_les)^2, altitudes_LES/1000., psym=5, thick=0.05
;oplot, 2.5*buoyancy1_les_ude(*,lt_plotindex_les), altitudes_LES/1000., psym=5, thick=0.05
;oplot, 2.5*buoyancy1_les_ude(*,lt_plotindex_les)-0.0015*w_mean1(*,lt_plotindex_les)^2, altitudes_LES/1000., psym=4, thick=0.3
;oplot, 2.5*buoyancy1_les(*,lt_plotindex_les), altitudes_LES/1000., psym=5, thick=0.05
;oplot, 2.*buoyancy1_les(*,lt_plotindex_les)-0.0012*w_mean1(*,lt_plotindex_les)^2, altitudes_LES/1000., psym=4, thick=0.3
oplot, 2.5*buoyancy1_les(*,lt_plotindex_les)-0.0015*w_mean1(*,lt_plotindex_les)^2, altitudes_LES/1000., psym=4, thick=0.3
;oplot, A[0]*buoyancy1_les(*,lt_plotindex_les)+A[1]*w_mean1(*,lt_plotindex_les)^2, altitudes_LES/1000., thick=0.3
;oplot, B[0]*buoyancy1_les(*,lt_plotindex_les)+B[1]*w_mean1(*,lt_plotindex_les)^2, altitudes_LES/1000., thick=0.3

; Rio et al 2010 :
oplot, A[0]*sm_buoyancy1_les(*)+A[1]*sm_w_mean1(*)^2, altitudes_LES/1000., thick=0.3
oplot, B[0]*sm_buoyancy1_les(*)+B[1]*sm_w_mean1(*)^2, altitudes_LES/1000., thick=0.3
oplot, C[0]*sm_buoyancy1_les(*)+C[1]*sm_w_mean1(*)^2, altitudes_LES/1000., thick=0.3, linestyle=2

; Gregory et al 2001 :
;oplot, A[0]*sm_buoyancy1_les(*)+A[1]*sm_full_d_rate_ude(*)*sm_w_mean1(*)^2, altitudes_LES/1000., thick=0.3
;oplot, B[0]*sm_buoyancy1_les(*)+B[1]*sm_full_d_rate_ude(*)*sm_w_mean1(*)^2, altitudes_LES/1000., thick=0.3
;oplot, C[0]*sm_buoyancy1_les(*)+C[1]*sm_full_d_rate_ude(*)*sm_w_mean1(*)^2, altitudes_LES/1000., thick=0.3, linestyle=2


PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

ENDFOR

print, ' ~~ comparing fit approaches :'

a1_simple=[1.848,1.723,1.455,1.280,1.993,1.789,1.582,1.050]
b1_simple=[-0.000842,-0.000511,-0.000251,-0.000329,-0.001473,-0.001081,-0.000129,-0.000453]
a1a_double=[1.820,1.716,1.454,1.274,1.918,1.767,1.568,1.029]
b1a_double=[0.000015,-0.000331,-0.000202,-0.000158, 0.001490,-0.000450, 0.000330, 0.000500]
a1b_double=[3.145,2.067,1.222,1.792,4.506,3.071,2.581,0.679]
b1b_double=[-0.001913,-0.000716,-0.000227,-0.000519,-0.005111,-0.001961,-0.000667,-0.000788]

a1_simple=[1.993,1.870,1.789,1.811,1.582,1.565,1.050,1.607,1.494,1.415,1.508,1.523,1.267,1.266,2.008,1.840,1.683,1.732,1.676,1.577,1.348,1.849,1.933,1.757,1.767,1.659,1.608,1.529,1.992,1.918,1.748,1.737,1.605,1.384,0.083,1.848,1.745,1.723,1.598,1.455,1.473,1.280]

b1_simple=[-0.001473,-0.001211,-0.001081,-0.000492,-0.000129,-0.000343,-0.000453,-0.001194,-0.000856,-0.000386,-0.000072,-0.000187,-0.000292,-0.000806,-0.001287,-0.000947,-0.000587,-0.000326,-0.000759,-0.000560, 0.000108,-0.001072,-0.001034,-0.000792,-0.000685,-0.000468,-0.000724,-0.000208,-0.000995,-0.000851,-0.000374,-0.000493,-0.000140,-0.000291,-0.000450,-0.000842,-0.000569,-0.000511,-0.000496,-0.000251,-0.000276,-0.000329]

print, 'simple approach'
print, 'mean a1 and mean b1'
print, MEAN(a1_simple),MEAN(b1_simple)
print, 'sigma'
print, STDDEV(a1_simple),STDDEV(b1_simple)

print, 'double approach'
print, 'mean a1a,b1a,a1b,b1b'
print, MEAN(a1a_double),MEAN(b1a_double),MEAN(a1b_double),MEAN(b1b_double)
print, 'sigma'
print, STDDEV(a1a_double),STDDEV(b1a_double),STDDEV(a1b_double),STDDEV(b1b_double)

;aa1=2.5 & bb1=0.0015
aa1=MEAN(a1_simple) & bb1=abs(MEAN(b1_simple))


print, '........ BUOYANCY AND VERTICAL VELOCITY ENTRAINMENT RATE DEPENDENCY'

B_w2_trac2 = make_array(nz)

dwdz_trac1 = deriv(altitudes_LES,smoothed_w_mean1_les)
;dwdz_trac2 = deriv(altitudes_LES,smoothed_w_mean2_les)
full_dwdz_trac1 = make_array(nz,nttot)
full_dadz_trac1 = make_array(nz,nttot)
FOR l=0, nttot -1 DO BEGIN
	full_dwdz_trac1(*,l) = deriv(altitudes_LES,w_mean1(*,l))
	full_dadz_trac1(*,l) = deriv(altitudes_LES,alpha1out(*,l))
ENDFOR
;alpha = 0.

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
;        if(fm_trac1_les(k,l) ne 0.) then full_e1(k,l)=(e_trac1_les(k,l)+e1_term2(k,l)+e1_term3(k,l))/fm_trac1_les(k,l) else full_e1(k,l)=0.
        if(w_mean1(k,l) ne 0.) then full_bw2(k,l)=grav*(tplume1moy(k,l)/tenv1moy_ude(k,l) -1.)/(w_mean1(k,l)^2) else full_bw2(k,l)=0.
;        if(w_mean1(k,l) ne 0.) then full_bw2(k,l)=grav*(tplume1moy(k,l)/tenv1moy(k,l) -1.)/(w_mean1(k,l)^2) else full_bw2(k,l)=0.

;        if(w_mean1(k,l) ne 0.) then full_bw2(k,l)=0.5*(alpha*buoyancy1_les(k,l)/(w_mean1(k,l)^2) - full_dwdz_trac1(k,l)/w_mean1(k,l)) else full_bw2(k,l)=0.
ENDFOR
ENDFOR

;~~~~~~~~~~~~~~~ fit

lt_plotindex_les3=fix(lt_plotindex_les3(0))
lt_plotindex_les0=fix(lt_plotindex_les0(0))

offset_fits=0
IF (TestCase eq 'Case_Z') THEN BEGIN
offset_fits =-10
ENDIF

nttotfit = lt_plotindex_les3+offset_fits - lt_plotindex_les0

Gamma_w2=make_array(nz,nttotfit,value=0.)

FOR ttt=0,nttotfit-1 DO BEGIN
        www = where(w_mean1(*,ttt+lt_plotindex_les0) ne 0.)
        if (www(0) eq -1 ) then print, 'AIE AIE AIE!'
	Gamma_w2(www,ttt)=aa1*buoyancy1_les(www,ttt+lt_plotindex_les0)/(w_mean1(www,ttt+lt_plotindex_les0))^2 - bb1
ENDFOR

D_out=make_array(2,nttotfit)
;D_out=make_array(3,nttotfit)

print, 'begining fits :'
print, lt_plotindex_les0, nttotfit, lt_plotindex_les3+offset_fits

FOR ttt=0,nttotfit-1 DO BEGIN
www = where(Gamma_w2(*,ttt) gt 0.)
eee = where(full_e1(*,ttt+lt_plotindex_les0) gt 0.)
nw = n_elements(www)
nee = n_elements(eee)
if (nw gt nee) then begin
Y = make_array(nee)
X = make_array(nee)
smooth_t,input=full_e1(eee,*),nz=nee,ndt=6,t0=ttt+lt_plotindex_les0,output=Y
Y = reform(full_e1(eee,ttt+lt_plotindex_les0))
X = reform(Gamma_w2(eee,ttt))
endif else begin
Y = make_array(nw)
X = make_array(nw)
smooth_t,input=full_e1(www,*),nz=nw,ndt=6,t0=ttt+lt_plotindex_les0,output=Y
Y = reform(full_e1(www,ttt+lt_plotindex_les0))
X = reform(Gamma_w2(www,ttt))
endelse
D = [0.08,0.6]
;D = [0.08,0.6,0.1]
coefs = lmfit(X,Y,D,/DOUBLE,function_name = 'myfunct2', itmax = 500)
;coefs = lmfit(X,Y,D,/DOUBLE,function_name = 'myfunct3', itmax = 500)
D_out(*,ttt)=D
ENDFOR

openw, lun, "fit_power_epsilon_AB_thermiques_"+TestCase+SubCase+les_special, /get_lun
printf, lun, "   Ae    ","    Be"
close, lun

openw, lun, "fit_power_epsilon_AB_thermiques_"+TestCase+SubCase+les_special, /append
for l=0, nttotfit-1 do printf, lun, D_out(0,l),D_out(1,l), format='((2x,F8.6)(4x,F8.6))'
FREE_LUN, lun
close, lun

print, 'fits complete, output in fit_power_epsilon_AB_thermiques_'+TestCase+SubCase+les_special

print, 'power approach'
print, 'mean Ae and mean Be'
print, MEAN(D_out(0,*)),MEAN(D_out(1,*))
print, 'sigma'
print, STDDEV(D_out(0,*)),STDDEV(D_out(1,*))

;what_I_plot = smoothed_e_rate_ude_trac1_les
what_I_plot = smoothed_e_rate_trac1_les
labels=['e_rate trac1']
title_user = TestCase+SubCase+LayerCase+' LES UDE entrainment rate dep with B/w2, average over '+taverage+' mn,' 
;filename = TestCase+SubCase+LayerCase+string(alpha,format='(F3.1)')+'Gcm_Les_Comp_e_Bw2.ps'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_e_Bw2.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=0., XMAX=0.2, YMIN=0., YMAX=0.1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, Y=what_I_plot, X=B_w2_trac1, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30, SYM=5, /NOLINES
AXES, XSTEP = 0.05 , XTITLE='B/w2', YSTEP=0.05, YTITLE='Entrainement rate m-1',NDECS=4

;oplot, smoothed_e_rate_trac2_les, B_w2_trac2, psym=5
FOR l=lt_plotindex_les0, lt_plotindex_les3-1 DO BEGIN
	oplot, full_bw2(*,l),full_e1(*,l),thick=0.05,psym=1
ENDFOR
;mean_full_e1 = make_array(nz) & mean_full_bw2 = make_array(nz)
;FOR k=0, nz-1 DO BEGIN
;	mean_full_e1(k) = MEAN(reform(full_e1(k,*)))
;	mean_full_bw2(k) = MEAN(reform(full_bw2(k,*)))
;ENDFOR
;oplot, mean_full_e1, mean_full_bw2, thick=0.3, psym = 2,color=5
;oplot, theoretical_e_trac1_les, B_w2_trac1,psym=2,thick=0.8,color=7
oplot,B_w2_trac1,(B_w2_trac1)/2.2222 + 0.0005,thick=0.3,color=7
;oplot, 0.0118*(B_w2_trac1/0.043)^(1./1.65),B_w2_trac1,thick=0.3,color=7
FOR l=lt_plotindex_les0, lt_plotindex_les3+offset_fits-1 DO BEGIN
;oplot, full_bw2(*,l),0.012*(full_bw2(*,l)/0.048)^(1./1.60),thick=0.1,color=7,psym=1

;oplot, full_bw2(*,l),0.04*(2.5*full_bw2(*,l))^(0.5)-0.0015,thick=0.1,color=7,psym=1
;oplot, full_bw2(*,l),0.045*(aa1*full_bw2(*,l)-bb1)^(0.6),thick=0.1,color=7,psym=1     ;entrainment formulation,baseline for a=2.5 & b=0.0015

;oplot, full_bw2(*,l),0.06*(aa1*full_bw2(*,l)-bb1)^(0.6),thick=0.1,color=7,psym=1
oplot, full_bw2(*,l),0.06*(aa1*full_bw2(*,l))^(0.75),thick=0.1,color=6,psym=1     ;detrainment formulation

oplot, full_bw2(*,l),MEAN(D_out(0,*))*(aa1*full_bw2(*,l)-bb1)^(MEAN(D_out(1,*))),thick=0.1,color=7,psym=1


ENDFOR
beta1=0.15

;FOR l=0, nttot-1 DO BEGIN
;;oplot, full_bw2(*,l),beta1*(2.5*full_bw2(*,l) - 0.0015)/(1.+beta1*(1.-w_mean1_env(*,l)/w_mean1(*,l))),thick=0.1,color=5,psym=1
;oplot, full_bw2(*,l),beta1*(2.5*full_bw2(*,l) - 0.0015)/(1.+beta1),thick=0.1,color=5,psym=1    ; earth formulation
;ENDFOR

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
full_dSiebesma = make_array(nz,nttot)


FOR k=0, nz-1 DO BEGIN
FOR l=0, nttot-1 DO BEGIN
        if(fm_trac1_les(k,l) ne 0.) then full_d1(k,l)=(d1_term1_ude(k,l)+d1_term2_ude(k,l)+d1_term3_ude(k,l))/fm_trac1_les(k,l) else full_d1(k,l)=-0.
;        if(fm_trac1_les(k,l) ne 0.) then full_d1(k,l)=(d1_term1(k,l)+d1_term2(k,l)+d1_term3(k,l))/fm_trac1_les(k,l) else full_d1(k,l)=-0.
        if(w_mean1(k,l) ne 0.) then full_dSiebesma(k,l)=0.75*0.5*buoyancy1_les(k,l)/(w_mean1(k,l)^2) -1.5*full_dwdz_trac1(k,l)/w_mean1(k,l) - full_dadz_trac1(k,l)/alpha1out(k,l) else full_dSiebesma(k,l)=-0.
ENDFOR
ENDFOR

;what_I_plot = smoothed_d_rate_ude_trac1_les
what_I_plot = smoothed_d_rate_trac1_les
labels=['d_rate trac1']
title_user = TestCase+SubCase+LayerCase+' LES UDE detrainment rate dep with B/w2, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_d_Bw2.ps'
;filename = TestCase+SubCase+LayerCase+string(alpha,format='(F3.1)')+'Gcm_Les_Comp_d_Bw2.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.1, XMAX=0.1, YMIN=0., YMAX=0.1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, Y=what_I_plot, X=full_bw2(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30, SYM=5, /NOLINES
AXES, XSTEP = 0.05 , XTITLE='B/w2', YSTEP=0.05, YTITLE='Detrainment rate',NDECS=4

FOR l=0, nttot-1 DO BEGIN
        oplot, full_bw2(*,l),full_d1(*,l),thick=0.05,psym=1
ENDFOR
;oplot, theoretical_d_trac1_les, full_bw2(*,lt_plotindex_les),psym=2,thick=0.8,color=7
;oplot,B_w2_trac1/2.7 + 0.0002,B_w2_trac1,thick=0.3,color=7
oplot,B_w2_trac1,B_w2_trac1/2.222 + 0.0002,thick=0.3,color=7
FOR l=0, nttot-1 DO BEGIN
oplot, full_bw2(*,l),0.06*(2.5*full_bw2(*,l))^(0.75),thick=0.1,color=7,psym=1   ;detrainment formulation
;oplot, full_bw2(*,l),0.045*(2.5*full_bw2(*,l)-0.0015)^(0.6),thick=0.1,color=6,psym=1   ;entrainment formulation
ENDFOR
FOR l=0, nttot-1 DO BEGIN
oplot, full_bw2(*,l),0.06*(-2.5*full_bw2(*,l))^(0.75),thick=0.1,color=7,psym=1   ;detrainment formulation
ENDFOR
;a1=3.

;FOR l=0, nttot-1 DO BEGIN
;oplot, full_bw2(*,l),a1*(beta1/(1.+beta1))*full_bw2(*,l),thick=0.1,color=5,psym=1   ;earth formulation
;ENDFOR

;oplot, 0.0118*(B_w2_trac1/0.043)^(1./1.65),B_w2_trac1,thick=0.3,color=6
oplot, B_w2_trac1,0.0105*(B_w2_trac1/0.048)^(1./1.7),thick=0.3,color=6

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

print, 'Detrainment : new approach '

full_d1_v2 = make_array(nz,nttot)
df_dz_les_full = make_array(nz,nttot)
da_dt_les_full = make_array(nz,nttot)

FOR k=0, nz-1 DO BEGIN
        da_dt_les_full(k,*) = deriv(localtime,reform(alpha1out(k,*)))/3700.
ENDFOR

FOR l=0, nttot-1 DO BEGIN
        df_dz_les_full(*,l) = deriv(altitudes_LES,reform(fm_trac1_les(*,l)))
	FOR k=0, nz-1 DO BEGIN
        	if(fm_trac1_les(k,l) ne 0.) then full_d1_v2(k,l)=full_e1(k,l) - df_dz_les_full(k,l)/fm_trac1_les(k,l) -(rho(k,l)/fm_trac1_les(k,l))*da_dt_les_full(k,l) else full_d1_v2(k,l)=0.
;                if(fm_trac1_les(k,l) ne 0.) then full_d1_v2(k,l)=full_e1(k,l) - df_dz_les_full(k,l)/fm_trac1_les(k,l) else full_d1_v2(k,l)=0.
	ENDFOR
ENDFOR


;; ~~~~~~~~~~~~~~~~~~~~ DETRAINMENT FIT

;lt_plotindex_les3=fix(lt_plotindex_les3(0))
;lt_plotindex_les0=fix(lt_plotindex_les0(0))

offset_fits=0
IF (TestCase eq 'Case_Z') THEN BEGIN
offset_fits =-10
ENDIF

nttotfit = lt_plotindex_les3+offset_fits - lt_plotindex_les0

B_w2_fits=make_array(nz,nttotfit,value=0.)

FOR ttt=0,nttotfit-1 DO BEGIN
        www = where(w_mean1(*,ttt+lt_plotindex_les0) ne 0.)
        if (www(0) eq -1 ) then print, 'AIE AIE AIE!'
;        B_w2_fits(www,ttt)=aa1*buoyancy1_les(www,ttt+lt_plotindex_les0)/(w_mean1(www,ttt+lt_plotindex_les0))^2
         B_w2_fits(www,ttt)=buoyancy1_les(www,ttt+lt_plotindex_les0)/(w_mean1(www,ttt+lt_plotindex_les0))^2
ENDFOR

E_out=make_array(2,nttotfit)
print, 'detrainment: begining fits :'
print, lt_plotindex_les0, nttotfit, lt_plotindex_les3+offset_fits

FOR ttt=0,nttotfit-1 DO BEGIN
eee = where((full_d1_v2(*,ttt+lt_plotindex_les0) gt 0.) and (abs(B_w2_fits(*,ttt)) gt 0.001))
nee = n_elements(eee)

eee=eee(4:nee-1)
nee = n_elements(eee)

Y = make_array(nee)
X = make_array(nee)

smooth_t,input=full_d1_v2(eee,*),nz=nee,ndt=6,t0=ttt+lt_plotindex_les0,output=Y

;Y = reform(full_d1_v2(eee,ttt+lt_plotindex_les0))
X = reform(B_w2_fits(eee,ttt))

E = [-0.38,0.0005]
coefs = lmfit(X,Y,E,/DOUBLE,function_name = 'myfunct', itmax = 1000)
E_out(*,ttt)=E
ENDFOR

openw, lun, "fit_lin_delta_AB_thermiques_"+TestCase+SubCase+les_special, /get_lun
printf, lun, "   Ad    ","    Bd"
close, lun

openw, lun, "fit_lin_delta_AB_thermiques_"+TestCase+SubCase+les_special, /append
for l=0, nttotfit-1 do printf, lun, E_out(0,l),E_out(1,l), format='((2x,F9.6)(4x,F9.6))'
FREE_LUN, lun
close, lun

print, 'fits complete, output in fit_lin_delta_AB_thermiques_'+TestCase+SubCase+les_special

print, 'delta: lin approach'
print, 'mean Ad and mean Bd'
print, MEAN(E_out(0,*)),MEAN(E_out(1,*))
print, 'sigma'
print, STDDEV(E_out(0,*)),STDDEV(E_out(1,*))

;~~~~~~~~~~~


;what_I_plot = smoothed_d_rate_ude_trac1_les
what_I_plot = full_d1_v2(*,lt_plotindex_les)
labels=['d_rate trac1']
title_user = TestCase+SubCase+LayerCase+' LES UDE detrainment rate dep with B/w2, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_d2_Bw2.ps'
;filename = TestCase+SubCase+LayerCase+string(alpha,format='(F3.1)')+'Gcm_Les_Comp_d_Bw2.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.1, XMAX=0.1, YMIN=0., YMAX=0.1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, Y=what_I_plot, X=full_bw2(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30, SYM=5, /NOLINES
AXES, XSTEP = 0.01 , XTITLE='B/w2', YSTEP=0.01, YTITLE='Detrainment rate',NDECS=2

FOR l=lt_plotindex_les0, lt_plotindex_les3+offset_fits-1 DO BEGIN
        oplot, full_bw2(*,l),full_d1_v2(*,l),thick=0.05,psym=1
ENDFOR
;oplot, theoretical_d_trac1_les, full_bw2(*,lt_plotindex_les),psym=2,thick=0.8,color=7
;oplot,B_w2_trac1/2.7 + 0.0002,B_w2_trac1,thick=0.3,color=7
oplot,B_w2_trac1,B_w2_trac1/2.222 + 0.0002,thick=0.3,color=7
FOR l=lt_plotindex_les0, lt_plotindex_les3+offset_fits-1 DO BEGIN

oplot, full_bw2(*,l), MEAN(E_out(0,*))*full_bw2(*,l)+MEAN(E_out(1,*)),thick=0.1,color=7,psym=1    ;new detrainment formulation

oplot, full_bw2(*,l), -0.45*full_bw2(*,l)+0.0005,thick=0.1,color=8,psym=1    ;new detrainment formulation

;oplot, full_bw2(*,l),0.06*(2.5*full_bw2(*,l))^(0.75),thick=0.1,color=7,psym=1   ;detrainment formulation (classical)
;oplot, full_bw2(*,l),0.045*(2.5*full_bw2(*,l)-0.0015)^(0.6),thick=0.1,color=6,psym=1   ;entrainment formulation
ENDFOR
;FOR l=0, nttot-1 DO BEGIN
;oplot, full_bw2(*,l),0.06*(-2.5*full_bw2(*,l))^(0.75),thick=0.1,color=7,psym=1   ;detrainment formulation
;ENDFOR
;a1=3.

;FOR l=0, nttot-1 DO BEGIN
;oplot, full_bw2(*,l),a1*(beta1/(1.+beta1))*full_bw2(*,l),thick=0.1,color=5,psym=1   ;earth formulation
;ENDFOR

;oplot, 0.0118*(B_w2_trac1/0.043)^(1./1.65),B_w2_trac1,thick=0.3,color=6
oplot, B_w2_trac1,0.0105*(B_w2_trac1/0.048)^(1./1.7),thick=0.3,color=6

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename



;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ mass flux derivative fits


;offset_fits=0
;IF (TestCase eq 'Case_Z') THEN BEGIN
;offset_fits =-10
;ENDIF

;nttotfit = lt_plotindex_les3+offset_fits - lt_plotindex_les0

full_d1_fits=make_array(nz,nttotfit,value=0.)


full_1_f_dfdz = make_array(nz,nttot)

FOR l=0, nttot-1 DO BEGIN
FOR k=0,nz-1 DO BEGIN
       IF (fm_trac1_les(k,l) GT 0.) THEN full_1_f_dfdz(k,l)=df_dz_les_full(k,l)/fm_trac1_les(k,l) ELSE full_1_f_dfdz(k,l)=0.
ENDFOR
ENDFOR	

e_from_fit = make_array(nz,nttotfit)
FOR l=0, nttotfit-1 DO BEGIN
FOR k=0,nz-1 DO BEGIN
	e_from_fit(k,l) = MEAN(D_out(0,*))*(MAX([0.,aa1*full_bw2(k,l+lt_plotindex_les0)-bb1]))^(MEAN(D_out(1,*)))
ENDFOR
ENDFOR


FOR ttt=0,nttotfit-1 DO BEGIN
        www = where(e_from_fit(*,ttt) gt full_1_f_dfdz(*,ttt+lt_plotindex_les0) )
        if (www(0) eq -1 ) then print, 'AIE AIE AIE!'
        full_d1_fits(www,ttt)=e_from_fit(www,ttt) - full_1_f_dfdz(www,ttt+lt_plotindex_les0)
ENDFOR

F_out=make_array(3,nttotfit)
print, 'mass flux derivative: begining fits :'
print, lt_plotindex_les0, nttotfit, lt_plotindex_les3+offset_fits

FOR ttt=0,nttotfit-1 DO BEGIN
;FOR ttt=3,nttotfit-4 DO BEGIN

eee = where((full_d1_fits(*,ttt) gt 0.) and (altitudes_LES(*) gt 500.) and (B_w2_fits(*,ttt) gt 0.))
fff = where((full_d1_fits(*,ttt) gt 0.) and (altitudes_LES(*) gt 500.) and (B_w2_fits(*,ttt) lt 0.))
nee = n_elements(eee)
nff = n_elements(fff)

;eee=eee(floor(nee/4.):nee-1)
;nee = n_elements(eee)

Y1 = make_array(nee)
X1 = make_array(nee)

Y2 = make_array(nff)
X2 = make_array(nff)

;smooth_t,input=full_d1_v2(eee,*),nz=nee,ndt=6,t0=ttt+lt_plotindex_les0,output=Y

Y1 = reform(full_d1_fits(eee,ttt))
X1 = reform(B_w2_fits(eee,ttt))


Y2 = reform(full_d1_fits(fff,ttt))
X2 = reform(B_w2_fits(fff,ttt))

;smooth_t,input=full_d1_fits(eee,*),nz=nee,ndt=4,t0=ttt,output=Y
;smooth_t,input=B_w2_fits(eee,*),nz=nee,ndt=4,t0=ttt,output=X

;F = [-0.38,0.0001]
F = -0.5
coefs = lmfit(X2,Y2,F,/DOUBLE,function_name = 'myfunct4', itmax = 1000)

F_out(0,ttt)=MEAN(reform(Y1))
;F_out(1:2,ttt)=F
F_out(1,ttt)=F

ENDFOR

;openw, lun, "fit_lin_delta_on_f_AB_thermiques_"+TestCase+SubCase+les_special, /get_lun
;printf, lun, "   Ada    ","    Ad    ","    Bd"
;close, lun
;
;openw, lun, "fit_lin_delta_on_f_AB_thermiques_"+TestCase+SubCase+les_special, /append
;for l=0, nttotfit-1 do printf, lun, F_out(0,l),F_out(1,l),F_out(2,l), format='((2x,F9.6)(4x,F9.6)(4x,F9.6))'
;FREE_LUN, lun
;close, lun

openw, lun, "fit_lin_delta_on_f_AB_thermiques_"+TestCase+SubCase+les_special, /get_lun
printf, lun, "   Ad    ","    Bd"
close, lun

openw, lun, "fit_lin_delta_on_f_AB_thermiques_"+TestCase+SubCase+les_special, /append
for l=0, nttotfit-1 do printf, lun, F_out(0,l),F_out(1,l), format='((2x,F9.6)(4x,F9.6))'
FREE_LUN, lun
close, lun


print, 'delta:lin approach'
print, 'mean Ada and mean Ad, Bd'
print, MEAN(F_out(0,*)),MEAN(F_out(1,*)) ;,MEAN(F_out(2,*))
print, 'sigma'
print, STDDEV(F_out(0,*)),STDDEV(F_out(1,*)) ;,STDDEV(F_out(2,*))

what_I_plot = full_1_f_dfdz(*,lt_plotindex_les)
labels=['TH mass flux vertical derivative']
title_user = TestCase+SubCase+LayerCase+' mass flux vertical derivative comparison, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_1_f_dfdz.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.02, XMAX=0.02, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(1)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.005 , XTITLE='m-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=6

df_dz_from_param = make_array(nz)
;FOR k=0, nz-1 DO BEGIN
;	df_dz_from_param(k) = MEAN(D_out(0,*))*(MAX([0.,aa1*full_bw2(k,lt_plotindex_les)-bb1]))^(MEAN(D_out(1,*))) - MEAN(F_out(0,*))*full_bw2(k,lt_plotindex_les)+MEAN(F_out(1,*))
;ENDFOR
FOR k=0, nz-1 DO BEGIN
       IF (full_bw2(k,lt_plotindex_les) gt 0.) THEN BEGIN
       df_dz_from_param(k) = MEAN(D_out(0,*))*(MAX([0.,aa1*full_bw2(k,lt_plotindex_les)-bb1]))^(MEAN(D_out(1,*))) - MEAN(F_out(0,*))
       ENDIF ELSE BEGIN
;       df_dz_from_param(k) = MEAN(D_out(0,*))*(MAX([0.,aa1*full_bw2(k,lt_plotindex_les)-bb1]))^(MEAN(D_out(1,*))) - (MEAN(F_out(1,*))*full_bw2(k,lt_plotindex_les)+MEAN(F_out(2,*)))
       df_dz_from_param(k) = MEAN(D_out(0,*))*(MAX([0.,aa1*full_bw2(k,lt_plotindex_les)-bb1]))^(MEAN(D_out(1,*))) - (MEAN(F_out(1,*))*full_bw2(k,lt_plotindex_les))
       ENDELSE
ENDFOR

oplot, df_dz_from_param, altitudes_LES/1000., psym=5
;oplot, df_dz_les1/smoothed_fm_trac1_les, altitudes_LES/1000., psym=4
;oplot, df_dz_les2, altitudes_LES/1000., psym=5
;oplot, df_dz_param, altitudes_LES/1000., psym=6, color=5
;oplot, df_dz_param2, altitudes_LES/1000., psym=6, color=6
;print, "fm*(e-d)"
;print, smoothed_fm_trac1_les*(0.045*(2.5*B_w2_trac1-0.0015)^(0.6) - 0.06*(-2.5*B_w2_trac1)^(0.75))
;print, smoothed_fm_trac1_les
;print, B_w2_trac1

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename




print, '........ BUOYANCY AND VERTICAL VELOCITY VERTICAL VELOCITY DERIVATIVE DEPENDANCY'

full_dw2_w2 = make_array(nz,nttot)
w2_mean1 = make_array(nz,nttot)
dw2_dz = make_array(nz,nttot)
w2_param = make_array(nz,nttot)
FOR l=0, nttot-1 DO BEGIN
	w2_mean1(*,l) = w_mean1(*,l)^2.
	dw2_dz(*,l) = deriv(altitudes_LES,w2_mean1(*,l))
ENDFOR
FOR k=0, nz-1 DO BEGIN
FOR l=0, nttot-1 DO BEGIN
       if(w2_mean1(k,l) ne 0.) then begin 
;       if(full_bw2(k,l) ge 0.) then full_dw2_w2(k,l) = 0.5*dw2_dz(k,l)/w2_mean1(k,l) + (1.-w_mean1_env_ude(k,l)/w_mean1(k,l))*0.012*(full_bw2(k,l)/0.048)^(1./1.6) else full_dw2_w2(k,l)=0.5*dw2_dz(k,l)/w2_mean1(k,l)
;       if(full_bw2(k,l) ge 0.) then full_dw2_w2(k,l) = 0.5*dw2_dz(k,l)/w2_mean1(k,l) + (1.-w_mean1_env(k,l)/w_mean1(k,l))*0.009*(full_bw2(k,l)/0.048)^(1./1.9) else full_dw2_w2(k,l)=0.5*dw2_dz(k,l)/w2_mean1(k,l)
        full_dw2_w2(k,l) = 0.5*dw2_dz(k,l)/w2_mean1(k,l); + (1.-w_mean1_env_ude(k,l)/w_mean1(k,l))*full_e1(k,l)
       endif else begin 
              full_dw2_w2(k,l)= 0.
       endelse

       if(w2_mean1(k,l) ne 0.) then begin
;       if(full_bw2(k,l) ge 0.) then w2_param(k,l) = full_bw2(k,l)/6. - (0.014*(full_bw2(k,l)/0.05)^(1./1.35))/2. else w2_param(k,l) = full_bw2(k,l)/6. - full_d1(k,l)/2. 
        w2_param(k,l) = 2.5*full_bw2(k,l) - 0.0015
       endif else begin
       w2_param(k,l)=0.
       endelse
;       if(w2_mean1(k,l) ne 0.) then w2_param(k,l) = full_bw2(k,l) else w2_param(k,l)=0.
ENDFOR
ENDFOR

what_I_plot = full_dw2_w2(*,lt_plotindex_les)

labels=['0.5*(dw2/dz)/w2']
title_user = TestCase+SubCase+LayerCase+' LES w2 equation vs B/w2, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_w2_Bw2.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, YMIN=-0.1, YMAX=0.1, XMIN=-0.1, XMAX=0.1, TITLE=title_user
cols=INDGEN(1)+2
GPLOT,Y=what_I_plot, X=w2_param(*,lt_plotindex_les), /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30, SYM=5, /NOLINES
AXES, XSTEP = 0.1 , XTITLE='aB/w2 - b', YSTEP=0.1, YTITLE='0.5*(dw2/dz)/w2 + E/fm',NDECS=4

FOR l=0, nttot-1 DO BEGIN
        oplot, w2_param(*,l), full_dw2_w2(*,l),thick=0.05,psym=1
ENDFOR
oplot, w2_param(*,lt_plotindex_les),w2_param(*,lt_plotindex_les),thick=0.05,color=7

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
		teta_est(k) = (f_star(k)*teta_p(k-1)+a_star(k)*0.25*(teta_les(k,lt_plotindex_les) + teta_p(k-1)))/(f_star(k) + a_star(k))
	endif else begin
		teta_est(k) = teta_p(k-1)
	endelse
	zbuoy_est(k) = grav*(teta_est(k)/teta_les(k,lt_plotindex_les) -1.)
	zw2fact=fact_epsilon*2.*approx_zdz_les(k)/(1.+betalpha)
	zdw2=afact*zbuoy_est(k)/fact_epsilon
;        w_est(k+1) = MAX([0.0001,exp(-zw2fact)*(w_est(k)-zdw2)+zdw2])
        w_est(k+1) = MAX([0.0001,w_est(k)*(1.-2.*approx_zdz_les(k)*0.5)+2.*1.*approx_zdz_les(k)*zbuoy_est(k)])
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
;        zw2(k+1) = MAX([0.0001,exp(-zw2fact)*(zw2(k)-zdw2)+zdw2])
        zw2(k+1) = MAX([0.0001,zw2(k)*(1.-2.*approx_zdz_les(k)*0.5)+2.*1.*approx_zdz_les(k)*zbuoy_est(k)])
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

smoothed_rho = make_array(nz)
FOR t=-ns,ns DO BEGIN
        smoothed_rho = smoothed_rho + REFORM(rho(*,lt_plotindex_les+t))
ENDFOR
smoothed_rho = smoothed_rho/nstot

;continuity1 = smoothed_rho*smoothed_da_dt + df_dz_les1 - smoothed_e_rate_ude_trac1_les*smoothed_fm_trac1_les + smoothed_d_rate_ude_trac1_les*smoothed_fm_trac1_les
continuity1 = smoothed_rho*smoothed_da_dt + df_dz_les1 - smoothed_e_rate_trac1_les*smoothed_fm_trac1_les + smoothed_d_rate_trac1_les*smoothed_fm_trac1_les

print, '........ CONTINUITY CHECK'

;what_I_plot = [[continuity1],[smoothed_rho*smoothed_da_dt],[df_dz_les1],[-smoothed_e_rate_ude_trac1_les*smoothed_fm_trac1_les],[smoothed_d_rate_ude_trac1_les*smoothed_fm_trac1_les]]
what_I_plot = [[continuity1],[smoothed_rho*smoothed_da_dt],[df_dz_les1],[-smoothed_e_rate_trac1_les*smoothed_fm_trac1_les],[smoothed_d_rate_trac1_les*smoothed_fm_trac1_les]]
labels=['total continuity','rho*da/dt','df/dz','-E','D']
title_user = TestCase+SubCase+LayerCase+' LES UDE continuity check, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_continuity.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-0.0001, XMAX=0.0001, YMIN=0, YMAX=10, TITLE=title_user
cols=INDGEN(5)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30 
AXES, XSTEP = 0.00001 , XTITLE='kg.m-2.s-1', YSTEP=1, YTITLE='Altitude (km)',NDECS=4
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

print, '........ UDE TRANSPORT : Q'

what_I_plot = [[wq(*,lt_plotindex_les)],[wq_updraft(*,lt_plotindex_les)],[wq_downdraft(*,lt_plotindex_les)]]

labels=['total turbulent transport','within updraft turbulent transport','within downdraft turbulent transport']
title_user = TestCase+SubCase+LayerCase+' LES UDE turbulence Q transport, average over '+taverage+' mn'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Qturbu_ude.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-1, XMAX=1.5, YMIN=0, YMAX=6, TITLE=title_user
cols=INDGEN(3)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=9, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 0.25 , XTITLE='m.K/s', YSTEP=0.5, YTITLE='Altitude (km)',NDECS=3
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

print, '........ BL TOP'

dteta_dz_les_full=make_array(nz,nttot)

FOR l=0, nttot-1 DO BEGIN
	dteta_dz_les_full(*,l)= deriv(altitudes_LES,reform(teta_les(*,l)))
ENDFOR

altindex=make_array(nttot,value=0)
altboundary=make_array(nttot)

altindex2=make_array(nttot,value=0)
altboundary2=make_array(nttot)

print, nttot
FOR l=0, nttot-1 DO BEGIN
 
        if (n_elements(where(wt(*,l) eq min(wt(*,l)))) eq 1) then begin
	altindex(l)=where(wt(*,l) eq min(wt(*,l)))
        endif
        if (altindex(l) ne -1) then altboundary(l)=altitudes_LES(altindex(l))/1000.
        if (l ne 0) then begin
            if ((altboundary(l) gt 1.5*altboundary(l-1)) and (localtime(l) gt 10.)) then altboundary(l)=altboundary(l-1)
        endif

        if (n_elements(where(dteta_dz_les_full(*,l) eq max(dteta_dz_les_full(*,l)))) eq 1) then begin
        altindex2(l)=where(dteta_dz_les_full(*,l) eq max(dteta_dz_les_full(*,l)))
        endif
        if (altindex2(l) ne -1) then altboundary2(l)=altitudes_LES(altindex2(l))/1000.

ENDFOR
FOR l=1, nttot-2 DO BEGIN
     if (localtime(l) gt 10.) then begin
     if ((altboundary2(l) gt 1.5*altboundary2(l-1)) and (altboundary2(l+1) lt 1.5*altboundary2(l-1))) then altboundary2(l)=(altboundary2(l-1)+altboundary2(l+1))/2.
     if (altboundary2(l) gt 1.5*altboundary2(l-1)) then altboundary2(l)=altboundary2(l-1)
     endif
     if (localtime(l) gt 10.) then begin
     if ((altindex2(l) gt 1.5*altindex2(l-1)) and (altindex2(l+1) lt 1.5*altindex2(l-1))) then altindex2(l)=(altindex2(l-1)+altindex2(l+1))/2.
     if (altindex2(l) gt 1.5*altindex2(l-1)) then altindex2(l)=altindex2(l-1)
     endif

ENDFOR


altboundary2 = SMOOTH(altboundary2,5,/EDGE_TRUNCATE)
altboundary2 = SMOOTH(altboundary2,10,/EDGE_TRUNCATE)
altindex2 = SMOOTH(altindex2,5,/EDGE_TRUNCATE)
altindex2 = SMOOTH(altindex2,10,/EDGE_TRUNCATE)

FOR l=0, nttot-1 DO BEGIN
     if (localtime(l) gt 16.) then altboundary2(l)=0.
     if ((localtime(l) lt 9.) or (localtime(l) gt 16.)) then altindex2(l)=0.
ENDFOR

what_I_plot = [[altboundary],[altboundary2]]

labels=['LES Boundary Layer top from Heat flux minimum','LES Boundary Layer top grom Teta gradient']

title_user = TestCase+SubCase+LayerCase+' LES boundary layer top'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Zi.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=0, YMAX=9, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time LES', YSTEP=1, YTITLE='Altitude (km)',NDECS=1

;oplot,localtime_gcm,zi_gcm/1000.,thick=0.4,color=5
oplot,localtime_gcm,altitudes_GCM(lmax_gcm)/1000.,thick=0.4,color=6


PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

what_I_plot = [[altindex],[altindex2]]
labels=['LES Boundary Layer top from Heat flux minimum','LES Boundary Layer top grom Teta gradient']
title_user = TestCase+SubCase+LayerCase+' LES boundary layer top'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_Ziindex.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=6, XMAX=18, YMIN=0, YMAX=400, TITLE=title_user
cols=INDGEN(2)+2
GPLOT, X=localtime, Y=what_I_plot, /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='local time LES', YSTEP=10, YTITLE='Altitude (km)',NDECS=1
PSCLOSE, /NOVIEW
spawn, 'ps2png '+filename

if (1 eq 1) then begin

openw, lun, "input_zipbl", /get_lun
for l=0, nttot-1 do printf, lun, altindex2(l),localtime(l), format='((2x,I0)(4x,F8.2))'
FREE_LUN, lun
close, lun

endif

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

print, '........ TRACER FLUX'

if (got_tracer_flux eq 'true') then begin

what_I_plot = [[alpha1out(*,lt_plotindex_les)*wq_updraft(*,lt_plotindex_les)],[beta1out(*,lt_plotindex_les)*wq_downdraft(*,lt_plotindex_les)],[(1.-alpha1out(*,lt_plotindex_les)-beta1out(*,lt_plotindex_les))*wq_env_ude(*,lt_plotindex_les)],[alpha1out(*,lt_plotindex_les)*w_mean1(*,lt_plotindex_les)*(q_mean_up(*,lt_plotindex_les) - q_mean(*,lt_plotindex_les))+beta1out(*,lt_plotindex_les)*w_mean1_down(*,lt_plotindex_les)*(q_mean_down(*,lt_plotindex_les) - q_mean(*,lt_plotindex_les))],[wq(*,lt_plotindex_les)]]

;labels=['TH updraft heat flux','Mellor and Yamada gcm heat flux','Total']
labels=['LES updraft tracer flux','LES downdraft tracer flux','LES env tracer flux','LES organized structures tracer flux','Total']

title_user = TestCase+SubCase+LayerCase+' LES vertical tracer flux'
filename = TestCase+SubCase+LayerCase+'Gcm_Les_Comp_WQ.ps'
PSOPEN, THICK=200, CHARSIZE=120, FILE = filename, FONT = 5, TFONT = 5
CS, SCALE=28
GSET, XMIN=-1, XMAX=4, YMIN=0, YMAX=5, TITLE=title_user
;cols=INDGEN(3)+2
cols=INDGEN(5)+2
GPLOT, X=what_I_plot, Y=altitudes_LES/1000., /LEGEND, LEGPOS=1, COL=cols, LABELS=labels, THICK = 30
AXES, XSTEP = 1, XTITLE='vertical turbulent tracer flux', YSTEP=1, YTITLE='Altitude (km)',NDECS=1

oplot, alpha1out(*,lt_plotindex_les)*w_mean1(*,lt_plotindex_les)*(q_mean_up(*,lt_plotindex_les) - q_mean(*,lt_plotindex_les)), altitudes_LES/1000., color=2,psym=2,thick=0.5
oplot, beta1out(*,lt_plotindex_les)*w_mean1_down(*,lt_plotindex_les)*(q_mean_down(*,lt_plotindex_les) - q_mean(*,lt_plotindex_les)), altitudes_LES/1000., color=3,psym=2,thick=0.5


;oplot, alpha1out(*,lt_plotindex_les)*(w_mean1(*,lt_plotindex_les) - w_mean1_full(*,lt_plotindex_les))*(tplume1moy(*,lt_plotindex_les) - tmoy_full(*,lt_plotindex_les)), altitudes_LES/1000., color=2

;oplot, beta1out(*,lt_plotindex_les)*(w_mean1_down(*,lt_plotindex_les) - w_mean1_full(*,lt_plotindex_les))*(tdown1moy(*,lt_plotindex_les) - tmoy_full(*,lt_plotindex_les)), altitudes_LES/1000.,color=6

;oplot, smoothed_hf1_ude_term4, altitudes_LES/1000., color = 5

;oplot, smoothed_hf1_term1, altitudes_LES/1000.,thick=0.1,LINESTYLE = 5
;oplot, smoothed_hf1_term2, altitudes_LES/1000.,color=2,thick=0.1,LINESTYLE = 5
;oplot, smoothed_hf1_term3, altitudes_LES/1000.,color=8,thick=0.1,LINESTYLE = 5

;oplot, wt(*,lt_plotindex_les), altitudes_LES/1000.

;oplot, smoothed_hf1_term3(*,lt_plotindex_les), altitudes_LES/1000., psym=7

PSCLOSE, /NOVIEW

spawn, 'ps2png '+filename

endif

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
;oplot, (altitudes_rel_LES(*,lt_plotindex_les)-0.075)/187.931 + 0.9977, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
oplot, (altitudes_rel_LES(*,lt_plotindex_les))/19.231 + 0.9938, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
oplot, (altitudes_rel_LES(*,lt_plotindex_les)-0.075)/187.931 + 0.9982, altitudes_rel_LES(*,lt_plotindex_les),thick=0.3,color=7
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
