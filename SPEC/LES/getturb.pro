pro getturb, saveps=saveps, noplot=noplot

;----------------------------------
; USE: getturb
;      getturb, saveps='false'     
;----------------------------------

history_interval_s = 400.
history_interval_s = 100.
smoothampl=3700/history_interval_s
;smoothampl=0.  ;; no smooth
;smoothampl=10.
smoothampl=2.
;smoothampl=5.

;;
;; constantes
;;
p0=610. & t0=220. & r_cp=1.0/4.4 & grav=3.72 & R=192.
;p0=610. & t0=220. & r_cp=1.0/3.9 & grav = 3.72 & R=191.
;print, 'ATTENTION ATTENTION R/cp !!!!', r_cp 

; INTERCOMP INTERCOMP
r_cp = 192./770.

;
; graphics definition
;
if (n_elements(saveps) eq 0) then saveps='true'
if (n_elements(noplot) eq 0) then noplot='false'
if (saveps eq 'false') then begin
   ;!p.multi=[0,3,2] 
   !P.CHARSIZE=2.
   WINDOW, /PIXMAP & WDELETE & DEVICE,BYPASS_TRANSLATION=0,DECOMPOSED=0,RETAIN=2
endif else begin
   ;PREF_SET, 'IDL_PATH', '/home/spiga/Save/SOURCES/IDL/fsc_psconfig:<IDL_DEFAULT>', /COMMIT
   ;PREF_SET, 'IDL_PATH', '/home/spiga/SVN/trunk/mesoscale/PLOT/MINIMAL:<IDL_DEFAULT>', /COMMIT
   PREF_SET, 'IDL_PATH', '/home/spiga/MODELES/MESOSCALE_DEV/PLOT/MINIMAL:<IDL_DEFAULT>', /COMMIT
endelse

;
; retrieve fields
;
openr,unit,'getturb.dat',/get_lun,error=err
IF (err ne 0) THEN BEGIN

;
; input files
;
OPENR, 22, 'input_coord' & READF, 22, lonu & READF, 22, latu & READF, 22, lsu & READF, 22, lctu & CLOSE, 22
OPENR, 23, 'input_more' & READF, 23, hgtu, tsurfu & CLOSE, 23

;
; get fields
;
domain='d01' & filesWRF = FindFile('wrfout_'+domain+'_????-??-??_??:??:??') & nf=n_elements(filesWRF)


; get dimensions
;
id=ncdf_open(filesWRF(0))
NCDF_DIMINQ, id, NCDF_DIMID(id, 'west_east'    ), toto, nx & NCDF_DIMINQ, id, NCDF_DIMID(id, 'south_north'  ), toto, ny
NCDF_DIMINQ, id, NCDF_DIMID(id, 'bottom_top'   ), toto, nz & NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), toto, nt
NCDF_CLOSE, id 
id=ncdf_open(filesWRF(nf-1))  ;; for interrupted runs
NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), toto, ntlast
NCDF_CLOSE, id

;
; prepare loop
;
;nloop1 = nf-1 & nloop2 = nt & yeye = 0 
yeye = 0 & nttot = (nf-1)*nt + ntlast

localtime = lctu + history_interval_s*findgen(nttot)/3700.

wt  = fltarr(nz,nttot) & tke = fltarr(nz,nttot) & ztke = fltarr(nz,nttot) & t = fltarr(nz,nttot)
p = fltarr(nz) & ph = fltarr(nz) & pht = fltarr(nz,nttot) & pt = fltarr(nz,nttot) & stst = fltarr(nz,nttot)
xtke = fltarr(nz,nttot) & ytke = fltarr(nz,nttot) & wmax = fltarr(nz,nttot) & wmin = fltarr(nz,nttot)
depressions = fltarr(nttot) & psmin = fltarr(nttot)

;
; loop loop
;
for loop  = 0, nf-1 do begin
                                          timetime = SYSTIME(1)
  if (loop ne nf-1) then nloop2=nt else nloop2=ntlast

for loop2 = 0, nloop2-1 do begin

   anomalt = 1. & anomalu = 1. & anomalv = 1.
	;print, loop, loop2

   ; ------------
   ; t' = t - <t>
   ; ------------
 tprime = getget(filesWRF(loop), 'T', anomaly=anomalt, count=[0,0,0,1], offset=[0,0,0,loop2])  ;; t' = t - <t>
 t(*,yeye) = t0 + TEMPORARY(anomalt)
   ; ------
   ; w' = w   
   ; ------
 wprime = getget(filesWRF(loop), 'W', count=[0,0,0,1], offset=[0,0,0,loop2])     
 	for toto = 0, nz-1 do begin
 		wmax(toto,yeye) = max(reform(wprime(*,*,toto,0)),min=tutu) & wmin(toto,yeye) = tutu 
 	endfor
;   ; ------
;   ; u' = u and v' = v   (car PAS de background wind !)
;   ; ------
; uprime = getget(filesWRF(loop), 'U', count=[0,0,0,1], offset=[0,0,0,loop2])
; vprime = getget(filesWRF(loop), 'V', count=[0,0,0,1], offset=[0,0,0,loop2])
   ; --------------------------------------------------------
   ; tke = 0.5 ( <u'^2> + <v'^2> + <w'^2> ) ; u' = u ; v' = v  
   ; --------------------------------------------------------
   ; ztke is 0.5 * <w'^2>/2 or 0.5 * sigma_w^2 
 ztke(*,yeye) = 0.5 * TOTAL(TOTAL(wprime^2,1),1) / float(nx) / float(ny)
 xtke(*,yeye) = 0.5 * TOTAL(TOTAL(getget(filesWRF(loop), 'U', anomaly=anomalu, count=[0,0,0,1], offset=[0,0,0,loop2])^2,1),1) / float(nx) / float(ny)
 ytke(*,yeye) = 0.5 * TOTAL(TOTAL(getget(filesWRF(loop), 'V', anomaly=anomalv, count=[0,0,0,1], offset=[0,0,0,loop2])^2,1),1) / float(nx) / float(ny)
; xtke(*,yeye) = 0.5 * TOTAL(TOTAL(uprime^2,1),1) / float(nx) / float(ny)
; ytke(*,yeye) = 0.5 * TOTAL(TOTAL(vprime^2,1),1) / float(nx) / float(ny)
  tke(*,yeye) = xtke(*,yeye) + $
                ytke(*,yeye) + $
                ztke(*,yeye)
   ; ------
   ; <w't'>
   ; ------
 wt(*,yeye)  = TOTAL(TOTAL(TEMPORARY(tprime)  * TEMPORARY(wprime),1),1) / float(nx) / float(ny)  
   ; ------
   ; p & ph
   ; ------
 ;if (loop + loop2 eq 0) then nloopbeware = nloop2-1 else nloopbeware = nloop2  ; 1ere valeur vaut 0
 nttotbeware = nttot-1 ; 1ere valeur vaut 0
 pht(*,yeye) = TOTAL(TOTAL(getget(filesWRF(loop), 'PHTOT',  count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny) / 1000. / 3.72
 pt(*,yeye)  = TOTAL(TOTAL(getget(filesWRF(loop), 'PTOT' ,  count=[0,0,0,1], offset=[0,0,0,loop2]),1),1) / float(nx) / float(ny)
 ph = TEMPORARY(ph) + pht(*,yeye) / nttotbeware ;/ nloopbeware / nloop1
 p  = TEMPORARY(p ) + pt(*,yeye) / nttot ;/ nloop2 / nloop1
   ; ----------------
   ; static stability
   ; ----------------
 stst(*,yeye) = DERIV( reform(pht(*,yeye)) - hgtu/1000. , reform(t(*,yeye)) * ( reform(pt(*,yeye)) /p0 )^r_cp ) + 1000.*grav / (R / r_cp)  ;; 4.9 dans Hinson
   ; ----------------
   ; surface pressure
   ; ----------------
 !QUIET=1
 ;psfc = getget(filesWRF(loop), 'PSFC', anomaly=1, count=[0,0,0,1], offset=[0,0,0,loop2])   ;; plus rapide que autre routine! pas grave avertissement
 psfc = getget(filesWRF(loop), 'PSFC', count=[0,0,0,1], offset=[0,0,0,loop2])
 !QUIET=0
; psfc = getget2d(filesWRF(loop), 'PSFC', anomaly=1, count=[0,0,1], offset=[0,0,loop2]) 
 psmin(yeye) = min(psfc) ;& print, psmin(yeye)
 w = where(TEMPORARY(psfc) lt -0.5) ;& print, n_elements(w)
 if (w(0) ne -1) then depressions(yeye) = 1000. * n_elements(w) / float(nx) / float(ny) else depressions(yeye) = 0.
   ; ----------
   ; loop count
   ; ----------
 yeye = TEMPORARY(yeye) + 1 ;& print, yeye & print, SYSTIME(1) - timetime, ' s'  

endfor
                                          print, 'file '+string(loop+1,'(I0)'), SYSTIME(1) - timetime, ' s'
endfor
h  = TEMPORARY(ph)  - hgtu/1000.  ;; altitude above ground
ht = TEMPORARY(pht) - hgtu/1000.  
;
; save
;
save, wt, tke, ztke, h, ht, t, p, pt, stst, localtime, xtke, ytke, wmax, wmin, depressions, psmin, filename='getturb.dat'
nz = n_elements(h)
nnn = n_elements(tke(0,*))

ENDIF ELSE BEGIN

print, 'OK, file is here'
restore, filename='getturb.dat'
nz = n_elements(h)
nnn = n_elements(tke(0,*))

OPENR, 23, 'input_more' & READF, 23, hgtu, tsurfu & CLOSE, 23

;@include_ustar.pro
;@include_hfx.pro
;@include_w.pro
;@include_tprime.pro

ENDELSE

;
; pbl height
;
pbl_height = fltarr(nnn) & hfsurf = fltarr(nnn) & w_star = fltarr(nnn) & hfbot = fltarr(nnn) & w_star_bot = fltarr(nnn) & hfbotmax = fltarr(nnn)
a_vel = fltarr(nz,nnn) & a_h = fltarr(nz,nnn) & a_wt = fltarr(nz,nnn) & a_vel_bot = fltarr(nz,nnn) & a_wt_bot = fltarr(nz,nnn)
tt_top = fltarr(nnn) & z_top = fltarr(nnn) & t_bot = fltarr(nnn) & p_bot = fltarr(nnn)
tt_bot = fltarr(nnn)
wmaxmax = fltarr(nnn)
wminmin = fltarr(nnn)
tkemax = fltarr(nnn)

for i=1, nnn-1 do begin    ;; ne pas traiter l'init 

	  ;
	  ; pbl height
	  ;
          z_in = reform(ht(*,i)) 
          profile_in = reform(stst(*,i))   ;; wt pas suffisamment regulier
          resol = 1000  ;200 un peu juste
         	  ;********************************************************************
	          ; ---- numerical recipes in C - section 3.3 ----
	          ;
	          ; Calculate interpolating cubic spline 
	                yspline = SPL_INIT(z_in,profile_in)
	          ; Define the X values P at which we desire interpolated Y values
	                xspline=min(z_in)+findgen(floor(max(z_in))*resol)/resol
	          ; Calculate the interpolated Y values
	                result=spl_interp(z_in,profile_in,yspline,xspline)
	          ;********************************************************************
	  w = where( (result ge 1.5) and (xspline gt 0.01) ) ;0.3) )   ;0.1 suffit si HR
	  ;if ( localtime(i) gt 16. ) then w = where( (result ge 1.5) and (xspline gt 2.) )

	if ( localtime(i) gt 15. ) then w = where( (result ge 1.5) and (xspline gt 1.) )  

	;;;; special tau = 5.
	;if ( localtime(i) ge 13. ) then w = where( (result ge 1.5) and (xspline gt 0.35) )
	;if ( localtime(i) ge 16.05 ) then w = [0]

	  pbl_height(i) = xspline(w(0)) ;& print, 'PBL height - mean profile', pbl_height(i)
	  z_top(i) = xspline(w(0)) + hgtu/1000.

          ;
          ; temperature @ pbl height
          ;
          z_in = reform(ht(*,i))
          profile_in = reform(t(*,i))   
	  ;profile_in = reform(t(*,i))*(reform(pt(*,i))/p0)^r_cp
	  profile_in2 = reform(pt(*,i))
          resol = 1000  ;200 un peu juste
                  ;********************************************************************
                  ; ---- numerical recipes in C - section 3.3 ----
                  ;
                  ; Calculate interpolating cubic spline 
                        yspline = SPL_INIT(z_in,profile_in)
			yspline2 = SPL_INIT(z_in,profile_in2)
                  ; Define the X values P at which we desire interpolated Y values
                        xspline=min(z_in)+findgen(floor(max(z_in))*resol)/resol
                  ; Calculate the interpolated Y values
                        result=spl_interp(z_in,profile_in,yspline,xspline)
			result2=spl_interp(z_in,profile_in2,yspline2,xspline)
                  ;********************************************************************
          w = where( abs(xspline-pbl_height(i)) eq min(abs(xspline-pbl_height(i))) ) & tt_top(i) = result[w(0)]*(result2[w(0)]/p0)^r_cp ;& print, z_top(i), t_top(i)
	  hhh = 1000. / 1000. ;100. / 1000.
	  w = where( abs(xspline-hhh) eq min(abs(xspline-hhh)) ) & t_bot(i) = result[w(0)] & p_bot(i) = result2[w(0)] ;& print, p_bot(i), t_bot(i)
	  tt_bot(i) = result[w(0)]*(result2[w(0)]/p0)^r_cp

	  ;
          ; heat flux at the "surface" 
	  ;
	  hfsurf(i) = wt(1,i)           ;; attention vaut 0 a la hauteur 0

          ;
          ; heat flux at the bottom of mixed layer = top of radiative layer
          ;
          z_in = reform(ht(*,i)) 
          profile_in = reform(wt(*,i))   ;; wt pas suffisamment regulier
          resol = 1000  ;200 un peu juste
                  ;********************************************************************
                  ; ---- numerical recipes in C - section 3.3 ----
                  ;
                  ; Calculate interpolating cubic spline 
                        yspline = SPL_INIT(z_in,profile_in)
                  ; Define the X values P at which we desire interpolated Y values
                        xspline=min(z_in)+findgen(floor(max(z_in))*resol)/resol
                  ; Calculate the interpolated Y values
                        result=spl_interp(z_in,profile_in,yspline,xspline)
                  ;********************************************************************
	  hhh = 300./1000. & yeah=where(abs(xspline-hhh) eq (min(abs(xspline-hhh)))) & nnzz=yeah(0) & hfbot(i)=result(nnzz)  ;; 300m semble OK, 500m super
	  w=where(xspline lt 1.5) & result=result[w] & xspline=xspline[w] & yeah=where(result eq max(result)) & nnzz=yeah(0) & hfbotmax(i)=result(nnzz)  
	  ;print, 'surf, 300m, max   ', hfsurf(i), hfbot(i), hfbotmax(i) ;& plot, result, xspline

result = reform(wmax(*,i)) & yeah=where(result eq max(result)) & nnzz=yeah(0) & wmaxmax(i)=result(nnzz)
result = reform(wmin(*,i)) & yeah=where(result eq min(result)) & nnzz=yeah(0) & wminmin(i)=result(nnzz)
result = reform(tke(*,i))  & yeah=where(result eq max(result)) & nnzz=yeah(0) & tkemax(i )=result(nnzz)

	  ;
	  ; free convection velocity scale
	  ;
;	  w_star_bot(i) = ( 1000.* grav * pbl_height(i) * hfbot(i) / t(0,i) ) ^ (1./3.)   
          w_star_bot(i) = ( 1000.* grav * pbl_height(i) * hfbotmax(i) / t(0,i) ) ^ (1./3.)
          w_star(i) 	= ( 1000.* grav * pbl_height(i) * hfsurf(i) / t(0,i) ) ^ (1./3.)     

	  ;
	  ; dimensional quantities
   	  ;
	  a_vel(*,i) 		= reform( 2.*ztke(*,i) / w_star(i)^2 )
          a_vel_bot(*,i) 	= reform( 2.*ztke(*,i) / w_star_bot(i)^2 )
	  a_h(*,i) 	= reform( ht(*,i) / pbl_height(i) )
	  a_wt(*,i) 	= reform( wt(*,i) / hfsurf(i) )		;reform( wt(*,i) / w_star(i) )
	  a_wt_bot(*,i) = reform( wt(*,i) / hfbotmax(i) )	;reform( wt(*,i) / w_star_bot(i) )

endfor
	  ;
	  ; free convection time scale
	  ;
	  fcts = 1000. * pbl_height / w_star / 3700.
	  ;
	  ; mixed layer temperature scale
	  ;
 	  mlts = hfsurf / w_star 

;pbl_height 	= SMOOTH(TEMPORARY(pbl_height), [smoothampl], /EDGE_TRUNCATE)
;hfsurf 	= SMOOTH(TEMPORARY(hfsurf), 	[smoothampl], /EDGE_TRUNCATE)
;w_star 	= SMOOTH(TEMPORARY(w_star), 	[smoothampl], /EDGE_TRUNCATE)
;hfbot 		= SMOOTH(TEMPORARY(hfbot), 	[smoothampl], /EDGE_TRUNCATE)
;w_star_bot 	= SMOOTH(TEMPORARY(w_star_bot), [smoothampl], /EDGE_TRUNCATE)
;hfbotmax 	= SMOOTH(TEMPORARY(hfbotmax), 	[smoothampl], /EDGE_TRUNCATE)
;a_vel 		= SMOOTH(TEMPORARY(a_vel), 	[0,smoothampl], /EDGE_TRUNCATE)
;a_h 		= SMOOTH(TEMPORARY(a_h), 	[0,smoothampl], /EDGE_TRUNCATE)
;a_wt 		= SMOOTH(TEMPORARY(a_wt), 	[0,smoothampl], /EDGE_TRUNCATE)
;a_vel_bot 	= SMOOTH(TEMPORARY(a_vel_bot), 	[0,smoothampl], /EDGE_TRUNCATE)
;a_wt_bot 	= SMOOTH(TEMPORARY(a_wt_bot), 	[0,smoothampl], /EDGE_TRUNCATE)
;tt_top		= SMOOTH(TEMPORARY(tt_top), 	[smoothampl], /EDGE_TRUNCATE)
;z_top          = SMOOTH(TEMPORARY(z_top),   	[smoothampl], /EDGE_TRUNCATE)
;t_bot          = SMOOTH(TEMPORARY(t_bot),   	[smoothampl], /EDGE_TRUNCATE)
;p_bot          = SMOOTH(TEMPORARY(p_bot),   	[smoothampl], /EDGE_TRUNCATE)

;
; save for comparison sake
;
save, pbl_height, hfbot, hfsurf, w_star, w_star_bot, a_vel_bot, a_vel, a_h, a_wt_bot, a_wt, fcts, mlts, hfbotmax, tt_top, z_top, t_bot, p_bot, tt_bot, wmaxmax, wminmin, tkemax, filename='compturb.dat'
if (noplot ne 'false') then stop

;
; smooth smooth
;
wt  = SMOOTH(TEMPORARY(wt),  [0,smoothampl], /EDGE_TRUNCATE)
tke = SMOOTH(TEMPORARY(tke), [0,smoothampl], /EDGE_TRUNCATE)
ztke = SMOOTH(TEMPORARY(ztke), [0,smoothampl], /EDGE_TRUNCATE)
xtke = SMOOTH(TEMPORARY(xtke), [0,smoothampl], /EDGE_TRUNCATE)
ytke = SMOOTH(TEMPORARY(ytke), [0,smoothampl], /EDGE_TRUNCATE)
wmax = SMOOTH(TEMPORARY(wmax), [0,smoothampl], /EDGE_TRUNCATE)
wmin = SMOOTH(TEMPORARY(wmin), [0,smoothampl], /EDGE_TRUNCATE)
depressions = SMOOTH(TEMPORARY(depressions), [smoothampl], /EDGE_TRUNCATE)
psmin = SMOOTH(TEMPORARY(psmin), [smoothampl], /EDGE_TRUNCATE)


;*******************;
;*******************;
; PLOTS PLOTS PLOTS ;
;*******************;
;*******************;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  a_vel_bot 
;;zey           =  a_h  
;;set_name      =  'A_VEL.ps'
;;set_title     =  "Scaled vertical velocity variance (m.s!U-1!N)"
;;set_titlex    =  'Scaled vertical velocity variance (m.s!U-1!N)'
;;set_titley    =  'Scaled height (km)'
;;set_subtitle  =  ''
;;set_xrange    =  [0.,1.2]
;;set_yrange    =  [0.,1.2]
;;set_tickx     =  0.1
;;set_ticky     =  0.1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zesim = 1.8 * a_h^(2./3.) * ( 1. - 0.8 * a_h )^2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;localtimes = localtime[where( (localtime ge 11) and (localtime le 16) )]
;;;localtimes = localtime
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;altlin=0
;;user_lt=localtimes[0] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;;plot, zefield(*,nntt), zey(*,nntt), xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin, /ISOTROPIC
;;for ll = 1, n_elements(localtimes)-1 do begin
;;  user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;;  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), linestyle=altlin
;;endfor
;;loadct, 33 & oplot, zesim, zey, psym=5, color=255 & loadct, 0 
;;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  localtime
zey           =  pbl_height
set_name      =  'HEIGHT.ps'
set_title     =  'Boundary layer height (km)'
set_titlex    =  'Local Time (h)'
set_titley    =  'Boundary layer height (km)'
set_subtitle  =  '' ;'Criterion is : static stability > 1.5 K.m!U-1!N'
set_xrange    =  [7.,19.] ;[11.,17.] ;[8.,17.]
set_yrange    =  [0.,9.]
set_tickx     =  1.
set_ticky     =  1.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;oplot, zefield, zey, psym=5
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  localtime
;;zey           =  w_star_bot
;;set_name      =  'W_STAR.ps'
;;set_title     =  "Free convection velocity scale (m.s!U-1!N)"
;;set_titlex    =  'Local Time (h)'
;;set_titley    =  'Free convection velocity scale (m.s!U-1!N)'
;;set_subtitle  =  ''
;;set_xrange    =  [8.,18.]
;;set_yrange    =  [0.,6.]
;;set_tickx     =  1.
;;set_ticky     =  0.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;;oplot, zefield, zey, psym=5
;;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  localtime
;;zey           =  hfbot
;;zey2	      =  hfsurf	
;;zey3          =  hfbotmax
;;set_name      =  'HFBOT.ps'
;;set_title     =  "Heat flux at the top of the radiative layer (K.m.s!U-1!N)"
;;set_titlex    =  'Local Time (h)'
;;set_titley    =  'Heat flux at the top of the radiative layer (K.m.s!U-1!N)'
;;set_subtitle  =  ''
;;set_xrange    =  [8.,18.]
;;set_yrange    =  [0.,2.5]
;;set_tickx     =  1.
;;set_ticky     =  0.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;oplot, zefield, zey2, linestyle=1
;;oplot, zefield, zey3, linestyle=2
;;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(tke) 
zex           =  localtime 
zey           =  h
set_name      =  'TKE.ps'
set_title     =  "Turbulent Kinetic Energy (m!U2!N.s!U-2!N)" ;0.5[<u'!U2!N>+<v'!U2!N>+<w'!U2!N>]
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  '' ;'Mean over the simulation domain'
set_xrange    =  [7.,19.] ;[8.,17.]
set_yrange    =  [0.,9.]  ;; [0.,200.]
set_tickx     =  1.
set_ticky     =  1. ;; 50.
minval        =  0.
maxval        =  20.
nlev          =  maxval-minval
pal           =  22
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0 
  ;; 2. color field
  loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b) 
              contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  2.*transpose(ztke)
;;set_title     =  "Vertical wind variance (m!U2!N.s!U-2!N)"
;;minval        =  -6.
;;maxval        =  12.
;;pal           =  0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(ztke)
zex           =  localtime
zey           =  h
set_name      =  'zTKE.ps'
set_title     =  "Vertical Turbulent Kinetic Energy (m!U2!N.s!U-2!N)" ;0.5[<w'!U2!N>]
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  '' ;'Mean over the simulation domain'
set_xrange    =  [7.,19.]
set_yrange    =  [0.,9.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  0.
maxval        =  8.
nlev          =  maxval-minval
pal           =  22
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;; 2. color field
loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;restore, filename='addturb.dat'
;;modvar = SMOOTH(TEMPORARY(modvar),  [0,smoothampl], /EDGE_TRUNCATE)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  transpose(modvar)
;;zex           =  localtime
;;zey           =  h
;;set_name      =  'modvar.ps'
;;set_title     =  "Horizontal wind speed variance (m!U2!N.s!U-2!N)"
;;set_titlex    =  'Local Time (h)'
;;set_titley    =  'Altitude above surface (km)'
;;set_subtitle  =  '' ;'Mean over the simulation domain'
;;set_xrange    =  [7.,17.]
;;set_yrange    =  [0.,7.]
;;set_tickx     =  1.
;;set_ticky     =  1.
;;minval        =  0.
;;maxval        =  10.
;;nlev          =  maxval-minval
;;pal           =  22
;;rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;;; 0. levels
;;lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;;;; 1. background
;;loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;;;; 2. color field
;;loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
;;            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;;;; 3. contour field
;;loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;;;; 4. choose output
;;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(wt)
zex           =  localtime
zey           =  h
set_name      =  'HF.ps'
set_title     =  "Turbulent Heat Flux (K.m.s!U-1!N)"  ;"Vertical Eddy Heat Flux <w'T'>" ;<w'!7h!3'>
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  '' ;'Mean over the simulation domain'
set_xrange    =  [7.,19.]
set_yrange    =  [0.,9.]  ;; [0.,200.]
set_tickx     =  1.
set_ticky     =  1. ;; 50.
minval        =  -1.5 ;-2.
maxval        =  2.5  ;2.
nlev          =  floor(maxval-minval)*10
pal           =  33  
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;minval        =  -0.8
;maxval        =  1.2
;pal           =  0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;; 2. color field
loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
            ;;;--------------------------------------------------------------------------------------------------------------------------------
            ;;; WHITE ZONE - 1. get location of interval in the CT - 2. change the CT to have a white zone
            ulim=0.09 & dlim=-0.09 & w=where(lev le dlim) & n1=w[n_elements(w)-1] & w=where(lev ge ulim) & n2=w[0] & yy=BYTSCL(lev) & nd=yy[n1] & nu=yy[n2]-5
            nu = nd + (nu-nd)/2  ;; otherwise the interval is too large (because we removed 0)
            TVLCT, r, g, b, /Get & r[nd:nu]=255 & g[nd:nu]=255 & b[nd:nu]=255 & TVLCT, r, g, b
            ;;;--------------------------------------------------------------------------------------------------------------------------------
            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;restore, filename='tpot_profB'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  t
;;zey           =  ht 
;;set_name      =  'T.ps'
;;set_title     =  "Potential Temperature (K)"
;;set_titlex    =  'Potential Temperature (K)'
;;set_titley    =  'Altitude above surface (km)'
;;set_subtitle  =  '' ;'Mean over the simulation domain'
;;set_xrange    =  [190.,230.] ;[min(t),max(t)]
;;set_yrange    =  [0.,7.]
;;set_tickx     =  5.
;;set_ticky     =  1.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;localtimes = [9,10,11,12,13,14,15,16,17,18]
;;localtimes = [7,9,11,13,15,17]
;;;localtimes = [17.0,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18.0]
;;;localtimes = [15.0]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;;altlin=0 & loadct, 0
;;altlin=1 & loadct, 0
;;user_lt=localtimes[0] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;;plot, zefield(*,nntt), zey(*,nntt), xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;;plot, zefield(*,nntt), les_column, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;for ll = 1, n_elements(localtimes)-1 do begin
;;  CASE altlin OF
;;  0: altlin=1
;;  1: altlin=0
;;  ENDCASE  
;;  user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;;  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), linestyle=altlin
;;;  if (nntt ne -1) then oplot, zefield(*,nntt), les_column, linestyle=altlin
;;endfor
;;;oplot, ro_tpot, ro_column, psym=7
;;;xyouts, 192.0, 0.25, '09:00'
;;;xyouts, 205.5, 0.25, '11:00'
;;;xyouts, 214.5, 0.25, '13:00'
;;;xyouts, 219.5, 0.25, '17:00'
;;;xyouts, 220.5, 2.30, '17:00 [RO]'
;;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  stst
;;zey           =  ht 
;;set_name      =  'STST.ps'
;;set_title     =  "Static stability (K.m!U-1!N)"
;;set_titlex    =  'Static stability (K.m!U-1!N)'
;;set_titley    =  'Altitude above surface (km)'
;;set_subtitle  =  'Mean over the simulation domain'
;;set_xrange    =  [-5.,5.]
;;set_yrange    =  [0.,7.]
;;set_tickx     =  1.
;;set_ticky     =  1.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;localtimes = [9,11,13,15,17]
;;localtimes = [12,13,14,15,16]
;;;localtimes = [17.0,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18.0]
;;;localtimes = [15.0]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;altlin=0 & loadct, 0
;;user_lt=localtimes[0] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;;plot, zefield(*,nntt), zey(*,nntt), xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;oplot, zefield(*,nntt), zey(*,nntt), psym=5
;;for ll = 1, n_elements(localtimes)-1 do begin
;;  CASE altlin OF
;;  0: altlin=1 
;;  1: altlin=0
;;  ENDCASE
;;  user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;;  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), linestyle=altlin
;;  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), psym=5
;;endfor
;;oplot, 0.*zefield(*,nntt) + 1.5, zey(*,nntt), linestyle=2  
;;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;if (n_elements(xtke) eq 0) then stop
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;zefield       =  2.*transpose(xtke)
;;;minval        =  -4
;;;maxval        =  8. ;12.
;;;pal           =  0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zex           =  localtime
;;zey           =  h
;;set_name      =  'xTKE.ps'
;;set_title     =  "Horizontal Turbulent Kinetic Energy (m!U2!N.s!U-2!N)" ;0.5[<u'!U2!N>]
;;set_titlex    =  'Local Time (h)'
;;set_titley    =  'Altitude above surface (km)'
;;set_subtitle  =  '';'Mean over the simulation domain'
;;set_xrange    =  [8.,17.]
;;set_yrange    =  [0.,7.] 
;;set_tickx     =  1.
;;set_ticky     =  1.
;;minval        =  0.
;;maxval        =  8. 
;;nlev          =  maxval-minval
;;pal           =  22
;;rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;;; 0. levels
;;lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;;;; 1. background
;;loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;;;; 2. color field
;;loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
;;            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;;;; 3. contour field
;;loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;;;; 4. choose output
;;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zefield       =  transpose(ytke) + transpose(xtke)
;zex           =  localtime
;zey           =  h
;set_name      =  'hTKE.ps'
;set_title     =  "Horizontal Turbulent Kinetic Energy (m!U2!N.s!U-2!N)" ;0.5[<v'!U2!N>]
;set_titlex    =  'Local Time (h)'
;set_titley    =  'Altitude above surface (km)'
;set_subtitle  =  '' ;'Mean over the simulation domain'
;set_xrange    =  [8.,17.]
;set_yrange    =  [0.,1.] ;[0.,7.] ;; [0.,200.]
;set_tickx     =  1.
;set_ticky     =  0.1 ;1. ;50.
;minval        =  0.
;maxval        =  12. ;10.
;nlev          =  maxval-minval
;pal           =  22
;rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;!P.Charsize = 1.2
;;; 0. levels
;lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;;; 1. background
;loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;;; 2. color field
;loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
;            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;;; 3. contour field
;loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;;; 4. choose output
;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (n_elements(wmax) eq 0) then stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(wmax)
zex           =  localtime
zey           =  h
set_name      =  'Wmax.ps'
set_title     =  "Maximum updraft speed (m.s!U-1!N)"
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  ''
set_xrange    =  [11.,17.] ;[8.,18.]
set_yrange    =  [0.,8.] ;[0.,8.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  0.
maxval        =  18. ;16. ;15.;14.;13.
nlev          =  maxval-minval
pal           =  22
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
levhalf = minval + (maxval-minval)*findgen((nlev/2)+1)/float(nlev/2) & if (minval ne 0.) then levhalf = levhalf[where(levhalf ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
  ;; 2. color field
  loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
              contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=levhalf, c_labels=findgen(n_elements(levhalf))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (levhalf LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  abs(transpose(wmin))
zex           =  localtime
zey           =  h
set_name      =  'Wmin.ps'
set_title     =  "Maximum downdraft speed (m.s!U-1!N)"
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  ''
set_xrange    =  [11.,17.] ;[8.,18.]
set_yrange    =  [0.,8.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  0.
maxval        =  10.;12.
nlev          =  maxval-minval
pal           =  22
rrr           =  'no' ;'yes'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
levhalf = minval + (maxval-minval)*findgen((nlev/2)+1)/float(nlev/2) & if (minval ne 0.) then levhalf = levhalf[where(levhalf ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;; 2. color field
loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=levhalf, c_labels=findgen(n_elements(levhalf))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (levhalf LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


restore, filename='addturb.dat'
velmax = SMOOTH(TEMPORARY(velmax),  [0,smoothampl], /EDGE_TRUNCATE)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(velmax)
zex           =  localtime
zey           =  h
set_name      =  'velmax.ps'
set_title     =  "Maximum horizontal wind speed (m.s!U-1!N)"
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  '' ;'Mean over the simulation domain'
set_xrange    =  [07.,19.] ;[8.,17.]
set_yrange    =  [0.,4.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  0.
maxval        =  30. ;12.
nlev          =  maxval-minval
pal           =  22
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels 
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
;;; 2. color field
loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
            contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  localtime
;;zey           =  depressions 
;;set_name      =  'DEPR.ps'
;;set_title     =  "Percentage of depressions in the domain"
;;set_titlex    =  'Local Time (h)'
;;set_titley    =  'Percentage'
;;set_subtitle  =  ''
;;set_xrange    =  [8.,18.]
;;set_yrange    =  [0.,3.]
;;set_tickx     =  1.
;;set_ticky     =  1.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;oplot, zefield, zey, psym=5
;;	oplot, localtime, abs(psmin), linestyle=1
;;if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
goto, no_staticstab
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then begin
  device, /close
  set_plot, 'ps' & device, filename='plot/staticstab.ps'
endif

user_lt=17. & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))

   ;;; recompute height @ given local time 
   caca=heightp+hgtu/1000.
   height=reform(ph(*,0,*,*))
   height=total(height,1)/n_elements(height(*,0,0))
   height=reform(height)
   height=reform(height(*,yeah(0)))
   height=height/1000./3.72
   heightp=height(0:n_elements(height(*))-2) 
   print, 'new minus old', heightp - caca
   ;;; recompute height @ given local time 

        press=reform(p(*,0,*,*))
        press=total(press,1)/n_elements(press(*,0,0))
        press=reform(press)
        press=reform(press(*,yeah(0)))
        press=press(0:n_elements(press(*))-2)

staticstab = DERIV( heightp , reform(what_I_plot5(*,yeah(0))) ) + 1000.*3.72/844.6  ;; 4.9 dans Hinson
   staticstab = staticstab(0:n_elements(staticstab)-5)
   heightp = heightp(0:n_elements(heightp)-5)

plot, $
        staticstab, $
        heightp,$
        xtitle='Static stability (K/km)',$
        xrange=[-1.,5.], $
        xtickinterval=0.5, $
        ytitle='Altitude above surface (km)', $
        yrange=[min(heightp),max(heightp)], $
        ytickinterval=1., $
        title="LMD LES Static stability @ LT 17h (K/km)", $
        subtitle='zonal average at lat. '+latwrite
oplot, $
        staticstab, $
        heightp,$
        psym=5
oplot, $
       findgen(n_elements(heightp))*0. + 1.,$
       heightp,$
       linestyle=2
oplot, $
       findgen(n_elements(heightp))*0. + 2.,$
       heightp,$
       linestyle=2

;;;;;;;;;;
w = where((staticstab gt 1.5) and ((heightp- hgtu/1000.) gt 1.))
t_top_plus = what_I_plot5(w(0),yeah(0))
z_top_plus = heightp(w(0))
p_top_plus = press(w(0))
t_top_moins = what_I_plot5(w(0)-1,yeah(0))
z_top_moins = heightp(w(0)-1)
p_top_moins = press(w(0)-1)
pbl_depth = (z_top_plus*alog(p_top_plus) + z_top_moins*alog(p_top_moins))/(alog(p_top_plus) + alog(p_top_moins)) - hgtu/1000.
xyouts, 3., 1.5 + (max(heightp) + min(heightp)) / 3., 'Ls = '+string(lsu,'(F5.1)')+'!Uo!N', CHARSIZE=1
xyouts, 3., 1. + (max(heightp) + min(heightp)) / 3., 'Lat = '+string(latu,'(F5.1)')+'!Uo!N', CHARSIZE=1
xyouts, 3., 0.5 + (max(heightp) + min(heightp)) / 3., 'LonE = '+string(lonu,'(F6.1)')+'!Uo!N', CHARSIZE=1
xyouts, 3., (max(heightp) + min(heightp)) / 3., 'T!Dt!N = '+string(t_top_plus,'(I0)')+'/'+string(t_top_moins,'(I0)')+' K ', CHARSIZE=1
xyouts, 3., -0.5 + (max(heightp) + min(heightp)) / 3., 'p!Dt!N = '+string(p_top_plus,'(I0)')+'/'+string(p_top_moins,'(I0)')+' Pa', CHARSIZE=1
xyouts, 3., -1. + (max(heightp) + min(heightp)) / 3., 'z!Dt!N = '+string(z_top_plus,'(F4.1)')+'/'+string(z_top_moins,'(F4.1)')+' km', CHARSIZE=1
xyouts, 3., -1.5 + (max(heightp) + min(heightp)) / 3., 'z!Ds!N = '+string(hgtu/1000.,'(F4.1)')+' km', CHARSIZE=1
xyouts, 3., -2. + (max(heightp) + min(heightp)) / 3., 'D = '+string(pbl_depth,'(F4.1)')+' km', CHARSIZE=1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
no_staticstab:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then begin
  device, /close
endif
stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


user_lt=13.
yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt))))

user_h=0.1
walt=where(abs(heightp-user_h) eq (min(abs(heightp-user_h))))

;mapfield=reform(t[*,*,walt(0),w(0)])
;help, mapfield
;contour, mapfield

section=reform(w[*,0,*,yeah(0)])

section=reform(w[*,0,*,160])

lev=[-12.,-8.,-4.,4.,8.,12.]
lev=[-5.,-4.,-3.,-2.,-1.,1.,2.,3.,4.,5.]
contour, $
	section, $
	(lon(*,0)-lon(0,0))*59., $
	heightp, $
	levels=lev , $
	C_LINESTYLE = (lev LT 0.0)

device, /close

end

