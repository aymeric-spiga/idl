pro compturb
saveps='true' & p0=610. & t0=220. & r_cp=1/4.4 & grav=3.72 & R=192.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; grep frames_per_outfile */*/namelist.input
history_interval_s = 100.
smoothampl=3700/history_interval_s
;smoothampl=10.
;smoothampl=20.
;smoothampl=4.
smoothampl=0.
smoothampl=5.
;;;;;;;
;
; ATTENTION IL FAUT CHANGER LES VARIABLES A DEUX ENDROITS
;
;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;
;;; SET HERE WHAT YOU WANT TO COMPARE
;;;
;@article_qjrms.inc
@report_tasi.inc
;;;
;;;
;;;


;;;
;;; TWEAKS
;;;
;SPAWN, 'cp compturb.pro '+saveplot
;SPAWN, 'gqview '+saveplot+' &'
;goto, quick


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; pour article, decommenter ici et commenter l'appel au deuxieme fichier
;;;com
;;!p.multi = [0,2,1]
;
;        path = paths(0) & print, path 
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;
;;; COEFFICIENT POUR LOIS DE SIMILITUDE
;;;
;;@coeffhf
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zefield       =  a_wt_bot 
;zey           =  a_h
;set_name      =  saveplot+'comp_A_HFLUX.ps'
;set_title     =  "Similarity: vert. eddy heat flux on Mars"
;set_titlex    =  "Dimensionless vertical heat flux <w'T'>/<w'T'>!Dmax!N"
;set_titley    =  'Dimensionless height z/z!Di!N'
;set_subtitle  =  ''
;set_xrange    =  [-0.5,1.]
;set_yrange    =  [0.,1.5]
;set_tickx     =  0.25
;set_ticky     =  0.1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zesim = 1 - 1.2 * a_h
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;localtimes = localtime[where( (localtime ge 11) and (localtime le 16) )]
;;;localtimes = localtime[where( (localtime ge 13) and (localtime le 16) )]
;;;localtimes = localtime[where( (localtime ge 12) and (localtime le 15) )]
;;localtimes = localtime[where( (localtime ge 13) and (localtime le 15) )]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;!P.Charsize = 1.2
;user_lt=localtimes[0] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;plot, zefield(*,nntt), zey(*,nntt), xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, psym=3;, /ISOTROPIC;, POS=[0.12, 0.22, 0.42, 1.0] ;color=0, linestyle=altlin
;for ll = 1, n_elements(localtimes)-1 do begin
;  user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), psym=3 ;linestyle=altlin
;endfor
;;wherex = set_xrange(0)
;;wherey = set_yrange(1)
;;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + dddx, wherey, path
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zecolor = 0. & loadct, 33
;for i=1, n_elements(paths)-1 do begin
;;        zecolor = 50. + zecolor
;zecolor = 0. & loadct, 0
;        path = paths(i) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        zefield = a_wt_bot
;        zey = a_h
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        for ll = 0, n_elements(localtimes)-1 do begin
;          user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;          if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), psym=3, color=zecolor
;        endfor
;        ;xyouts, wherex + dddx, wherey - float(i)*dddy, path, color=zecolor
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;endfor
;loadct, 0
;;oplot, [0.,0.], set_yrange, linestyle=0, color=0
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	;ye, ryeah, tyeah & oplot, ryeah, tyeah, linestyle=0, color=255, /NOCLIP
;	;oplot, set_xrange, [0.,0.], linestyle=0, color=0
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_End, /PNG
;;com
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;	path = paths(0) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;
;;; COEFFICIENT POUR LOIS DE SIMILITUDE
;;;
;;@coeffw
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zefield       =  a_vel_bot 
;zey           =  a_h  
;set_name      =  saveplot+'comp_A_VEL.ps'
;set_title     =  "Similarity: vert. vel. variance on Mars"
;set_titlex    =  "Dimensionless vert. vel. variance <w'!U2!N>/w!D*!N!U2!N"
;set_titley    =  'Dimensionless height z/z!Di!N'
;set_subtitle  =  ''
;set_xrange    =  [0.,1.2]
;set_yrange    =  [0.,1.5]
;set_tickx     =  0.25 ;0.1
;set_ticky     =  0.1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zesim = 1.8 * a_h^(2./3.) * ( 1. - 0.8 * a_h )^2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;localtimes = localtime[where( (localtime ge 11) and (localtime le 16) )]
;;localtimes = localtime[where( (localtime ge 13) and (localtime le 15) )]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;com
;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;!P.Charsize = 1.2
;user_lt=localtimes[0] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;plot, zefield(*,nntt), zey(*,nntt), xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, psym=3;, /ISOTROPIC;, POS=[0.62, 0.22, 0.92, 1.0] ;color=0, linestyle=altlin
;for ll = 1, n_elements(localtimes)-1 do begin
;  user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), psym=3 ;linestyle=altlin
;endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zecolor = 0. & loadct, 33
;for i=1, n_elements(paths)-1 do begin
;zecolor = 0. & loadct, 0
;;	zecolor = 50. + zecolor
;	path = paths(i) & print, path
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	file = path+'/compturb.dat' & restore, filename=file
;	file = path+'/getturb.dat'  & restore, filename=file
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	zefield = a_vel_bot 
;	zey = a_h
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	for ll = 0, n_elements(localtimes)-1 do begin
;	  user_lt=localtimes[ll] & yeah=where(abs(localtime-user_lt) eq (min(abs(localtime-user_lt)))) & nntt = yeah(0)
;	  if (nntt ne -1) then oplot, zefield(*,nntt), zey(*,nntt), psym=3, color=zecolor
;	endfor
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;endfor
;loadct, 0
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	;ye2, ryeah, tyeah & oplot, ryeah, tyeah, linestyle=0, color=255, /NOCLIP
;	;oplot, set_xrange, [0.,0.], linestyle=0, color=0
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_End, /PNG
;
;;;com
;;stop
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;        path = paths(0) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zefield       =  localtime
;zey           =  SMOOTH(w_star_bot, [smoothampl], /EDGE_TRUNCATE)
;set_name      =  saveplot+'comp_W_STAR.ps'
;set_title     =  '';"Free convection velocity scale W!D*!N (m.s!U-1!N)"
;set_titlex    =  'Local Time (h)'
;set_titley    =  'Free convection velocity scale (m s!U-1!N)'
;set_subtitle  =  ''
;set_xrange    =  [08.,17.]
;set_yrange    =  [00.,06.]
;set_tickx     =  1.0
;set_ticky     =  0.5
;wherex	      =  8.5;11.0
;wherey	      =  5.5;1.75
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;fixfix = where((w_star_bot le 0.) and (localtime ge 15.)) & if (fixfix(0) ne -1) then zey(fixfix(0):n_elements(zey)-1) = !Values.F_NAN
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;!P.Charsize = 1.2
;altlin=0
;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + 2.*dddx, wherey, tttitle(0);path
;        plots,  wherex + 2.*dddx/1.1,  wherey
;        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zecolor = 0. & loadct, 33
;for i=1, n_elements(paths)-1 do begin
;        ;zecolor = float(255)/n_elements(paths) + zecolor
;	zecolor = 0
;        path = paths(i) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        zefield = localtime
;        zey = SMOOTH(w_star_bot, [smoothampl], /EDGE_TRUNCATE);w_star_bot
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	fixfix = where((w_star_bot le 0.) and (localtime ge 15.)) & if (fixfix(0) ne -1) then zey(fixfix(0):n_elements(zey)-1) = !Values.F_NAN
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
;        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, tttitle(i), color=zecolor
;        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
;        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;endfor
;loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_End, /PNG
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        path = paths(0) & print, path
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        file = path+'/compturb.dat' & restore, filename=file
        file = path+'/getturb.dat'  & restore, filename=file
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  localtime
zey           =  SMOOTH(pbl_height, [smoothampl], /EDGE_TRUNCATE)
set_name      =  saveplot+'comp_HEIGHT.ps'
set_title     =  '';"Boundary layer height z!Di!N (km)"; [definition : static stability > 1.5 K.m!U-1!N] "
set_titlex    =  'Local Time (h)'
set_titley    =  'Boundary layer height (km)'
set_subtitle  =  ''
set_xrange    =  [08.,18.]
set_yrange    =  [00.,08.] ;[00.,07.]
set_tickx     =  1. 
set_ticky     =  1.
wherex        =  09.;13.25 ;12.75
wherey        =  06.5;2.5 ;2.0
;	set_xrange    =  [13.5,18.] & set_yrange    =  [06.4,08.0] & set_tickx     =  0.5 & set_ticky     =  0.2 
;       wherex=15.30
;       wherey=06.90
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
altlin=0
plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + 2.*dddx, wherey, tttitle(0)
        plots,  wherex + 2.*dddx/1.1,  wherey
        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;w = where((localtime ge 17.) and (localtime le 17.5)) & if (w(0) ne -1) then print, 'BL height 17:00-17:30 : ', mean(pbl_height[w])
w = where((localtime ge 14.) and (localtime le 14.5)) & if (w(0) ne -1) then print, 'BL height 14:00-14:30 : ', mean(pbl_height[w])
;print, max(SMOOTH(tke, [0,37], /EDGE_TRUNCATE))
;print, max(SMOOTH(wmax, [0,37], /EDGE_TRUNCATE))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zecolor = 0. & loadct, 33
for i=1, n_elements(paths)-1 do begin
        ;zecolor = float(255)/n_elements(paths) + zecolor
	zecolor=0
        path = paths(i) & print, path
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        file = path+'/compturb.dat' & restore, filename=file
        file = path+'/getturb.dat'  & restore, filename=file
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        zefield = localtime
        zey = SMOOTH(pbl_height, [smoothampl], /EDGE_TRUNCATE);pbl_height 
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, tttitle(i), color=zecolor
        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;w = where((localtime ge 17.) and (localtime le 17.5)) & if (w(0) ne -1) then print, 'BL height 17:00-17:30 : ', mean(pbl_height[w])
        w = where((localtime ge 14.) and (localtime le 14.5)) & if (w(0) ne -1) then print, 'BL height 14:00-14:30 : ', mean(pbl_height[w])
	;print, max(SMOOTH(tke, [0,37], /EDGE_TRUNCATE))
	;print, max(SMOOTH(wmax, [0,37], /EDGE_TRUNCATE))
endfor
loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_End, /PNG

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;        path = paths(0) & print, path
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        file = path+'/compturb.dat' & restore, filename=file
;;        file = path+'/getturb.dat'  & restore, filename=file
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  localtime
;;zey           =  smooth(DERIV(localtime, smooth(pbl_height,smoothampl,/EDGE_TRUNCATE)),smoothampl,/EDGE_TRUNCATE)  
;;set_name      =  saveplot+'comp_dHEIGHT.ps'
;;set_title     =  "Boundary layer height growth dz!Di!N/dt (km/hour)"
;;set_titlex    =  'Local Time (h)'
;;set_titley    =  'Boundary layer height growth (km/hour)'
;;set_subtitle  =  ''
;;set_xrange    =  [9.,18.]
;;set_yrange    =  [-0.4,1.6] 
;;set_tickx     =  1.
;;set_ticky     =  0.2
;;wherex        =  10.30
;;wherey        =  0.1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;altlin = 0
;;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + 2.*dddx, wherey, path
;;        plots,  wherex + 2.*dddx/1.1,  wherey
;;        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;zecolor = 0. & loadct, 33
;;for i=1, n_elements(paths)-1 do begin
;;        ;zecolor = float(255)/n_elements(paths) + zecolor
;;	zecolor = 0
;;        path = paths(i) & print, path
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        file = path+'/compturb.dat' & restore, filename=file
;;        file = path+'/getturb.dat'  & restore, filename=file
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        zefield = localtime
;;        zey = smooth(DERIV(localtime, smooth(pbl_height,smoothampl,/EDGE_TRUNCATE)),smoothampl,/EDGE_TRUNCATE) 
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
;;        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, path, color=zecolor
;;        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
;;        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;endfor
;;loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_End, /PNG
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;        path = paths(0) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zefield       =  localtime
;zey           =  SMOOTH(wmaxmax, [smoothampl], /EDGE_TRUNCATE)
;zey2          =  SMOOTH(wminmin, [smoothampl], /EDGE_TRUNCATE)
;set_name      =  saveplot+'comp_WMAX.ps'
;set_title     =  '';"Maximum and minimum vertical velocity (m.s!U-1!N)"
;set_titlex    =  'Local Time (h)'
;set_titley    =  'Maximum and minimum vertical velocity (m s!U-1!N)'
;set_subtitle  =  ''
;;set_xrange    =  [09.0,18.0]
;set_yrange    =  [-12.0,18.0];[-15.0,20.0]
;set_tickx     =  1.
;set_ticky     =  2.
;wherex        =  8.5;11.00
;wherey        =  16.5;5.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;!P.Charsize = 1.2
;altlin = 0
;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;oplot, zefield, zey2, color=0, linestyle=altlin
;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + 2.*dddx, wherey, tttitle(0)
;        plots,  wherex + 2.*dddx/1.1,  wherey
;        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zecolor = 0. & loadct, 33
;for i=1, n_elements(paths)-1 do begin
;        ;zecolor = float(255)/n_elements(paths) + zecolor
;        zecolor = 0.
;        path = paths(i) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        zefield = localtime
;        zey = SMOOTH(wmaxmax, [smoothampl], /EDGE_TRUNCATE)
;	zey2 = SMOOTH(wminmin, [smoothampl], /EDGE_TRUNCATE) 
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
;        oplot, zefield, zey2, color=zecolor, linestyle=altlin+i
;        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, tttitle(i), color=zecolor
;        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
;        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;endfor
;loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_End, /PNG
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;        path = paths(0) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zefield       =  localtime
;zey           =  SMOOTH(tkemax, [smoothampl], /EDGE_TRUNCATE);hfbotmax
;set_name      =  saveplot+'comp_TKE.ps'
;set_title     =  '';"Maximum turbulent kinetic energy (m!U2!N.s!U-2!N)" ; 0.5[<u'!U2!N>+<v'!U2!N>+<w'!U2!N>]!Dmax!N (m!U2!N.s!U-2!N)"
;set_titlex    =  'Local Time (h)'
;set_titley    =  'Maximum turbulent kinetic energy (m!U2!N s!U-2!N)'
;set_subtitle  =  ''
;;set_xrange    =  [09.0,18.0]
;set_yrange    =  [00.0,18.0];[00.0,26.0]
;set_tickx     =  1.
;set_ticky     =  2.
;wherex        =  8.5;11.70
;wherey        =  16.;5.4
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;!P.Charsize = 1.2
;altlin = 0
;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/25. & xyouts, wherex + 2.*dddx, wherey, tttitles(0)
;        plots,  wherex + 2.*dddx/1.1,  wherey
;        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zecolor = 0. & loadct, 33
;for i=1, n_elements(paths)-1 do begin
;        ;zecolor = float(255)/n_elements(paths) + zecolor
;        zecolor = 0.
;        path = paths(i) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        zefield = localtime
;        zey = SMOOTH(tkemax, [smoothampl], /EDGE_TRUNCATE);hfbotmax 
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
;        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, tttitles(i), color=zecolor
;        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
;        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;endfor
;loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_End, /PNG
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;        path = paths(0) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zefield       =  localtime
;zey           =  SMOOTH(hfbotmax, [smoothampl], /EDGE_TRUNCATE);hfbotmax
;set_name      =  saveplot+'comp_HFBOT.ps'
;set_title     =  '';"Maximum heat flux <w'T'>!Dmax!N (K.m.s!U-1!N)"
;set_titlex    =  'Local Time (h)'
;set_titley    =  'Maximum vertical eddy heat flux (K m s!U-1!N)'
;set_subtitle  =  ''
;;set_xrange    =  [09.0,18.0]
;set_yrange    =  [00.0,02.0];[00.0,02.8]
;set_tickx     =  1.
;set_ticky     =  0.2
;wherex        =  8.5;10.30
;wherey        =  1.8;0.7
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_Start, FILENAME=set_name 
;!P.Charsize = 1.2
;altlin = 0
;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + 2.*dddx, wherey, tttitle(0)
;        plots,  wherex + 2.*dddx/1.1,  wherey
;        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zecolor = 0. & loadct, 33
;for i=1, n_elements(paths)-1 do begin
;        ;zecolor = float(255)/n_elements(paths) + zecolor
;	zecolor = 0.
;        path = paths(i) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/getturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        zefield = localtime
;        zey = SMOOTH(hfbotmax, [smoothampl], /EDGE_TRUNCATE);hfbotmax 
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
;        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, tttitle(i), color=zecolor
;        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
;        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;endfor
;loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_End, /PNG 
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;        path = paths(0) & print, path
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        file = path+'/compturb.dat' & restore, filename=file
;;        file = path+'/getturb.dat'  & restore, filename=file
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  localtime
;;zey           =  t_bot;SMOOTH(t_bot, [smoothampl], /EDGE_TRUNCATE);tt_top
;;zey	      =  SMOOTH(DERIV(localtime, SMOOTH(t_bot,smoothampl/2.,/EDGE_TRUNCATE)),smoothampl,/EDGE_TRUNCATE)
;;;zey           =  SMOOTH(DERIV(localtime, SMOOTH(t_bot/((p0 / p_bot)^r_cp),smoothampl/2.,/EDGE_TRUNCATE)),smoothampl,/EDGE_TRUNCATE)
;;set_name      =  saveplot+'comp_TdtBOT.ps';saveplot+'comp_TBOT.ps'
;;set_title     =  "Potential temperature time derivative dT!Db!N/dt (K.hour!U-1!N)"
;;set_titlex    =  'Local Time (h)'
;;set_titley    =  'Potential temperature time derivative (K.hour!U-1!N)'
;;set_subtitle  =  ''
;;set_xrange    =  [ 09., 18.] 
;;set_yrange    =  [-01., 06.] 
;;set_tickx     =  1.
;;set_ticky     =  1. 
;;wherex        =  13.25 
;;wherey        =  5.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;altlin = 0
;;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + 2.*dddx, wherey, tttitle(0)
;;        plots,  wherex + 2.*dddx/1.1,  wherey
;;        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;zecolor = 0. & loadct, 33
;;for i=1, n_elements(paths)-1 do begin
;;        ;zecolor = float(255)/n_elements(paths) + zecolor
;;        zecolor = 0.
;;        path = paths(i) & print, path
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        file = path+'/compturb.dat' & restore, filename=file
;;        file = path+'/getturb.dat'  & restore, filename=file
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        zefield = localtime
;;        zey = t_bot;SMOOTH(t_bot, [smoothampl], /EDGE_TRUNCATE);tt_top 
;;		;zey = DERIV(localtime, t_bot)
;;		zey = SMOOTH(DERIV(localtime, SMOOTH(t_bot,smoothampl/2.,/EDGE_TRUNCATE)),smoothampl,/EDGE_TRUNCATE)
;;		;zey = SMOOTH(DERIV(localtime, SMOOTH(t_bot/((p0 / p_bot)^r_cp),smoothampl/2.,/EDGE_TRUNCATE)),smoothampl,/EDGE_TRUNCATE)
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
;;        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, tttitle(i), color=zecolor
;;        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
;;        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;endfor
;;loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;oplot, set_xrange, [0,0], linestyle=2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_End, /PNG
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;        path = paths(0) & print, path
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        file = path+'/compturb.dat' & restore, filename=file
;;        file = path+'/getturb.dat'  & restore, filename=file
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zefield       =  localtime
;;zey           =  SMOOTH(tt_top, [smoothampl], /EDGE_TRUNCATE);tt_top
;;set_name      =  saveplot+'comp_TTOP.ps'
;;set_title     =  "Temperature at the top of the boundary layer (K)"
;;set_titlex    =  'Local Time (h)'
;;set_titley    =  'Temperature at the top of the boundary layer (K)'
;;set_subtitle  =  ''
;;set_xrange    =  [09.,18.]
;;set_yrange    =  [195.,225.]
;;set_tickx     =  1.
;;set_ticky     =  5.
;;wherex        =  13.00
;;wherey        =  223.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_Start, FILENAME=set_name
;;!P.Charsize = 1.2
;;altlin = 0
;;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + 2.*dddx, wherey, path
;;        plots,  wherex + 2.*dddx/1.1,  wherey
;;        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;zecolor = 0. & loadct, 33
;;for i=1, n_elements(paths)-1 do begin
;;        ;zecolor = float(255)/n_elements(paths) + zecolor
;;	zecolor = 0.
;;        path = paths(i) & print, path
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        file = path+'/compturb.dat' & restore, filename=file
;;        file = path+'/getturb.dat'  & restore, filename=file
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        zefield = localtime
;;        zey = SMOOTH(tt_top, [smoothampl], /EDGE_TRUNCATE);tt_top 
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
;;        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, path, color=zecolor
;;        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
;;        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
;;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;endfor
;;loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;if (saveps eq 'true') then PS_End, /PNG
;
;
;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;
;lev = 4
;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;
;
;        path = paths(0) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/addturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zefield       =  localtime
;zey 	      =  SMOOTH(reform(velmax(lev,*)), [smoothampl], /EDGE_TRUNCATE)
;set_name      =  saveplot+'comp_VELMAX.ps'
;set_title     =  '';"Maximum horizontal velocity "+string(h(lev)*1000.,'(I0)')+"m above ground (m.s!U-1!N)"
;set_titlex    =  'Local Time (h)'
;set_titley    =  "Maximum horizontal velocity "+string(h(lev)*1000.,'(I0)')+"m above ground (m s!U-1!N)"
;set_subtitle  =  ''
;;set_xrange    =  [09.0,18.0]
;set_yrange    =  [00.0,14.0]
;set_tickx     =  1.
;set_ticky     =  2.
;wherex        =  8.5;10.30
;wherey        =  12.;0.7
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;print, h(lev)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_Start, FILENAME=set_name 
;!P.Charsize = 1.2
;altlin = 0
;plot, zefield, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
;dddx = (set_xrange(1) - set_xrange(0))/20. & dddy = (set_yrange(1) - set_yrange(0))/20. & xyouts, wherex + 2.*dddx, wherey, tttitle(0)
;        plots,  wherex + 2.*dddx/1.1,  wherey
;        plots,  wherex + 2.*dddx/20.,  wherey, linestyle=altlin, /continue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;zecolor = 0. & loadct, 33
;for i=1, n_elements(paths)-1 do begin
;        ;zecolor = float(255)/n_elements(paths) + zecolor
;        zecolor = 0.
;        path = paths(i) & print, path
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        file = path+'/compturb.dat' & restore, filename=file
;        file = path+'/addturb.dat'  & restore, filename=file
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        zefield = localtime
;	zey = SMOOTH(reform(velmax(lev,*)), [smoothampl], /EDGE_TRUNCATE)
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        oplot, zefield, zey, color=zecolor, linestyle=altlin+i
;        xyouts, wherex + 2.*dddx,      wherey - float(i)*dddy, tttitle(i), color=zecolor
;        plots,  wherex + 2.*dddx/1.1,  wherey - float(i)*dddy
;        plots,  wherex + 2.*dddx/20.,  wherey - float(i)*dddy, linestyle=altlin+i, /continue
;        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;endfor
;loadct, 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if (saveps eq 'true') then PS_End, /PNG






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


goto, no_therm
quick:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
cp = R / r_cp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x  = fltarr(n_elements(paths)) & y  = fltarr(n_elements(paths))
x2 = fltarr(n_elements(paths)) & y2 = fltarr(n_elements(paths))
for i=0, n_elements(paths)-1 do begin
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	path = paths(i) & print, path
        file = path+'/compturb.dat' & restore, filename=file
        file = path+'/getturb.dat'  & restore, filename=file
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	yeye = 17. 
	;yeye = 15.
        ;yeye = 13.
        ;yeye = 09.
	w = where(localtime eq yeye) ;& print, w


        	;yeye = 08. 
        	;w2 = where(localtime eq yeye) ;& print, w
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	x(i) = (p0 / p_bot[w])^r_cp
	y(i) = t_bot[w];-t_bot[w2] & print, y(i)
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	y2(i) = z_top[w]
	;y2(i) = pbl_height[w]
	x2(i) = t_bot[w]
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;w2 = where( abs(h-0.01) eq min(abs(h-0.01)) )
;w2(0)=0
;w2 = where( abs(h-0.1) eq min(abs(h-0.1)) )
;w2 = where( abs(h-1.5) eq min(abs(h-1.5)) )
;w2 = where( abs(h-1.) eq min(abs(h-1.)) )
;print, w2(0), w(0)
;x(i) = (p0 / pt[w2(0),w(0)])^r_cp
;y(i) = t[w2(0),w(0)]
;print, h[w2]
;stop

;ye=smooth(DERIV(localtime, smooth(t_bot,smoothampl,/EDGE_TRUNCATE)),smoothampl,/EDGE_TRUNCATE)
;ye=DERIV(localtime, t_bot)
;w3 = where(localtime eq 13.) ;;11. flat!
;y(i) = ye[w3]

endfor
print, x2
print, y2
fit  = LINFIT(x ,y )
print, fit
fit2 = LINFIT(x2,y2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
file='./data_hinson' & header='' & nlines_header=0 & ncol = 7
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nlines = FILE_LINES(file)-nlines_header & data=FLTARR(ncol,nlines)
OPENR, lun, file, /GET_LUN & READF, lun, data & CLOSE, lun & FREE_LUN, lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
alt = reform(data(5,*)) - reform(data(3,*))		;; altitude
xdata = (exp(alt/10.))^r_cp				;; p/p0 with assumed scale height
ydata = reform(data(6,*))				;; potential temperature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ydata2 = reform(data(5,*));reform(data(3,*));reform(data(5,*))  				;; altitude of the top of the boundary layer zt
xdata2 = reform(data(6,*))	                        ;; potential temperature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data = 0.
fitdata  = LINFIT(xdata ,ydata )
fitdata2 = LINFIT(xdata2,ydata2)
;fitdata2 = POLY_FIT(xdata2,ydata2,2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
set_name      =  saveplot+'comp_EXN.ps'
set_title     =  ""
set_titlex    =  'Base Exner function P!U-1!N=(p!D0!N/p!Db!N)!UR/cp!N'
set_titley    =  'Base potential temperature T!Db!N (K)'
set_subtitle  =  ''
set_xrange    =  [0.90,1.25]
set_yrange    =  [210.,265.]
;set_xrange    =  [0.90,1.10]
;set_yrange    =  [210.,250.]
	;set_xrange    =  [0.80,1.20]
	;set_yrange    =  [2.,5.]
        ;set_xrange    =  [0.80,1.20]
        ;set_yrange    =  [180.,230.]
set_tickx     =  0.05
set_ticky     =  5.
wherex	      =  1.00
wherey        =  220.	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
altlin=0
dddy = (set_yrange(1) - set_yrange(0))/20.
dddx = (set_xrange(1) - set_xrange(0))/20.

!p.multi=[0,2,1]

set_title='LES results. T!Db!N = '+string(fit(1), '(F5.1)')+' P!U-1!N + '+string(fit(0), '(F5.1)')
                plot, $
set_xrange, fit[1]*set_xrange + fit[0], linestyle=0, $
                xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, $
                ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, $
                 title=set_title, subtitle=set_subtitle, color=0
for i=0, n_elements(paths)-1 do xyouts, x(i), y(i), shorttt(i)
;		xyouts, $
;wherex, wherey, 'squares : data. T!Db!N = '+string(fitdata(1), '(F6.2)')+' P!U-1!N + '+string(fitdata(0), '(F6.2)')
;                xyouts, $
;wherex, wherey-dddy, 'triangles : LES. T!Db!N = '+string(fit(1), '(F6.2)')+' P!U-1!N + '+string(fit(0), '(F6.2)')

set_title='RO data. T!Db!N = '+string(fitdata(1), '(F5.1)')+' P!U-1!N + '+string(fitdata(0), '(F5.1)')
                plot, $
set_xrange, fitdata[1]*set_xrange + fitdata[0], $
                xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, $
                ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, $
                 title=set_title, subtitle=set_subtitle, color=0, linestyle=1
                oplot, $
xdata, ydata, psym=6
;                oplot, $
;x, y, psym=5

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_End, /PNG


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
fit 	= fit2
x	= x2
y	= y2
fitdata	= fitdata2
xdata 	= xdata2 
ydata 	= ydata2  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
set_name      =  saveplot+'comp_ZTOP.ps'
set_title     =  "Absolute altitude of the BL top z!Dt!N (km)"
set_titlex    =  'Base potential temperature T!Db!N (K)'
set_titley    =  'Absolute altitude of the BL top (km)'
set_subtitle  =  ''
set_xrange    =  [150.,280.] ;[200.,265.]
set_yrange    =  [-8.,25.]   ;[-4.,20.]
set_tickx     =  10.
set_ticky     =  2.
wherex        =  205.
wherey        =  17.0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
altlin=0
dddy = (set_yrange(1) - set_yrange(0))/20.
dddx = (set_xrange(1) - set_xrange(0))/20.
;yeye = set_xrange[0] + (set_xrange[1]-set_xrange[0])*findgen(100.)/99.
                plot, $
set_xrange, fitdata[1]*set_xrange + fitdata[0], $
;yeye, fitdata[2]*yeye^2 + fitdata[1]*yeye + fitdata[0], $
                xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, $
                ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, $
                 title=set_title, subtitle=set_subtitle, color=0, linestyle=altlin
                oplot, $
xdata, ydata, psym=6
                oplot, $
x, y, psym=5
                oplot, $
set_xrange, fit[1]*set_xrange + fit[0], linestyle=1
                xyouts, $
wherex, wherey, 'squares : data.  z!Dt!N  =  ('+strtrim(string(fitdata(1)),1)+')  T!Db!N + ('+strtrim(string(fitdata(0)),1)+')'
                xyouts, $
wherex, wherey-dddy, 'triangles : model. z!Dt!N  =  ('+strtrim(string(fit(1)),1)+')  T!Db!N + ('+strtrim(string(fit(0)),1)+')'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_End, /PNG


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


no_therm:

end
