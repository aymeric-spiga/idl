pro wind_analysis

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
retrieve='true'
;retrieve='false'
file='zefolder/wrfout_d01_2024-03-21_00:00:00'
start_time=0.
interv=1.
landing_site=[353.87-360.,-1.88]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
; graphics definition
;
PS_Start, filename='wind.ps'
!P.Charsize = 1.2

;!p.charthick = 2.0
;!p.thick = 3.0
;!x.thick = 2.0
;!y.thick = 2.0

if (retrieve eq 'true') then begin

;
; load WRF coordinates
;
getcdf, file=file, charvar='XLAT', invar=lat
getcdf, file=file, charvar='XLONG', invar=lon
getcdf, file=file, charvar='HGT', invar=elevation
;
; find simulated lander grid point
;
findxy, $
	x=lon, $
	y=lat, $
        point=landing_site, $
	ind=indp, $
	flag=2
;
; get field
;
getcdf, file=file, charvar='U', invar=u
getcdf, file=file, charvar='V', invar=v
;;field=sqrt(u^2+v^2) & u=0. & v=0.
;;what_I_plot=reform(field(indp[0],indp[1],*,*))
  ;;; BEWARE: arrays not the same size... because staggering... u^2+v^2 is wrong and IDL does not complain...
field=sqrt( u(indp[0],indp[1],*,*)^2 + v(indp[0],indp[1],*,*)^2 ) & u=0. & v=0.
what_I_plot=reform(field)

;
; local time
;
localtime = float(start_time) + float(interv)*findgen(n_elements(field(0,0,0,*))) ;+ lon(indp[0],indp[1])/15. 
;
; get model height
;
getcdf, file=file, charvar='PHTOT', invar=ph
height=reform(ph(indp[0],indp[1],*,*))
height=total(height,2)/n_elements(height(0,*))
height=height/1000./3.72
heightp=height(0:n_elements(height(*))-2)-elevation(indp[0],indp[1])/1000.
;
; save
;
save, what_I_plot, localtime, heightp, filename='yeahyeah' 

endif else begin

restore, filename='yeahyeah'	
	
endelse	

;yeah=where(localtime lt 0) & if (yeah(0) ne -1) then localtime[yeah]=24+localtime[yeah]
;;; 24h issue
;localtime=localtime MOD 24 & yeah=where(localtime eq 0) & if (yeah(0) ne -1) then localtime[yeah]=24
;localtime=24.+localtime

map_latlon, $
        transpose(what_I_plot), $      ; 2D field
        localtime, $                   ; 1D latitude
        heightp, $                      ; 1D longitude
	overcontour=0
;        overcontour=transpose(what_I_plot)     ; another 2D field to overplot with contour lines (=0: no)

PS_End, /PNG

end
