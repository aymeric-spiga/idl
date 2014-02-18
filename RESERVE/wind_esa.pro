pro wind_esa

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filename='ok_mesoscale_profile'
retrieve='true'
;retrieve='false'
file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_68.5/jour20/wrfout_d01_2024-03-20_06:00:00'
start_time=6.
interv=1.
landing_site=[353.87-360.,-1.88]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
;height=total(height,2)/n_elements(height(0,*))
height=height/1000./3.72
;heightp=height(0:n_elements(height(*))-2)-elevation(indp[0],indp[1])/1000.
heightp = height - elevation(indp[0],indp[1]) / 1000.
;
; save
;
save, what_I_plot, localtime, heightp, filename='yeahyeah' 

endif else begin

restore, filename='yeahyeah'	
	
endelse	

lctime=[7.,8.,9.,10.,11.]
for ll=0,n_elements(lctime)-1 do begin
  w = where(localtime eq lctime(ll))
  print, localtime[w(0)]
  tab1 = reform(heightp(*,w(0)))
  tab2 = reform(what_I_plot(*,w(0)))
  tab1 = tab1[0:n_elements(tab1)-14]
  tab2 = tab2[0:n_elements(tab2)-14]
  make_ascii, filename+'_LT'+string(lctime(ll),'(I0)'), tab1, tab2
  plot, tab1, tab2
endfor


end
