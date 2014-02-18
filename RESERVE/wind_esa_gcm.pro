pro wind_esa_gcm

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filename='gcm_profile'
filetruc='/d6/vblmd/MERIDIANI_EXOMARS/etude_GCM/wrfinput_d01_2024-03-20_'
start_time=7.
landing_site=[353.87-360.,-1.88]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


for i=1,4 do begin

if (start_time lt 10) then yeye='0' else yeye=''
file = filetruc + yeye + string(start_time,'(I0)') + ':00:00'

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
  ;;; it is actually possible to do better with staggering... use XLAT_V and XLONG_U
field=sqrt( u(indp[0],indp[1],*)^2 + v(indp[0],indp[1],*)^2 ) & u=0. & v=0.
what_I_plot=reform(field)
;
; get model height
;
getcdf, file=file, charvar='PHTOT', invar=ph
height=reform(ph(indp[0],indp[1],*))
height=height/1000./3.72
heightp = height - elevation(indp[0],indp[1]) / 1000.

        ;getcdf, file=file, charvar='PTOT', invar=p
        ;getcdf, file=file, charvar='PSFC', invar=psfc
        ;hh = reform(p(0,0,*,0))
        ;p0=610. & t0=220. & r_cp=1/4.4 & grav=3.72 & R=192.
        ;getcdf, file=file, charvar='T', invar=tpot
        ;tpot = tpot + t0
        ;t = tpot * ( p / p0 ) ^ r_cp
       ;hh(0) = - R * t(indp[0],indp[1],0) / grav / 1000. * alog(psfc(indp[0],indp[1]) / p(indp[0],indp[1],0))
        ;for i=1,n_elements(hh)-1 do hh(i)    = hh(i-1) - R * t(indp[0],indp[1],i) / grav / 1000. * alog(p(indp[0],indp[1],i) / p(indp[0],indp[1],i-1))
        ;;heightp = hh
        ;print, heightp-hh  ;; PAREIL ! pb de density global au modele

;
; ASCII table
;
tab1 = heightp
tab2 = what_I_plot
tab1 = tab1[0:n_elements(tab1)-14]
tab2 = tab2[0:n_elements(tab2)-14]
make_ascii, filename+'_LT'+string(start_time,'(I0)'), tab1, tab2
plot, tab1, tab2
;
;
;
start_time = start_time+1
print, start_time

endfor

end
