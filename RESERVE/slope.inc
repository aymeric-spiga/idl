;
; horizontal coordinates (meters)
;
rad=3397000.
deg_m=2*rad*!pi/360.
dlon=deg_m*cos(lat*!pi/180)
dlat=deg_m
x=(lon-min(lon))*dlon
y=(lat-min(lat))*dlat
;print, x/1000.
;print, y/1000.

;;----------------------------
;; methode directe
;;----------------------------
;; var is a 2D field
;; must be HGT_M
grad, hx, var*1000., x, y, 1
grad, hy, var*1000., x, y, 2
;hx=hx*100. ;pente en %

;;-------------------------------
;; methode "confiance en geogrid"
;;-------------------------------
	;varid=ncdf_varid(cdfid,'MAPFAC_M')
        ;	ncdf_varget, cdfid, varid, msf


varid=ncdf_varid(cdfid,'SLPX')
	ncdf_varget, cdfid, varid, hxwrf
varid=ncdf_varid(cdfid,'SLPY')
	ncdf_varget, cdfid, varid, hywrf

;print, 100.*mean(abs(hx-hxwrf*1000.)/hx)
;print, 100.*mean(abs(hy-hywrf*1000.)/hy)

hx=hxwrf*1000.
hy=hywrf*1000.



;; SLOPE ANGLE vs HORIZONTAL
;;---------------------------	
theta=atan(sqrt(hx^2+hy^2))
theta=180.*theta/!pi
var=round(theta)
title='Slope angle (deg)'
format='(I0)'
pal=16

;goto, squeeze
;; SLOPE ORIENTATION
;;--------------------
;; cas où hx < 0
	psi=-!pi/2 - atan(hy/hx)
;; cas où hx > 0
	w=where(hx ge 0)
	psi[w]=psi[w]-!pi
;; SOUTH 0 & WEST 90
	psi=2*!pi+psi
	psi=180.*psi/!pi
;; NORTH 0 & EAST 90 (as MONTECARLO & PARAM)
psi=(180.+psi) MOD 360

var=round(psi)
title='Slope orientation (deg)'
format='(I0)'
pal=6
;;
squeeze:


;
; latitude
;
varid=ncdf_varid(cdfid,'XLAT_M')
        ncdf_varget, cdfid, varid, latyeah


;
; SCATTERED FLUX RATIO
;

nx=n_elements(latyeah(*,0))
ny=n_elements(latyeah(0,*))

ratio_s=dblarr(nx,ny)
ratio_d=dblarr(nx,ny)

count=0

for i=0,nx-1 do begin
for j=0,ny-1 do begin

	param_slope, $
        	ls=set_ls, $            ; < input *
        	lct=set_lct, $          ; < input *
        	lat=latyeah(i,j), $     ; < input *
        	tau=set_tau, $         	; < input *
        	spz=theta(i,j), $       ; < input *
        	spa=psi(i,j), $         ; < input *
        	r_s=tmps, $     ; > output
        	r_d=tmpd     	; > output

	;count=count+1
	;print, count, tmps, tmpd

	ratio_s(i,j)=tmps
	ratio_d(i,j)=tmpd

endfor
endfor

pal=0	
format='(I0)'
var=ratio_s*100. & title='Ratio S/S0 (%)'

;var=ratio_d*100. & title='Ratio Fdir/Fdir0 (%)'

title='Ls='+string(set_ls,'(I0)')+' LT='+string(set_lct,'(I0)')+' tau='+string(set_tau,'(F3.1)')
minfield_init=30.
maxfield_init=170.

var=round(var)
