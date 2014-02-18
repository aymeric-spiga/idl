pro map_static
;
;
;
what_I_plot=0.
overcontour=0.
@map_static_inc.pro
SPAWN, '\rm param_plot.idl ; cp map_static_inc.pro param_plot.idl'
if (n_elements(coord2d) eq 0) then coord2d='false'
;
;
;
  zefile=save_ps
  PS_Start, filename=zefile+'.ps'
  print, zefile+'.ps'
  !P.Charsize = 1.2
  !p.charthick = 2.0
  !p.thick = 2.0
  !x.thick = 2.0
  !y.thick = 2.0
;
;
;
if (n_elements(field1) ne 0) then getcdf, file=filename, charvar=field1, invar=cfield1
getcdf, file=filename, charvar='XLONG_M', invar=longi
getcdf, file=filename, charvar='XLAT_M', invar=lati
getcdf, file=filename, charvar='HGT_M', invar=hgt
getcdf, file=filename, charvar='SLPX', invar=hxwrf
getcdf, file=filename, charvar='SLPY', invar=hywrf

	hx=hxwrf*1000. & hy=hywrf*1000.
	;; SLOPE ANGLE vs HORIZONTAL
	theta=atan(sqrt(hx^2+hy^2)) & theta=180.*theta/!pi ;& theta=round(theta)
;
;
;
        hgt = smooth(hgt, [2,2], /EDGE_TRUNCATE) ;;; truc
;
;
;
sp = 5 ;; relaxation width
nx = n_elements(longi(0,*))
ny = n_elements(longi(*,0))
if (n_elements(field1) ne 0) then begin
	cfield1	= cfield1	(sp:nx-sp-1,sp:ny-sp-1,*,*)
endif else begin
	cfield1 = theta         (sp:nx-sp-1,sp:ny-sp-1,*,*)
endelse
longi		= longi		(sp:nx-sp-1,sp:ny-sp-1,*)
lati           	= lati 	        (sp:nx-sp-1,sp:ny-sp-1,*)
hgt		= hgt		(sp:nx-sp-1,sp:ny-sp-1,*)
nx = n_elements(longi(0,*))
ny = n_elements(longi(*,0))
;
;
;
overcontour = reform(hgt(*,*))
lon = reform(longi(*,*))
lat = reform(lati(*,*))
overvector_x = 0
overvector_y = 0
what_I_plot = reform(cfield1(*,*))
help, what_I_plot, lon, lat
;
;
;
map_latlon, $
        what_I_plot, $                          ; 2D field
        lon, $                                  ; 1D latitude
        lat, $                                  ; 1D longitude
        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
        overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
        overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
        overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
        ct=pal, $                               ; color table (33-rainbow is default)
        colors=colors, $                        ; number of colors/levels (32 is default)
        title=title_user, $                     ; title of the plot ('' is default)
        format=format                           ; format of colorbar annotations ('(F6.2)' is default)
;
;
;
PS_End, /PNG
end
