pro map_uvt
;
;
;
what_I_plot=0.
overcontour=0.
@map_uvt_inc.pro
SPAWN, '\rm param_plot.idl ; cp map_uvt_inc.pro param_plot.idl'
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
getcdf, file=filename, charvar='Um', invar=u
getcdf, file=filename, charvar='Vm', invar=v
getcdf, file=filename, charvar='XLONG_M', invar=longi
getcdf, file=filename, charvar='XLAT_M', invar=lati
getcdf, file=filename, charvar='HGT_M', invar=hgt
;
;
;
sp = 5 ;; relaxation width
nx = n_elements(longi(0,*))
ny = n_elements(longi(*,0))
if (n_elements(field1) ne 0) then begin
	cfield1	= cfield1	(sp:nx-sp-1,sp:ny-sp-1,*,*)
endif
u 		= u		(sp:nx-sp-1,sp:ny-sp-1,*,*) 
v       	= v       	(sp:nx-sp-1,sp:ny-sp-1,*,*)
longi		= longi		(sp:nx-sp-1,sp:ny-sp-1,*)
lati           	= lati 	        (sp:nx-sp-1,sp:ny-sp-1,*)
hgt		= hgt		(sp:nx-sp-1,sp:ny-sp-1,*)
nx = n_elements(longi(0,*))
ny = n_elements(longi(*,0))
;
;
;
overcontour = reform(hgt(*,*,ntime))
lon = reform(longi(*,*,ntime))
lat = reform(lati(*,*,ntime))
overvector_x = reform(u(*,*,nlevel,ntime))
overvector_y = reform(v(*,*,nlevel,ntime))
if (n_elements(field1) eq 0) then begin
        print, 'field1: horizontal velocity'
	zevel = overvector_x^2 + overvector_y^2  ;; attention il faut que les tableaux soient de la meme taille
	what_I_plot = sqrt(zevel)
endif else begin
;	what_I_plot = reform(cfield1(*,*,nlevel,ntime))
what_I_plot = reform(cfield1(*,*))
endelse
;
;
;
minlat=min(lat) & maxlat=max(lat) & minlon=min(lon) & maxlon=max(lon)
if (coord2d eq 'true') then begin
	npoints=n_elements(lon(*,0)) + n_elements(lon(0,*))  ;; trop de points, mais au moins on ne perd rien 
	TRIANGULATE, lon, lat, tr 
	what_I_plot  = GRIDDATA( lon, lat, what_I_plot,  /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
	overvector_x = GRIDDATA( lon, lat, overvector_x, /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
	overvector_y = GRIDDATA( lon, lat, overvector_y, /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
        overcontour  = GRIDDATA( lon, lat, overcontour,  /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
			; sale sale sale
			if (minlat lt min(lat(*,0))) then overvector_y=-overvector_y
			if (minlon lt min(lon(0,*))) then overvector_x=-overvector_x
        lon =  minlon + (maxlon - minlon)*findgen(npoints)/float(npoints-1)
        lat =  minlat + (maxlat - minlat)*findgen(npoints)/float(npoints-1)
endif else begin
	;npoints=n_elements(lon(*,0)) + n_elements(lon(0,*))
        ;what_I_plot = REBIN( what_I_plot, npoints, npoints )
        ;overvector_x = REBIN( overvector_x, npoints, npoints )
	;overvector_y = REBIN( overvector_y, npoints, npoints )
	;overcontour = REBIN( overcontour, npoints, npoints )
        ;lon =  minlon + (maxlon - minlon)*findgen(npoints)/float(npoints-1)
        ;lat =  minlat + (maxlat - minlat)*findgen(npoints)/float(npoints-1)
endelse

;lon = reform(lon(*,0))
;lat = reform(lat(0,*))

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
