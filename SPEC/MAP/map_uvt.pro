pro map_uvt
;
;
;
what_I_plot=0.
overcontour=0
no3d='true'
@map_uvt_inc.pro
SPAWN, '\rm param_plot.idl ; cp map_uvt_inc.pro param_plot.idl ; cp -f map_uvt_inc.pro '+save_ps+'.map_uvt_inc.pro'
if (n_elements(coord2d) eq 0) then coord2d='false'
if (n_elements(overvector_x) eq 0) then overvector_x = 1  ;; initialisation (!! map_uvt_inc.pro utilise dans map_uvt et map_latlon) 
;
;
;
if (n_elements(field1) ne 0) then getcdf, file=filename, charvar=field1, invar=cfield1
if (n_elements(field2) ne 0) then getcdf, file=filename, charvar=field2, invar=cfield2
if ( ( overvector_x ne 0 ) or ( n_elements(field1) eq 0 ) ) then begin
   u = getget(filename, 'Um', count=[0,0,1,0], offset=[0,0,nlevel,0])
   v = getget(filename, 'Vm', count=[0,0,1,0], offset=[0,0,nlevel,0])
   help, u
endif
if (n_elements(fieldscoord) eq 0) then begin
 getcdf, file=filename, charvar='XLONG', invar=longi
 getcdf, file=filename, charvar='XLAT', invar=lati
 getcdf, file=filename, charvar='HGT', invar=hgt
endif else begin
 print, 'coordinates are ', fieldscoord
 getcdf, file=filename, charvar=fieldscoord(0), invar=longi
 getcdf, file=filename, charvar=fieldscoord(1), invar=lati
 getcdf, file=filename, charvar=fieldscoord(2), invar=hgt
endelse

;
;
;
sp = 5 ;; relaxation width
nx = n_elements(longi(0,*))
ny = n_elements(longi(*,0))
if (n_elements(field1) ne 0) then begin
	if (no3d ne 'true') then cfield1 = cfield1 (sp:nx-sp-1,sp:ny-sp-1,*,*) else cfield1 = cfield1 (sp:nx-sp-1,sp:ny-sp-1,*)
endif
if (n_elements(field2) ne 0) then begin
        if (no3d ne 'true') then cfield2 = cfield2 (sp:nx-sp-1,sp:ny-sp-1,*,*) else cfield2 = cfield2 (sp:nx-sp-1,sp:ny-sp-1,*)
endif
if ( ( overvector_x ne 0 ) or ( n_elements(field1) eq 0 ) ) then begin
  u 		= reform (u		(sp:nx-sp-1,sp:ny-sp-1,0,*) )
  v       	= reform (v       	(sp:nx-sp-1,sp:ny-sp-1,0,*) )
endif
longi		= longi		(sp:nx-sp-1,sp:ny-sp-1,*)
lati          	= lati 	        (sp:nx-sp-1,sp:ny-sp-1,*)
hgt		= hgt		(sp:nx-sp-1,sp:ny-sp-1,*)
nx = n_elements(longi(*,0,0))
ny = n_elements(longi(0,*,0))
nt = n_elements(longi(0,0,*))
if ( overvector_x ne 0 ) then begin  
  overvector_x = u
  overvector_y = v
endif
;
;
;
                if (ntime eq 99) then begin
                        PRINT, '-- ALL TIME STEPS', nt-1 
                        ntstart = 0  & ntend = nt-1 
                endif else begin
                        PRINT, '-- ONLY TIME STEP ', string(ntime,'(I0)')
                        ntstart = ntime & ntend = ntime
                endelse
                for ntime = ntstart,ntend do begin
help, u
;
;
;
  zefile=save_ps+string(100+ntime,'(I0)')
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
lon = reform(longi(*,*,ntime))
lat = reform(lati(*,*,ntime))
if (n_elements(field1) eq 0) then begin
        print, 'field1: horizontal velocity'
		;;;attention attention overvector_x et overvector_y sont vides
		;zevel = overvector_x^2 + overvector_y^2  ;; attention il faut que les tableaux soient de la meme taille, OK si uvmet
        zevel = (reform(u(*,*,ntime)))^2 + (reform(v(*,*,ntime)))^2 ;; attention il faut que les tableaux soient de la meme taille, OK si uvmet
	what_I_plot = sqrt(zevel)
endif else begin
	if (no3d ne 'true') then what_I_plot = reform(cfield1(*,*,nlevel,ntime)) else what_I_plot = reform(cfield1(*,*,ntime))
        if (n_elements(u) ne 0) then overvector_x = reform(u(*,*,ntime))   ;; ne pas utiliser test overvector_x a cause de la boucle temps
        if (n_elements(v) ne 0) then overvector_y = reform(v(*,*,ntime))   ;; ne pas utiliser test overvector_y a cause de la boucle temps
endelse
if (n_elements(field2) eq 0) then begin
        print, 'field2: topography'
        overcontour = reform(hgt(*,*,ntime))
endif else begin
        if (no3d ne 'true') then overcontour = reform(cfield2(*,*,nlevel,ntime)) else overcontour = reform(cfield2(*,*,ntime))
endelse
;
;
;
if ( coord2d eq 'polar' ) then begin
   print, 'OK YOU USE MAP_SET with POLAR PROJECTION. VECTORS ARE NOT SUPPORTED. USE polar_uv OR ADAPT THIS SCRIPT.'
   overvector_x = 0 
   overvector_y = 0
   if (n_elements(windowx) ne 0) then begin
       latmin = windowy(0) & latmax = windowy(1) & lonmin = windowx(0) & lonmax = windowx(1)
   endif else begin
       latmin = 65. & latmax = 90. & lonmin = -180. & lonmax = 180.
   endelse
   print, 'latmin,lonmin,latmax,lonmax', latmin,lonmin,latmax,lonmax
   map_set, latmax, 0., /isotropic, /azimuthal, /noborder, limit=[latmin,lonmin,latmax,lonmax], title=title_user, /advance
endif else begin
  minlat=min(lat) & maxlat=max(lat) & minlon=min(lon) & maxlon=max(lon)
  if (coord2d eq 'true') then begin
        ;;;; CECI EST DESORMAIS FAIT DANS MAP_LATLON  
	;;npoints=n_elements(lon(*,0)) + n_elements(lon(0,*))  ;; trop de points, mais au moins on ne perd rien 
	;;TRIANGULATE, lon, lat, tr 
	;;what_I_plot  = GRIDDATA( lon, lat, what_I_plot,  /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
        ;;if ( overvector_x ne 0 ) then begin
	;;  overvector_x = GRIDDATA( lon, lat, overvector_x, /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
	;;  overvector_y = GRIDDATA( lon, lat, overvector_y, /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
        ;;endif
        ;;overcontour  = GRIDDATA( lon, lat, overcontour,  /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
        ;;                if ( overvector_x ne 0 ) then begin
	;;		  ; sale sale sale
	;;		  if (minlat lt min(lat(*,0))) then overvector_y=-overvector_y
	;;		  if (minlon lt min(lon(0,*))) then overvector_x=-overvector_x
        ;;                endif
        ;;lon =  minlon + (maxlon - minlon)*findgen(npoints)/float(npoints-1)
        ;;lat =  minlat + (maxlat - minlat)*findgen(npoints)/float(npoints-1)
  endif else begin
	lon = reform(lon(*,0))
	lat = reform(lat(0,*))
	;;npoints=n_elements(lon(*,0)) + n_elements(lon(0,*))
        ;;what_I_plot = REBIN( what_I_plot, npoints, npoints )
        ;;overvector_x = REBIN( overvector_x, npoints, npoints )
	;;overvector_y = REBIN( overvector_y, npoints, npoints )
	;;overcontour = REBIN( overcontour, npoints, npoints )
        ;;lon =  minlon + (maxlon - minlon)*findgen(npoints)/float(npoints-1)
        ;;lat =  minlat + (maxlat - minlat)*findgen(npoints)/float(npoints-1)
		;lon = findgen( n_elements(what_I_plot(*,0)) )
		;lat = findgen( n_elements(what_I_plot(0,*)) )
		;print, 'ATTENTION!!!! NOUS SOMMES BIEN D ACCORD QUE VOUS NE CHARGEZ PAS UMET et VMET ??'
  endelse
                        ;;;;; trouve dans polar_uv.pro
			;;if (coord2d eq 'regular') then begin ;; carte avec mercator
			;;	lon = reform(lon(*,0))
			;;	lat = reform(lat(0,*))
			;;endif
endelse
help, what_I_plot, lon, lat
print, 'COLOR PLOT min/max ', min(what_I_plot), max(what_I_plot)
print, 'CONTOUR PLOT min/max ', min(overcontour), max(overcontour)
print, 'VECTOR PLOT min/max ', min(overvector_x), min(overvector_y)
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
if ( coord2d eq 'polar' ) then begin    ;;; pourrait aller dans map_latlon
  loadct, 0
  MAP_GRID, $
          CHARSIZE = 1., $
          COLOR    = 0,   $
          ;lats=-60, $
          LABEL    = 1,   $   ;; /LABEL or LABEL=2 (one label any 2 grid lines)
          LATDEL   = intervaly, $   ;;5 10
          LONDEL   = intervalx, $   ;;15
          ;LONLAB   = latmin + intervaly/2., $ ;5. + (latmin+latmax)/2., $ ;0.
          LONLAB   = (latmin+latmax)/2., $
          LATLAB   = -0.001, $
          GLINESTYLE = 2, $
          GLINETHICK = 0.3
          ;LONALIGN = 0., $
          ;LATALIGN = 1.
endif
;
;
;
PS_End, /PNG
		endfor
end
