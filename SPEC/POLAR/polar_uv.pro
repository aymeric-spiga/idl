pro polar_uv
;
;
;
what_I_plot=0.
overcontour=0.
@polar_inc_uv.pro
SPAWN, '\rm param_plot.idl ; cp polar_inc_uv.pro param_plot.idl ; cp -f polar_inc_uv.pro '+save_ps+'_polar_inc_uv.pro'
if (n_elements(coord2d) eq 0) then coord2d='true'
if (n_elements(model) eq 0) then model='LMD_MMM'
;
;
;
if (n_elements(field1) ne 0) then getcdf, file=filename, charvar=field1, invar=cfield1
getcdf, file=filename, charvar='Um', invar=u
getcdf, file=filename, charvar='Vm', invar=v
case model of
'LMD_MMM': begin
  getcdf, file=filename, charvar='XLONG', invar=longi
  getcdf, file=filename, charvar='XLAT', invar=lati
  getcdf, file=filename, charvar='HGT', invar=hgt
end
'MRAMS': begin
  getcdf, file=filename, charvar='topo', invar=hgt
  getcdf, file=filename, charvar='glat', invar=lati
  getcdf, file=filename, charvar='glon', invar=longi
end
endcase
if (n_elements(field2) ne 0) then getcdf, file=filename, charvar=field2, invar=cfield2
;
;
;
	;sp = 5 ;; relaxation width
	;nx = n_elements(longi(0,*))
	;ny = n_elements(longi(*,0))
	;if (n_elements(field1) ne 0) then begin
	;	cfield1	= cfield1	(sp:nx-sp-1,sp:ny-sp-1,*,*)
	;endif
	;if (n_elements(field2) ne 0) then begin
	;        cfield2 = cfield2       (sp:nx-sp-1,sp:ny-sp-1,*,*)
	;endif
	;u 		= u		(sp:nx-sp-1,sp:ny-sp-1,*,*) 
	;v       	= v       	(sp:nx-sp-1,sp:ny-sp-1,*,*)
	;longi		= longi		(sp:nx-sp-1,sp:ny-sp-1,*)
	;lati          	= lati 	        (sp:nx-sp-1,sp:ny-sp-1,*)
	;hgt		= hgt		(sp:nx-sp-1,sp:ny-sp-1,*)
	;nx = n_elements(longi(0,*))
	;ny = n_elements(longi(*,0))
;
;
;
;latmin = -90. & latmax = -50.0 & lonmin = -180. & lonmax = 180.
;map_set, -90., 0., /isotropic, /azimuthal, /noborder, limit=[latmin,lonmin,latmax,lonmax],title=title_user,/advance 
;latmin = 75. & latmax = 90.0 & lonmin = -180. & lonmax = 180.
;map_set, 90., 0., /isotropic, /azimuthal, /noborder, limit=[latmin,lonmin,latmax,lonmax],title=title_user,/advance 
;
;
;
		;if (n_elements(ntime) eq 0 or ntime eq 99) then begin
		if (ntime eq 99) then begin
			PRINT, '-- ALL TIME STEPS'
			ntstart = 0
			ntend = n_elements(reform(u(0,0,0,*)))-1 & print, ntend 
		endif else begin
			PRINT, '-- ONLY TIME STEP ', string(ntime,'(I0)')
			ntstart = ntime
			ntend = ntime
		endelse
                for ntime = ntstart,ntend do begin

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
overvector_x = reform(u(*,*,nlevel,ntime))
overvector_y = reform(v(*,*,nlevel,ntime))
if (n_elements(field1) eq 0) then begin
        print, 'field1: horizontal velocity'
	zevel = overvector_x^2 + overvector_y^2  ;; attention il faut que les tableaux soient de la meme taille
	what_I_plot = sqrt(zevel)
endif else begin
	if (no3d eq 'true') then what_I_plot = reform(cfield1(*,*,ntime)) else what_I_plot = reform(cfield1(*,*,nlevel,ntime))
endelse
if (n_elements(field2) eq 0) then overcontour = reform(hgt(*,*,ntime)) else overcontour = reform(cfield2(*,*,ntime))
;
;
;
if (coord2d eq 'false') then begin ;; SI PAS UMET et VMET
	;npoints=n_elements(lon(*,0)) + n_elements(lon(0,*))
        ;what_I_plot = REBIN( what_I_plot, npoints, npoints )
        ;overvector_x = REBIN( overvector_x, npoints, npoints )
	;overvector_y = REBIN( overvector_y, npoints, npoints )
	;overcontour = REBIN( overcontour, npoints, npoints )
        ;lon =  minlon + (maxlon - minlon)*findgen(npoints)/float(npoints-1)
        ;lat =  minlat + (maxlat - minlat)*findgen(npoints)/float(npoints-1)
lon = findgen( n_elements(what_I_plot(*,0)) )
lat = findgen( n_elements(what_I_plot(0,*)) )
print, 'ATTENTION!!!! NOUS SOMMES BIEN D ACCORD QUE VOUS NE CHARGEZ PAS UMET et VMET ??'
endif
if (coord2d eq 'regular') then begin ;; carte avec mercator
lon = reform(lon(*,0))
lat = reform(lat(0,*))
endif
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
;;
;;
;;
;loadct, 0
;contour, reform(field3(*,*,ntime)), $
;	lon,lat, $
;        /overplot, $ 
;	levels=[1,100,200], $
;	c_labels=[0,0,0], $
;        c_thick=1.5, $
;	color=255
;


PS_End, /PNG
		endfor
end



;
; SI MAP_SET EST REGLE....
;
;; marche pas avec les vecteurs


loadct, 0
MAP_GRID, $
          CHARSIZE = 1., $
          COLOR    = 0,   $
        ;  LABEL    = 2,   $   ;; /LABEL or LABEL=2 (one label any 2 grid lines)
	;          LATDEL   = 10., $   ;;5
	;          LONDEL   = 30., $   ;;15
	;lats=-60, $
	;          LONLAB   = -70., $ ;(latmin+latmax)/2., $
	;          LATLAB   = (lonmin+lonmax)/2., $
          LABEL    = 1,   $   ;; /LABEL or LABEL=2 (one label any 2 grid lines)
          LATDEL   = 10., $   ;;5
          LONDEL   = 15., $   ;;15
          LONLAB   = -0., $ ;(latmin+latmax)/2., $
          LATLAB   = -0.001, $
          GLINESTYLE = 2, $
          GLINETHICK = 0.3
	  ;LONALIGN = 0., $
	  ;LATALIGN = 1.

;;
;;
PS_End, /PNG
end
