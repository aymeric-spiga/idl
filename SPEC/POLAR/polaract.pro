pro polaract
;
;
;
what_I_plot=0.
overcontour=0.
@polaract_inc.pro
SPAWN, '\rm param_plot.idl ; cp polar_inc.pro param_plot.idl'
;
;
;
  zefile=save_ps
  PS_Start, filename=zefile+'.ps'
  print, zefile+'.ps'
  !P.Charsize = 0.6
  !P.multi=[0,2,1]
	;!p.charthick = 2.0
	;!p.thick = 2.0
	;!x.thick = 2.0
	;!y.thick = 2.0
;
;
;
jojo = ['gw','gwsmooth']
	;jojo = ['gwplus','gwsmooth']
	;jojo = ['gwhires','gw']
filename_sav = filename
for i=0,n_elements(jojo)-1 do begin
;
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS//OUTILS_CONVERSION/'+jojo(i)+'/'+filename_sav
print, filename
;
getcdf, file=filename, charvar=cfield1, invar=field1
getcdf, file=filename, charvar=cfield2, invar=field2
getcdf, file=filename, charvar='XLONG', invar=longi
getcdf, file=filename, charvar='XLAT', invar=lati
;
;
;
if (jojo(i) eq 'gwhires') then begin
	sp = 6;5 ;; relaxation width
	nx = n_elements(longi(0,*))
	ny = n_elements(longi(*,0))
	field1 		= field1        (sp:nx-sp-1,sp:ny-sp-1,*,*)
	field2          = field2        (sp:nx-sp-1,sp:ny-sp-1,*)
	longi           = longi         (sp:nx-sp-1,sp:ny-sp-1,*)
	lati            = lati          (sp:nx-sp-1,sp:ny-sp-1,*)
endif
;
;
;
latmin = -90. & latmax = -50.0 & lonmin = -180. & lonmax = 180.
	;lonmin = 0. & lonmax = 90. ;; pour des tranches
map_set, -90., 0., /isotropic, /azimuthal, /noborder, limit=[latmin,lonmin,latmax,lonmax],title=title_user,/advance ;, position=[0.10, 0.12, 0.90, 0.92]
;
;
;
  ;; average on time and height
dbl_integral = TOTAL(TOTAL(field1*field1,3),3) 
dbl_integral = TEMPORARY(dbl_integral) / n_elements(field1(0,0,*,0)) / n_elements(field1(0,0,0,*))
print, n_elements(field1(0,0,*,0)), n_elements(field1(0,0,0,*))
  ;; smooth a little bit
smoothampl=2 ;0;3
overcontour = smooth(TEMPORARY(dbl_integral),[smoothampl,smoothampl],/EDGE_TRUNCATE)
  ;; plot topography
what_I_plot = reform(field2(*,*,0))
  ;; coordinates
lon = reform(longi(*,*,0))
lat = reform(lati (*,*,0))
overvector_x=0.
overvector_y=0.
;
;
;
SPAWN, '\rm param_plot.idl ; cp polaract_inc.pro param_plot.idl'
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
loadct, 0
MAP_GRID, CHARSIZE = 0.6, $
          COLOR    = 0,   $
          LABEL    = 1,   $   (one label any 1 grid lines)
          LATDEL   = 05., $   
          LONDEL   = 45., $   
	  lats=-60, $
          GLINESTYLE = 1, $ 
          GLINETHICK = 0.1, $ 
          LONLAB   = -70., $ 
          LATLAB   = (lonmin+lonmax)/2.
;
;
;
endfor
PS_End, /PNG
end
