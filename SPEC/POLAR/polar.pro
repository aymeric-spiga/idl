pro polar
;
;
;
what_I_plot=0.
overcontour=0.
@polar_inc.pro
SPAWN, '\rm param_plot.idl ; cp polar_inc.pro param_plot.idl'
;
;
;
  ;zefile=save_ps
  ;PS_Start, filename=zefile+'.ps'
  ;print, zefile+'.ps'
  ;;!p.multi=[0,3,2]
  ;;!p.multi=[0,3,1] 
  ;!p.multi=[0,2,2]
  ;!P.Charsize = 0.6
;
;
;
        jojo = ['gw','gwsmooth']
        ;jojo = ['gwhires','gw']
        jojo = ['gw/30km']
filename_sav = filename
for i=0,n_elements(jojo)-1 do begin
;
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS//OUTILS_CONVERSION/'+jojo(i)+'/'+filename_sav
print, filename
;
getcdf, file=filename, charvar=cfield1, invar=field1
getcdf, file=filename, charvar=cfield2, invar=field2
getcdf, file=filename, charvar=cfield3, invar=field3
getcdf, file=filename, charvar='XLONG', invar=longi
getcdf, file=filename, charvar='XLAT', invar=lati
;
;
;

for ntime = 0, n_elements(field1(0,0,0,*))-1 do begin
;  for ntime = 11, n_elements(field1(0,0,0,*))-1, 5 do begin
;for ntime = 16, n_elements(field1(0,0,0,*))-1, 5 do begin

;for ntime = 06, n_elements(field1(0,0,0,*))-1, 5 do begin

	;;; MOVIE MOVIE MOVIE
	zefile=save_ps+string(100+ntime,'(I0)')
	PS_Start, filename=zefile+'.ps'
	;!P.Charsize = 1.2
	;!p.charthick = 2.0
	;!p.thick = 2.0
	;!x.thick = 2.0
	;!y.thick = 2.0
;
;
;
latmin = -90. & latmax = -50.01 & lonmin = -180. & lonmax = 180.
map_set, -90., 0., /isotropic, /azimuthal, /noborder, limit=[latmin,lonmin,latmax,lonmax],title=title_user,/advance ;, position=[0.10, 0.12, 0.90, 0.92]
;
;
;
what_I_plot = reform(field1(*,*,nlevel,ntime))
overcontour = reform(field2(*,*,ntime))
lon = reform(longi(*,*,ntime))
lat = reform(lati (*,*,ntime))
overvector_x=0.
overvector_y=0.
;
;
;
	;loadct, 0
	;w=where(overcontour lt -5000.) & overcontour[w]=-5000.
	;contour, overcontour, lon, lat, /cell_fill, nlevels=32, /overplot ;& overcontour=0
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
loadct, 0
contour, reform(field3(*,*,ntime)), $
	lon,lat, $
        /overplot, $ 
	levels=[1,100,200,300], $
	c_labels=[0,1,0,1], $
        c_thick=1.5, $ 
	;c_linestyle=2,$
	color=0;255
;
;
;
loadct, 0
MAP_GRID, CHARSIZE = 0.6, $
          COLOR    = 0,   $
          LABEL    = 1,   $   ;(one label any 1 grid lines)
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
PS_End, /PNG
 endfor
  endfor
;PS_End, /PNG
end
