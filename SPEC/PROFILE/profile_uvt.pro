pro profile_uvt
;
;
;
what_I_plot=0.
overcontour=0
mention=''
@map_uvt_inc.pro
SPAWN, '\rm param_plot.idl ; cp map_uvt_inc.pro param_plot.idl ; cp -f map_uvt_inc.pro '+save_ps+'.map_uvt_inc.pro'
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
!p.multi=[0,2,2]
!P.Charsize = 1.0
;
;
;
if (n_elements(field1) ne 0) then getcdf, file=filename, charvar=field1, invar=cfield1
u  = getget(filename, 'Um',   count=[1,1,0,1], offset=[zex, zey, 0,ntime])
v  = getget(filename, 'Vm',   count=[1,1,0,1], offset=[zex, zey, 0,ntime])
w  = getget(filename, 'WAVE', count=[1,1,0,1], offset=[zex, zey, 0,ntime]) 
u2 = getget(filename, 'Um',   count=[1,1,0,1], offset=[zex2,zey2,0,ntime])
v2 = getget(filename, 'Vm',   count=[1,1,0,1], offset=[zex2,zey2,0,ntime])
w2 = getget(filename, 'WAVE', count=[1,1,0,1], offset=[zex2,zey2,0,ntime])
tk = getget(filename, 'tk',   count=[1,1,0,1], offset=[zex, zey, 0,ntime])
tk2 = getget(filename, 'tk',   count=[1,1,0,1], offset=[zex2, zey2, 0,ntime])
tpot = getget(filename, 'tpot',   count=[1,1,0,1], offset=[zex, zey, 0,ntime])
tpot2 = getget(filename, 'tpot',   count=[1,1,0,1], offset=[zex2, zey2, 0,ntime])
getcdf, file=filename, charvar='XLONG', invar=longi
getcdf, file=filename, charvar='XLAT', invar=lati
getcdf, file=filename, charvar='HGT', invar=hgt
getcdf, file=filename, charvar='vert', invar=vert
getcdf, file=filename, charvar='TSURF', invar=tsurf
getcdf, file=filename, charvar='USTM', invar=ustar
print, longi(zex ,zey ,0), lati(zex ,zey ,0)
print, longi(zex2,zey2,0), lati(zex2,zey2,0)
;
;
;
loadct, 0
contour, reform(hgt(*,*,ntime)), $
reform(longi(*,*,ntime)), $
reform(lati(*,*,ntime)), $
nlevels=20, $
xtitle='Longitude', $
ytitle='Latitude', $
;xrange=[-146.,-126.], $
;yrange=[11.,27.], $
xrange=[-141.,-133.], $
yrange=[18.,26.], $
xtickinterval=1., $
ytickinterval=1., $
max_value=22000., $
min_value=-4000., $
/cell_fill
xyouts, longi(zex ,zey ,0), lati(zex ,zey ,0), '+ full', color=255, charsize=1.5
xyouts, longi(zex2 ,zey2 ,0), lati(zex2 ,zey2 ,0), '+ dashed', color=255, charsize=1.5
;
;
;
	;plot, tsurf(zex ,zey ,*), yrange=[150, 300]
	;oplot, tsurf(zex2 ,zey2 ,*), linestyle=1
;
;
;
what_I_plot = tk & column = vert
;yeye = tsurf(zex,zey,ntime) & what_I_plot = [yeye,reform(what_I_plot)] & column = [0.01, vert]
overplot = tk2 & overplot_column = vert
;yeye = tsurf(zex2,zey2,ntime) & overplot = [yeye,reform(overplot)] & overplot_column = [0.01, vert]
print, min(what_I_plot), max(what_I_plot)
print, min(overplot), max(overplot)
mention='three'
;
;
;
profile, $
        what_I_plot, $                          ; 1D vertical profile
        column, $                               ; altitudes     
        alt=alt, $                              ; altitude range [altmin, altmax]
        minfield=minfield_init2, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init2, $               ; maximum value of plotted field (=0: calculate)
        inprofile=overplot, $                   ; another vertical profile to overplot
        incolumn=overplot_column, $             ; altitudes of the other vertical profile (in case /= column)
        discrete=discrete, $                    ; show the profile points (= type of points in !psym)
        title_plot=title_user, $                ; title of the plot ('Profile' is default)
        title_axis=title_axis2, $                ; title of the [x,y] axis (['Field','Altitude'] is default)
        mention=mention                         ; add text precision within the plot window (default is nothing or '')
xyouts, 210.,20.,'T!Ds!N = '+string(tsurf(zex,zey,ntime),'(I0)')+' K'
xyouts, 163.,20.,'T!Ds!N = '+string(tsurf(zex2,zey2,ntime),'(I0)')+' K'
;
;
;
what_I_plot = w & column = vert
overplot = w2 & overplot_column = vert
;
;
;
profile, $
        what_I_plot, $                          ; 1D vertical profile
        column, $                               ; altitudes     
        alt=alt, $                              ; altitude range [altmin, altmax]
        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
        inprofile=overplot, $                   ; another vertical profile to overplot
        incolumn=overplot_column, $             ; altitudes of the other vertical profile (in case /= column)
        discrete=discrete, $                    ; show the profile points (= type of points in !psym)
        title_plot=title_user, $                ; title of the plot ('Profile' is default)
        title_axis=title_axis, $                ; title of the [x,y] axis (['Field','Altitude'] is default)
        mention=mention                         ; add text precision within the plot window (default is nothing or '')
;
;
;
	;plot, tpot, vert
	;oplot, tpot2, vert, linestyle=2
	;PS_End, /PNG
	;stop
what_I_plot = u*u + v*v & what_I_plot = sqrt(what_I_plot) & column = vert
overplot = u2*u2 + v2*v2 & overplot = sqrt(overplot) & overplot_column = vert
print, min(what_I_plot), max(what_I_plot)
print, min(overplot), max(overplot)
mention='two'
;
;
;
profile, $
        what_I_plot, $                          ; 1D vertical profile
        column, $                               ; altitudes     
        alt=alt, $                              ; altitude range [altmin, altmax]
        minfield=minfield_init2, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init2, $               ; maximum value of plotted field (=0: calculate)
        inprofile=overplot, $                   ; another vertical profile to overplot
        incolumn=overplot_column, $             ; altitudes of the other vertical profile (in case /= column)
        discrete=discrete, $                    ; show the profile points (= type of points in !psym)
        title_plot=title_user, $                ; title of the plot ('Profile' is default)
        title_axis=title_axis2, $                ; title of the [x,y] axis (['Field','Altitude'] is default)
        mention=mention                         ; add text precision within the plot window (default is nothing or '')
xyouts, 23.5,20.,'u!D*!N = '+string(ustar(zex,zey,ntime),'(F4.1)')+' m s!U-1!N'
xyouts, 5.,20.,'u!D*!N = '+string(ustar(zex2,zey2,ntime),'(F4.1)')+' m s!U-1!N'
;
;
;
PS_End, /PNG
end
