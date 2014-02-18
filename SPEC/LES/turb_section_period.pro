pro turb_section_period
;
;
;
what_I_plot=0.
overcontour=0.
SPAWN, '\rm param_plot.idl ; cp turb_section_period_inc.pro param_plot.idl'
;
;
;
filefile = "./MERIDIANI_tau050_wind10_section_plume.nc"
getcdf, file=filefile, charvar="W", invar=zew
getcdf, file=filefile, charvar="PHTOT", invar=zevar
zew   = reform(zew)
zevar = reform(zevar)

        nx = n_elements(zevar(*,0))
        nz = n_elements(zevar(0,*))

        xx = findgen(nx) * 100. / 1000.
        zz = reform(zevar(0,*)) / 3.72 / 1000. + ( 1473.35 / 1000. )

	zew = zew[*,0:nz-2]
	zz = zz[0:nz-2]


	PS_Start, filename='sectionLES.ps'
        !P.Charsize = 1.2
        !p.charthick = 2.0
        !p.thick = 2.0
        !x.thick = 2.0
        !y.thick = 2.0

	what_I_plot = zew
	lon = xx
	lat = zz

	map_latlon, $
	        what_I_plot, $                          ; 2D field
	        lon, $                                  ; 1D latitude OR 2D
	        lat, $                                  ; 1D longitude OR 2D
	        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
	        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
	        overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
	        overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
	        overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
	        ct=pal, $                               ; color table (33-rainbow is default)
	        colors=colors, $                        ; number of colors/levels (32 is default)
	        title=title_user, $                     ; title of the plot ('' is default)
	        format=format                           ; format of colorbar annotations ('(F6.2)' is default)

PS_End, /PNG
end
