pro turb_period
;
;
;
what_I_plot=0.
overcontour=0.
SPAWN, '\rm param_plot.idl ; cp turb_inc.pro param_plot.idl'
;
;
;
altiuser=0.5
start_lt = 8
zefile='turbh_period'
;
;
;
filefile = "/home/aymeric/Work/submitted/mars_journal/LMD_MM_MARS_ustar_LEScaseC.nc"
filefile = "/home/aymeric/extract.nc"
getcdf, file=filefile, charvar="USTM", invar=zevar
getcdf, file=filefile, charvar="U", invar=zeu
getcdf, file=filefile, charvar="V", invar=zev
getcdf, file=filefile, charvar="W", invar=zew
getcdf, file=filefile, charvar="PSFC", invar=zevar

        nx = n_elements(zevar(*,0))
        ny = n_elements(zevar(0,*))
        ;lon = findgen(nx) / 10.
        ;lat = findgen(ny) / 10.
	lon = findgen(2*nx+1) / 10.
	lat = findgen(2*ny+1) / 10.

alpha=1
nstart=195
nend=215
nstride=1
	;nstart=1
	;nend=n_elements(zevar(0,0,*))-1
	;nstride=10
for nt=nstart,nend,nstride do begin
print, nt
;
;
;
  PS_Start, filename=zefile+string(1000+alpha,'(I0)')+'.ps'
  print, zefile+'.ps'
  ;!p.multi=[0,3,2]
  ;!P.Charsize = 0.6
	!P.Charsize = 1.2
	!p.charthick = 2.0
	!p.thick = 2.0
	!x.thick = 2.0
	!y.thick = 2.0
;;
;;
;;
;restore, filename='../LES/MERIDIANI_tau100/savew'
;;restore, filename='../LES/MERIDIANI_tau100_wind10/savew'
;;restore, filename='../LES/MERIDIANI_tau100_wind20/savew'
;;;restore, filename='saveu'
;indind=where(abs(h - altiuser) eq min(abs(h - altiuser)))
;print, h(indind)
;;
;;
;;
;;what_I_plot = reform(wprime(*,*,indind))
;what_I_plot = wprime(*,*,indind)
;;what_I_plot = veltot(*,*,indind)
;;contour, wprime(*,*,indind), nlevels=30
;lon = findgen(n_elements(wprime(*,0,0))) / 10.
;lat = findgen(n_elements(wprime(0,*,0))) / 10.
;nx = n_elements(wprime(*,0,0))
;ny = n_elements(wprime(0,*,0))

	;;;
	what_I_plot = reform(zevar(*,*,nt))
        overvector_x = reform(zeu(0:nx-1,*,0,nt))	
	overvector_y = reform(zev(*,0:ny-1,0,nt))
	what_I_plot = reform(zew(*,*,0,nt))
	;what_I_plot = reform(zevar(*,*,nt)) - mean(reform(zevar(*,*,nt)))
	;;;

what_I_plot_period = fltarr(2*nx+1,2*ny+1)
what_I_plot_period(0:nx-1   ,0:ny-1   ) = what_I_plot(0:nx-1,0:ny-1)
what_I_plot_period(nx:2*nx-1,0:ny-1   ) = what_I_plot(0:nx-1,0:ny-1)
what_I_plot_period(0:nx-1   ,ny:2*ny-1) = what_I_plot(0:nx-1,0:ny-1)
what_I_plot_period(nx:2*nx-1,ny:2*ny-1) = what_I_plot(0:nx-1,0:ny-1)
what_I_plot_period(2*nx,*) = what_I_plot_period(0,*)
what_I_plot_period(*,2*ny) = what_I_plot_period(*,0)
what_I_plot = TEMPORARY(what_I_plot_period)

overvector_x_period = fltarr(2*nx+1,2*ny+1)
overvector_x_period(0:nx-1   ,0:ny-1   ) = overvector_x(0:nx-1,0:ny-1)
overvector_x_period(nx:2*nx-1,0:ny-1   ) = overvector_x(0:nx-1,0:ny-1)
overvector_x_period(0:nx-1   ,ny:2*ny-1) = overvector_x(0:nx-1,0:ny-1)
overvector_x_period(nx:2*nx-1,ny:2*ny-1) = overvector_x(0:nx-1,0:ny-1)
overvector_x_period(2*nx,*) = overvector_x_period(0,*)
overvector_x_period(*,2*ny) = overvector_x_period(*,0)
overvector_x = TEMPORARY(overvector_x_period)

overvector_y_period = fltarr(2*nx+1,2*ny+1)
overvector_y_period(0:nx-1   ,0:ny-1   ) = overvector_y(0:nx-1,0:ny-1)
overvector_y_period(nx:2*nx-1,0:ny-1   ) = overvector_y(0:nx-1,0:ny-1)
overvector_y_period(0:nx-1   ,ny:2*ny-1) = overvector_y(0:nx-1,0:ny-1)
overvector_y_period(nx:2*nx-1,ny:2*ny-1) = overvector_y(0:nx-1,0:ny-1)
overvector_y_period(2*nx,*) = overvector_y_period(0,*)
overvector_y_period(*,2*ny) = overvector_y_period(*,0)
overvector_y = TEMPORARY(overvector_y_period)


print, 'plot !!!!'
;
;
;
	;localtime_h = start_lt + 100*nt/3700  ;; division euclidienne
	;localtime_m = float(100*nt mod 3700)/3700.
	;localtime_m = 60.*localtime_m
	;title_user=string(localtime_h,'(I0)')+'h'+string(localtime_m,'(I0)')
title_user='local time = '+string(start_lt + 100.*float(nt)/3700.,'(F5.2)')
print, title_user
;
;
;
map_latlon, $
        what_I_plot, $                          ; 2D field
        lon, $                                  ; 1D latitude
        lat, $                                  ; 1D longitude
        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
;        overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
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
alpha=alpha+1
endfor
end
