pro mawd

;;;;;;;
experiment='a035'
;experiment='a040'
;experiment='a045'
;;f_user='/donnees/aslmd/SVN/trunk/mesoscale/TMPDIR/GCMINI/mawd.nc'
;;f_user='/tmp7/aslmd/mawd.nc'
f_user='/d5/aslmd/WATERCYCLE/mawd_'+experiment+'.nc'
v_user = 'mtot'
x_user = 'Time'
y_user = 'latitude'
transp = 'yes'
tes = 'no'
;;;;;;;
;f_user='/donnees/aslmd/TES/TES.SeasonalClimatology.nc'
;v_user = 'water'
;x_user = 'time'
;y_user = 'latitude'
;transp = 'yes'
;tes = 'yes'
;;;;;;;


what_I_plot=0.
overcontour=0
@mawd_inc.pro
print, lev
SPAWN, '\rm param_plot.idl ; cp mawd_inc.pro param_plot.idl'

PS_START, file='mawd_'+experiment+'.ps'
;  !P.Charsize = 1.2
;  !p.charthick = 2.0
;  !p.thick = 2.0
;  !x.thick = 2.0
;  !y.thick = 2.0

cdfid = ncdf_open(f_user)

varid=ncdf_varid(cdfid,v_user)
ncdf_varget, cdfid, varid, champ

varid=ncdf_varid(cdfid,y_user)
ncdf_varget, cdfid, varid, yy

varid=ncdf_varid(cdfid,x_user)
ncdf_varget, cdfid, varid, xx

if (tes eq 'no') then begin
  ;;; en precip-microns
  champ = champ * 1.e6 / 917.
  ;;;; entre 0 et 360
  ;xx = xx MOD 360
  xx = xx - 360.*3.
endif else begin
  xx = xx - 360.*2.
endelse

help, champ
if (tes eq 'no') then begin
  what_I_plot = champ
endif else begin
  what_I_plot = float(reform(champ(0,*,*)))
endelse
help, what_I_plot
overcontour = what_I_plot

if (tes eq 'yes') then begin
  w=where(what_I_plot eq -1.)
  what_I_plot[w] = !VALUES.F_NAN
  overcontour[w] = !VALUES.F_NAN
endif

if (transp eq 'yes') then begin
	what_I_plot = transpose(what_I_plot)
	overcontour = transpose(overcontour)
endif

what_I_plot = smooth(what_I_plot,[10,1])

map_latlon, $
        what_I_plot, $                          ; 2D field
        xx, $                                  ; 1D latitude OR 2D
        yy, $                                  ; 1D longitude OR 2D
;        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
;        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
        overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
;        overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
;        overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
        ct=pal, $                               ; color table (33-rainbow is default)
;        colors=colors, $                        ; number of colors/levels (32 is default)
;        title=title_user, $                     ; title of the plot ('' is default)
        format=format                           ; format of colorbar annotations ('(F6.2)' is default)


PS_END, /PNG

end
