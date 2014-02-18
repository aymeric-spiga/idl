pro script_con

;;;;;;;
file_user='/home/aslmd/EXOMARS/LMD_MM_MARS_Meridiani_S114_Ls247_LT0-23.nc'
var_user = 'TSURF'
var_user2 = 'HGT'
nt = 5
format = '(I0)'
pal=33
;;;;;;;

PS_START, file='toto.ps'


cdfid = ncdf_open(file_user)
varid=ncdf_varid(cdfid,var_user)
ncdf_varget, cdfid, varid, champ

cdfid = ncdf_open(file_user)
varid=ncdf_varid(cdfid,var_user2)
ncdf_varget, cdfid, varid, champ2

cdfid = ncdf_open(file_user)
varid=ncdf_varid(cdfid,'XLONG')
ncdf_varget, cdfid, varid, lon

cdfid = ncdf_open(file_user)
varid=ncdf_varid(cdfid,'XLAT')
ncdf_varget, cdfid, varid, lat



nx = 25
ny = 25
what_I_plot = reform(champ(nx,ny,*))

plot, what_I_plot

PS_END, /PNG

end



what_I_plot = reform(champ(*,*,nt))
overcontour = reform(champ2(*,*,nt))
longi = reform(lon(*,0))
lati = reform(lat(0,*))

map_latlon, $
        what_I_plot, $                          ; 2D field
        longi, $                                  ; 1D latitude OR 2D
        lati, $                                  ; 1D longitude OR 2D
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
