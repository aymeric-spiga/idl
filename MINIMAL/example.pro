
pro example

name='example'
PS_Start, filename=name+'.ps'
!P.Charsize = 1.2

restore, filename='example.idl'

map_latlon, $
            what_I_plot, $   			  ; 2D field
            x, $                                  ; 1D latitude
            y, $                                  ; 1D longitude
            minfield=0., $                	  ; minimum value of plotted field (=0: calculate)
            maxfield=0., $               	  ; maximum value of plotted field (=0: calculate)
            overcontour=overcontour, $            ; another 2D field to overplot with contour lines (=0: no)
            overvector_x=overvector_x, $          ; wind vector - x component (=0: no)
            overvector_y=overvector_y, $          ; wind vector - y component (=0: no)
            ct=pal, $                             ; color table (33-rainbow is default)
            colors=colors, $                      ; number of colors/levels (32 is default)
            title=title, $                        ; title of the plot ('' is default)
            format=format                         ; format of colorbar annotations ('(F6.2)' is default)

PS_End, /PNG

end
