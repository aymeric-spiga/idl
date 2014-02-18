pro plotplot


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
file='hi_res_run_swordfish/3830.002/m4/diagfi.nc'
set_name='psurf'
it=0
iz=1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

SPAWN, 'cp psurf_param_plot.idl param_plot.idl'

	;
	; graphics definition
	;
	if (n_elements(saveps) eq 0) then saveps='true'
	if (saveps eq 'false') then begin 
	   ;!p.multi=[0,3,2] 
	   !P.CHARSIZE=2.
	   WINDOW, /PIXMAP & WDELETE & DEVICE,BYPASS_TRANSLATION=0,DECOMPOSED=0,RETAIN=2
	endif else begin
	   PREF_SET, 'IDL_PATH', '/padata/beta/users/aspiga/Save/SOURCES/IDL/fsc_psconfig:<IDL_DEFAULT>', /COMMIT
	endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;
getcdf, file=file, charvar='ps', invar=ps
getcdf, file=file, charvar='u'    , invar=u
getcdf, file=file, charvar='v'    , invar=v
getcdf, file=file, charvar='lon'  , invar=lon
getcdf, file=file, charvar='lat'  , invar=lat
;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;
what_I_plot = reform(ps(*,*,it))
ox          = reform(u(*,*,iz,it))
oy          = reform(v(*,*,iz,it))
;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;
	if (saveps eq 'true') then PS_Start, FILENAME=set_name+'.ps'
	!P.Charsize = 1.2
map_latlon, what_I_plot, lon, lat, overvector_x=ox, overvector_y=oy
	if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;
cmd='cp plotplot.pro '+set_name+'_plotplot.pro ; cp param_plot.idl '+set_name+'_param_plot.idl'
SPAWN, cmd
;;;;;;;;;;;;;;;;;;;;;;;;;;;


end
