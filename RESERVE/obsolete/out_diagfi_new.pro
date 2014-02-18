pro out_diagfi, $
		plot=plot_user, $		; 'map', 'profile', 'hovmuller', 'meridional', 'zonal' (default is 'map')
	;----------------
	; input fields
	;----------------
		field1=field_user, $		; field to contour: 'U', 'V', 'tk', ... etc ...
		field2=fieldc_user, $		; field to overplot (contour lines, optional)
	;------------------------------------------------------
	; 'map' options (no effect when plot is 'profile')
	;------------------------------------------------------
		topo=topo, $			; overplot topography contours (!! must load HGT !! override field2)
		winds=winds, $			; option to overplot winds (e.g. =['U','V'])
		level=level, $			; vertical level subscript (starts at 1 ...)(if -1: all levels)
	;------------------------------------------------------
	; 'profile' options (no effect when plot is 'map')
	;------------------------------------------------------
		inprofile=inprofile, $		; a profile to compare with ... (override field2)
		incolumn=incolumn, $		; altitudes ... in case /= column
		;; NB: inprofile and incolumn may be used to output profiles
		;; ... reset with inprofile=0, incolumn=0
	;----------------
	; choose time
	;----------------
		when=when, $			; time subscript (starts at 1 ...)(if -1: first time step)
	;----------------
	; choose latitude/longitude
	;----------------
		nlat=nlat, $			; latitude subscript (starts at 1 ...)(if -1: all latitudes)(default is middle)
		nlon=nlon, $			; longitude subscript (starts at 1 ...)(if -1: all longitudes)(default is middle)
	;----------------
	; plot options
	;----------------
		range=range_user, $		; range=[min,max]
		colors=colors, $		; 16, 32, 64, ... etc
	;-------------------------
	; no useless data loading
	;-------------------------
		save_data=save_data		; to avoid loading data if already done once



;---------------------------------------------------------------------------
;
; A general routine to plot ARWpost outputs
;
; USE:
; - out_wrf	(will plot everything !)
; - out_wrf, field1='tk', range=[160,300]
;
; - out_wrf, field1='pressure', field2='tk'
; - out_wrf, field1=TSK, field2=HGT, winds=['U','V'] 
;
; - out_wrf, save_data=field
;   ... then next calls including 'save_data=field' won't include the data loading 
;
; - out_wrf, field1='tk', winds=['umet','vmet'], topo=1, save_data=yeah, level=1, when=2
;
; - out_wrf, field1='tk', plot='profile'
;
; - out_wrf, plot='profile', field1='tk', when=1, inprofile=inprofile, incolumn=incolumn
;   out_diagfi, plot='profile', field1='temp', when=112, inprofile=inprofile, incolumn=incolumn 	
;
; - out_wrf, plot='meridional',field1='U'
;
; - out_wrf, plot='hovmuller',field1='TSK'
;
; - out_wrf, plot='zonal',field1='V', nlat=-1
; 
;---------------------------------------------------------------------------
; A. Spiga, August 2007, adding plot keyword September 2007
;---------------------------------------------------------------------------


    


;-------------------------------------
; default settings and user settings
;-------------------------------------
diff_scale=0  ;;1
charend=''
format=''
if (n_elements(plot_user) eq 0) then plot_user='map'
if (n_elements(field_user) eq 0) then field_user='all'
if (n_elements(level) eq 0) then level=-1
if (n_elements(when) eq 0) then when=-1
;;if (n_elements(nlon) eq 0) then nlon=-1
;;if (n_elements(nlat) eq 0) then nlat=-1
if (n_elements(colors) eq 0) then colors=32
if (n_elements(minfield_init)*n_elements(maxfield_init) eq 0) then begin
	minfield_init = 0.
	maxfield_init = 0.
endif

if (n_elements(range_user) eq 2) then begin
	minfield_init=float(range_user(0))
	maxfield_init=float(range_user(1))
	print, 'setting user limits !', minfield_init, maxfield_init 
	diff_scale=1
endif

if (n_elements(topo) ne 0) then fieldc_user='HGT'	;case with topography



;******************************************************************************
;
; BEGIN: different from out_wrf
;
;******************************************************************************

;TEMP TEMP
vert_coord='sigma'
;vert_coord='pressure'
;vert_coord='height'
;TEMP TEMP

;----------------------
; set parameters
;----------------------
file='input.diagfi'
if (n_elements(topo) ne 0) then fieldc_user='phisinit'	;case with topography



;----------------------------------------------------
; Open diagfi file GCM and get coordinates & field
;----------------------------------------------------
cdfid = ncdf_open(file)

print, 'read coord ...'
varid=ncdf_varid(cdfid,'latitude')
        ncdf_varget, cdfid, varid, lat_gcm
varid=ncdf_varid(cdfid,'longitude')
        ncdf_varget, cdfid, varid, lon_gcm
varid=ncdf_varid(cdfid,'Time')
        ncdf_varget, cdfid, varid, hour
lon=lon_gcm
west_east=n_elements(lon)
lat=lat_gcm
south_north=n_elements(lat)
print, 'done !'
;;print, 'time ...', hour*24 MOD 24, ' hours and day number ', ceil(hour)

; field
if (n_elements(save_data) eq 0) then begin
	print, 'read data ...'
	varid=ncdf_varid(cdfid,field_user)
		ncdf_varget, cdfid, varid, var
		bottom_top=n_elements(var(0,0,*,0))
		time=n_elements(var(0,0,0,*))
;save_array=fltarr(3,west_east,south_north,bottom_top,time)
;save_array(0,*,*,*,*)=var(*,*,*,*)
;;**NB: ralentit beaucoup
	print, 'done !'
endif else begin	
	var=reform(save_data(0,*,*,*,*))
		bottom_top=n_elements(var(0,0,*,0))
		time=n_elements(var(0,0,0,*))
endelse

print, 'vertical levels in file ... ', bottom_top
print, 'time steps in file ...', time



	;-----------------------------
	; min/max on the whole field
	;-----------------------------
	w=where(abs(var) lt 1e30)
	if (w(0) ne -1) then begin
		minfield=min(var[w])
		maxfield=max(var[w])
			;;; manage unvalid points
			;;;w=where(var eq 1e37)
			;;;if (w(0) ne -1) then var[w]=!Values.F_NAN
	endif else begin
		var=0
		print, 'no valid values in this field', charvar 
	endelse



;---------------------
; Get altitude
;---------------------
if (vert_coord ne 'pressure') then begin
	print, 'read altitude ...'
	varid=ncdf_varid(cdfid,'aps')
        	ncdf_varget, cdfid, varid, aps
	varid=ncdf_varid(cdfid,'bps')
        	ncdf_varget, cdfid, varid, bps
	varid=ncdf_varid(cdfid,'ps')
        	ncdf_varget, cdfid, varid, ps
	varid=ncdf_varid(cdfid,'altitude')
        	ncdf_varget, cdfid, varid, z
	print, 'done !'
endif else begin
	print, 'read altitude ...'
	varid=ncdf_varid(cdfid,'altitude')
        	ncdf_varget, cdfid, varid, z
	z=z/100.   ;mbar
	print, 'done !'
endelse


;-------------
; Get winds
;-------------
if (n_elements(winds) ne 0) then begin
if (n_elements(save_data) eq 0) then begin
	print, 'read winds ...'
	varid=ncdf_varid(cdfid,winds(0))
		ncdf_varget, cdfid, varid, u
	varid=ncdf_varid(cdfid,winds(1))
		ncdf_varget, cdfid, varid, v
	;save_array(1,*,*,*,*)=u(*,*,*,*)
	;save_array(2,*,*,*,*)=v(*,*,*,*)
	print, 'done !'
endif else begin	
	u=reform(save_data(1,*,*,*,*))
	v=reform(save_data(2,*,*,*,*))
endelse
endif

;-------------------------------
; Get another field to contour
;-------------------------------
if (n_elements(fieldc_user) ne 0) then begin
	print, 'get alternate field ...'
	varid=ncdf_varid(cdfid,fieldc_user)
		ncdf_varget, cdfid, varid, varcontour
	print, 'done !'
endif


;;-------------------------------
;; Save data for further use ...
;;-------------------------------
;if (n_elements(save_data) eq 0) then begin
;	save_data=save_array
;	save_array=0. ;free 'save_array' array
;endif


;---------------------------------
; Some stuff for continuity
;---------------------------------
;;if (vert_coord ne 'pressure') then surfpres=reform(ps(*,*,nt))
;;A CORRIGER - A METTRE DANS LA BOUCLE
if (vert_coord ne 'pressure') then surfpres=reform(ps(*,*,0))
if (vert_coord eq 'sigma') then begin
	z=aps/surfpres(0,0)+bps
	z=z*100.
	vert_coord='model_level'
endif
charend='_diagfi'

case field_user of 
'temp': charvar='tk'
'u': charvar='U'
'v': charvar='V'
'ps': charvar='PSFC'
else:
endcase

;;
;;


;******************************************************************************
;
; END: different from out_wrf
;
;******************************************************************************



;----------------------
; initialize graphics
;----------------------

!p.charthick = 2.0
!p.thick = 3.0
!x.thick = 2.0
!y.thick = 2.0

if (when ne -1) then charend=charend+string(when, '(I0)')

set_plot, 'ps'
device, filename=plot_user+'_'+vert_coord+'_'+charvar+charend+'.ps', /color

case charvar of
'HGT': begin
        title='Topography (km)'
	pal=16
	level=1
	format='(F7.0)'
end
'tk': begin
	title='Temperature (K)'
	pal=3
	format='(F4.0)'
end
'TSK': begin
	title='Surface Temperature (K)'
	pal=3
	level=1
end
'U': begin
	title='Zonal wind (!Nm!N.s!U-1!N)'
	pal=33
end
'V': begin
	title='Meridional wind (!Nm!N.s!U-1!N)'
	pal=33
end
'pressure': begin
	title='Pressure (hPa)'
	pal=16
end
'height': begin
	title='Height (km)'
	pal=16
end
'T': title='Potential temperature departure from t0 (K)'
'W': begin
	title='Vertical wind (!Nm!N.s!U-1!N)'
	pal=16
	format='(F6.2)'
end
else: begin
	title=''
	pal=33
end
endcase
title_save=title



;---------------------------------
; choose positioning ...
;---------------------------------

case plot_user of
'map': begin
;
; 1. map: looping on vertical levels
;
	if (level eq -1) then begin
		theloopmin=0
		theloopmax=bottom_top-1
	endif else begin
		theloopmin=level-1 ;!! IDL arrays starts at 0
		theloopmax=level-1
	endelse	
	if (when eq -1) then nt = 0 else nt=when-1
end
'profile': begin
;
; 2. profile: looping on time
;
	;;discrete=4
	title_axis=strarr(2)
	title_axis(0)=title_save
	if (when eq -1) then begin
		theloopmin=0
		theloopmax=time-1
	endif else begin
		theloopmin=when-1 ;!! IDL arrays starts at 0
		theloopmax=when-1
	endelse	

	column=z
	case vert_coord of
		'height': title_axis(1)='Altitude (km)'
		'pressure': begin
			H=10.
			p0=6.1
			column=H*alog(p0/z)
			title_axis(1)='Log-Pressure altitude (km)'
		end
		'model_level': title_axis(1)='Model level'
	endcase

	displaylat=float(round(lat(nlat)*100))/100.
	displaylon=float(round(lon(nlon)*100))/100.
	description='WRF Profile at latitude '+string(displaylat,format='(1F6.2)')+$
	'!Uo'+'!N and longitude '+string(displaylon,format='(1F7.2)')+'!Uo!N '
	title=description
end
'hovmuller': begin
;
; 3. hovmuller: looping on vertical levels
;
	if (level eq -1) then begin
		theloopmin=0
		theloopmax=bottom_top-1
	endif else begin
		theloopmin=level-1 ;!! IDL arrays starts at 0
		theloopmax=level-1
	endelse	
end
'zonal': begin
;
; 4. zonal: looping on latitude
;
	if (nlat eq -1) then begin
		theloopmin=0
		theloopmax=south_north-1
	endif else begin
		theloopmin=nlat-1 ;!! IDL arrays starts at 0
		theloopmax=nlat-1
	endelse	
	if (when eq -1) then nt = 0 else nt=when-1
	case vert_coord of
		'height': charlev='Altitude (km)'
		'pressure': begin
			H=10.
			p0=6.1
			z=H*alog(p0/z)
			charlev='Log-Pressure altitude (km)'
		end
		'model_level': charlev='Model level'
	endcase
end
'meridional': begin
;
; 5. meridional: looping on longitude
;
	if (nlon eq -1) then begin
		theloopmin=0
		theloopmax=west_east-1
	endif else begin
		theloopmin=nlon-1 ;!! IDL arrays starts at 0
		theloopmax=nlon-1
	endelse	
	if (when eq -1) then nt = 0 else nt=when-1
	case vert_coord of
		'height': charlev='Altitude (km)'
		'pressure': begin
			H=10.
			p0=6.1
			z=H*alog(p0/z)
			charlev='Log-Pressure altitude (km)'
		end
		'model_level': charlev='Model level'
	endcase
end
else: stop
endcase


;-------------
; plot loop
;-------------	
for theloop=theloopmin,theloopmax do begin

case plot_user of
'map': begin
	;
	; case 1: map
	;	
	case vert_coord of
	'height': charlev=string(round(z(theloop)*1000.),'(I0)')+' m'
	'pressure': charlev=string(round(z(theloop)*100.),'(I0)')+' Pa'
	'model_level': charlev=string(round(z(theloop)),'(I0)')
	endcase
	title=title_save+' at level '+charlev
	print, title

	;-------------
	; draw map
	;-------------

	what_I_plot=var[*,*,theloop,nt]

	if (diff_scale eq 0) then begin
		minfield_init=minfield	;-same scale for all levels 
		maxfield_init=maxfield	;-same scale for all levels
	endif
	if (n_elements(fieldc_user) eq 0) then begin
		overcontour=0.
	endif else begin
		if (n_elements(topo) ne 0) then begin
			overcontour=varcontour[*,*,0,0] ; other levels than 0 are blank !
		endif else begin
			overcontour=varcontour[*,*,theloop,nt]
		endelse
	endelse
	if (n_elements(winds) eq 0) then begin
		overvector_x=0.
		overvector_y=0.
	endif else begin
		overvector_x=u[*,*,theloop,nt]
		overvector_y=v[*,*,theloop,nt]


;--------------------------------------------------------
; NB: you can use the winds fields to combine variables 
;--------------------------------------------------------
if (n_elements(winds) eq 2) then begin
	;** HYDROSTATIC REDUCTION TO LEVEL OF REFERENCE
	if ((winds(0) eq 'tk') and (winds(1) eq 'HGT')) then begin
	overvector_x=reform(u[*,*,6,nt])	;;T 1km
	overvector_y=reform(v[*,*,0,nt])	;;topo (m)
	zref=mean(overvector_y)
	;zref=median(overvector_y)
	;zref=(max(overvector_y)-min(overvector_y))/2.
	;zref=3000.
	print, 'zref', zref
	H=192.*overvector_x/3.72	;unit is m
	what_I_plot=what_I_plot*exp((overvector_y-zref)/H)
	overvector_x=0
	overvector_y=0
	endif
endif
;--------------------------------------------------------
; NB: you can use the winds fields to combine variables
;--------------------------------------------------------

	endelse

;;******************TEMP
;;pressure anomaly
;nday=3
;oneday=reform(var[*,*,theloop,0:nday-1])
;what_I_plot=what_I_plot-(total(oneday,3)/float(nday))
;;******************TEMP

	map_latlon, $
		what_I_plot, $
		lon, $
		lat, $
		ct=pal, $
		minfield=minfield_init, $
		maxfield=maxfield_init, $
		overcontour=overcontour, $
		overvector_x=overvector_x, $
		overvector_y=overvector_y, $
		colors=colors, $
		title=title, $
		format=format	
end
'profile': begin
	;
	; case 2: profile
	;
	print, title

	what_I_plot=var[nlon,nlat,*,theloop]
	if (diff_scale eq 0) then begin
		minfield_init=minfield	;-same scale for all levels 
		maxfield_init=maxfield	;-same scale for all levels
	endif
	if (n_elements(fieldc_user) eq 0) then begin
		overcontour=0.
		overcontourz=0.
	endif else begin
		overcontour=varcontour[nlon,nlat,*,theloop]
		overcontourz=0.
	endelse
	if (n_elements(inprofile) gt 1) then begin
		overcontour=inprofile
		overcontourz=0.
		if (n_elements(incolumn) gt 1) then overcontourz=incolumn
	endif

	profile, $
		what_I_plot, $				; 1D vertical profile
		column, $				; altitudes	
	;	alt=alt, $				; altitude range [altmin, altmax]
		minfield=minfield_init, $		; minimum value of plotted field (=0: calculate)
		maxfield=maxfield_init, $		; maximum value of plotted field (=0: calculate)
		inprofile=overcontour, $		; another vertical profile to overplot
		incolumn=overcontourz, $		; altitudes of the other vertical profile (in case /= column)
		title_plot=title, $			; title of the plot ('Profile' is default)
		title_axis=title_axis, $		; title of the [x,y] axis (['Field','Altitude'] is default)
	;	mention=mention, $			; add text precision within the plot window (default is nothing or '')
		discrete=discrete			; show the profile points (= type of points in !psym)
end
'hovmuller': begin
	;
	; case 3: hovmuller
	;
	; NB: lat/lon hovmuller 
	;
	what_I_plot=reform(var[*,nlat,theloop,*])
	space=lon
	temps=findgen(time)*2
	if (n_elements(fieldc_user) ne 0) then begin
		overcontour=reform(varcontour[*,nlat,theloop,*])
	endif else begin
		overcontour=0.
	endelse
	minspace=0.
	maxspace=0.

;;******************TEMP
;	minspace=-90.
;	maxspace=90.

if (n_elements(winds) ne 0) then begin		;;;winds sert à faire des calculs !!
overvector_x=reform(u[*,nlat,4,*])	;;T 1km
overvector_y=reform(v[*,nlat,0,*])	;;topo (m)
zref=median(overvector_y)
H=192.*overvector_x/3.72	;unit is m
;print, (overvector_y-zref)/H
what_I_plot=what_I_plot*exp((overvector_y-zref)/H)


endif 
;oneday=reform(var[*,nlat,0,0:11])
;slice=(total(oneday,2)/12.)
;for i=0, n_elements(var[*,0,0,0])-1 do begin
;	what_I_plot(i,*)=what_I_plot(i,*)-slice(i)
;endfor
;;******************TEMP





	hovmuller, $
		what_I_plot, $				; 2D field
		space, $				; space coordinate
		temps, $				; time coordinate
		minfield=minfield_init, $		; minimum value of plotted field (=0: calculate)
		maxfield=maxfield_init, $		; maximum value of plotted field (=0: calculate)
		minspace=minspace, $			; minimum value of space window (=0: calculate)
		maxspace=maxspace, $			; maximum value of space window (=0: calculate)
		overcontour=overcontour, $		; another 2D field to overplot with contour lines (=0: no)
		ct=pal, $				; color table (33-rainbow is default)
		colors=colors, $			; number of colors/levels (32 is default)
		title=title				; title of the plot ('' is default)
end
'zonal': begin
	;
	; case 4: zonal section
	; 
	displaylat=float(round(lat(theloop)*100))/100.
	title=title_save+ $
		' section at latitude '+string(displaylat,format='(1F6.2)')+ $
		'!Uo!N '
	print, title

	what_I_plot=reform(var[*,theloop,*,nt])
	space=lon
	altitude=z

	if (n_elements(fieldc_user) ne 0) then begin
		overcontour=reform(varcontour[*,theloop,*,nt])
	endif else begin
		overcontour=0.
	endelse

	if (n_elements(winds) eq 0) then begin
		overvector_x=0.
		overvector_y=0.
	endif else begin
		overvector_x=reform(u[*,theloop,*,nt])		;;horizontal component
		overvector_y=reform(v[*,theloop,*,nt])		;;vertical component
	endelse

	minspace=0.
	maxspace=0.

	title_user=title
	title_axis=['Longitude',charlev]

	section, $
		what_I_plot, $				; 2D field
		space, $				; horizontal coordinate
		altitude, $				; altitude coordinate
		minfield=minfield_init, $		; minimum value of plotted field (=0: calculate)
		maxfield=maxfield_init, $		; maximum value of plotted field (=0: calculate)
		minspace=minspace, $			; minimum value of space window (=0: calculate)
		maxspace=maxspace, $			; maximum value of space window (=0: calculate)
		overcontour=overcontour, $		; another 2D field to overplot with contour lines (=0: no)
		overvector_x=overvector_x, $		; wind vector - x component (=0: no)
		overvector_y=overvector_y, $		; wind vector - y component (=0: no)
		colors=colors, $			; number of colors/levels (32 is default)
		title_plot=title_user, $		; title of the plot ('Profile' is default)
		title_axis=title_axis, $		; title of the [x,y] axis (['Field','Altitude'] is default)
		ct=pal, $				; color table (33-rainbow is default)
		format=format
end	
'meridional': begin
	;
	; case 5: meridional section
	; 
	displaylon=float(round(lon(theloop)*100))/100.
	title=title_save+ $
		' section at longitude '+string(displaylon,format='(1F7.2)')+ $
		'!Uo!N '
	print, title

	what_I_plot=reform(var[theloop,*,*,nt])
	space=lat
	altitude=z

	if (n_elements(fieldc_user) ne 0) then begin
		overcontour=reform(varcontour[theloop,*,*,nt])
	endif else begin
		overcontour=0.
	endelse

	if (n_elements(winds) eq 0) then begin
		overvector_x=0.
		overvector_y=0.
	endif else begin
		overvector_x=reform(u[theloop,*,*,nt])		;;horizontal component
		overvector_y=reform(v[theloop,*,*,nt])		;;vertical component
	endelse

	minspace=0.
	maxspace=0.

	title_user=title
	title_axis=['Latitude',charlev]

	section, $
		what_I_plot, $				; 2D field
		space, $				; horizontal coordinate
		altitude, $				; altitude coordinate
		minfield=minfield_init, $		; minimum value of plotted field (=0: calculate)
		maxfield=maxfield_init, $		; maximum value of plotted field (=0: calculate)
		minspace=minspace, $			; minimum value of space window (=0: calculate)
		maxspace=maxspace, $			; maximum value of space window (=0: calculate)
		overcontour=overcontour, $		; another 2D field to overplot with contour lines (=0: no)
		overvector_x=overvector_x, $		; wind vector - x component (=0: no)
		overvector_y=overvector_y, $		; wind vector - y component (=0: no)
		colors=colors, $			; number of colors/levels (32 is default)
		title_plot=title_user, $		; title of the plot ('Profile' is default)
		title_axis=title_axis, $		; title of the [x,y] axis (['Field','Altitude'] is default)
		ct=pal					; color table (33-rainbow is default)
end	








else: 
endcase



endfor
device, /close

	; save the profiles for another plot (last of the time loop)
	if (n_elements(inprofile) le 1) and (plot_user eq 'profile') then begin
		inprofile=what_I_plot
		incolumn=column
	endif

;******************************************************************************
;
; BEGIN: different from out_wrf
;
;******************************************************************************
;;;skip:
;;;endfor

end
