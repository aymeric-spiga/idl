pro hodo

;;Procedure pour tracer des hodographes des vents ou des vents effectifs
	;;On trace en un point donne les vecteurs vents au cours du temps


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;Choix du fichier
	;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_68.5/jour20/wrfout_d01_2024-03-20_06:00:00'
	;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau2/wrfout_d01_2024-09-08_00:00:00'
	;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau05/wrfout_d01_2024-09-08_00:00:00'
	file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau5/wrfout_d01_2024-09-08_00:00:00'
	;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau1/wrfout_d01_2024-09-08_00:00:00'
	
;;Choix de la methode
	effective='true'
	;effective='false'
	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;Reglages

;;Nom du fichier
	nomfich='Hodo_landingsite_tau5'
	
;;Choix des heures (0 = 1ere heure du ncdf)
	start_hour=0.
	end_hour=23.	
	
	
;;Choix de l'eta-level (premier a 1)
	level=2
	
;;Choix de la position (en indices)
	pos=[72.,72.]
	
;;Choix de l'affichage des heures
	debut=0	;;heure reelle (a adapter a start_hour)
	pas=2	;;commenter si on ne veut rien
	
;;Sous titre (commentable)	
	subtitre='Tau = 5'	
	
;;Echelle (commenter pour echelle automatique)
	hrange=[-20.,15.]
	xdiv=7
	vrange=[-30.,5.]
	ydiv=7
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

getcdf, file=file, charvar='XLAT', invar=lat
getcdf, file=file, charvar='XLONG', invar=lon
	lon=reform(lon(*,0,0))		;;pour un [i,j] donne l'abscisse de lon donne la longitude
	lat=reform(lat(0,*,0))

if (effective eq 'false') then begin

	getcdf, file=file, charvar='U', invar=u
	getcdf, file=file, charvar='V', invar=v

	u=reform(u(*,*,level-1,start_hour:end_hour))
	v=reform(v(*,*,level-1,start_hour:end_hour))
	
	var=''
	lev='_lev'+STRTRIM(level,2)


endif else begin

	;;Vent effectif obtenu par integration sur l'altitude
	effwind, $
		'U', $	
		file, $		
		u, $	
		start_time=start_hour, $	
		end_time=end_hour		

	effwind, $
		'V', $	
		file, $		
		v, $
		start_time=start_hour, $
		end_time=end_hour

	var='Effective'
	lev=''
	
endelse
	
	u=reform(u(pos[0],pos[1],*))
	v=reform(v(pos[0],pos[1],*))
	

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PS_Start,filename=nomfich+var+'_pos'+STRTRIM(pos[0],2)+'-'+STRTRIM(pos[1],2)+lev+'.ps',xsize=8,ysize=8
	!P.Charsize = 1.2

;;Reglages graphiques

if (n_elements(hrange) eq 0.) then begin	
	if ( abs( max(u)-min(u)) ge  abs( max(v)-min(v))) then range=abs( max(u)-min(u)) else range=abs( max(v)-min(v))
	hrange=[min(u)-1,min(u)-1+ range +2] 
	vrange=[min(v)-1,min(v)-1 + range +2]
endif
	
	valon=lon(pos[0])
	valat=lat(pos[1])
	
;;Trace	
	plot, $
		u, $
		v, $
		xtitle=var+' U Wind Velocity (m/s)', $
		ytitle=var+' V Wind Velocity (m/s)', $
		title='Hodograph : lon '+string(valon, FORMAT='(F7.2)')+' | lat '+string(valat, FORMAT='(F7.2)'), $
		subtitle=subtitre, $
		PSYM=-4, $
		Xrange=hrange, $
		Xticks=xdiv, $
		Yrange=vrange, $
		Yticks=ydiv, $
		XTicklen=1.0, $
	 	YTicklen=1.0, $
		XGridStyle=1, $
		YGridStyle=1

	;;affichage de l'origine	
	oplot,[0],[0],PSYM=5,symsize=1.5

;;Reperage des heures 
if n_elements(pas) ne 0 then begin

	for t=0,end_hour-start_hour,pas do begin
	
	XYOUTS, u(t),v(t),debut+t,CHARSIZE=0.9
	
	endfor
	 
endif
	
;;Conversion en .png
	PS_End, /PNG

end
