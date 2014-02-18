pro windmeso

;;Procedure de trace de valeurs issues du Mesoscale en temps-altitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;Choix des champs a tracer, 3 cas :
;;Un champ de deux variables, dont on trace le module, entrer [X,Y]
;;Une variable X seule, dont on trace la valeur entrer  [X,X]
;;tk, entrer [T,PTOT]

;champ1='U' & champ2='V' & champlat = 'XLAT' & champlon = 'XLONG' & champhgt = 'HGT' & file='/home/aymeric/Science/INTERCOMP/LMD_MMM/LMD_MMM_d04_2024-05-57_06:00:00' & nomfich='nest4_MSL' & charphtot = 'PHTOT'

champ1='u_avg' & champ2='v_avg' & champlat = 'glat' & champlon = 'glon' & champhgt = 'topo' & file='/home/aymeric/Science/INTERCOMP/MRAMS/MRAMS_d04_2024-05-57_06:00:00' & nomfich='nest4_MSL_mrams' & charphtot = 'z_lyrmid_agl'

;nomfich=nomfich+'proch'

;champ1='T' & champ2='PTOT' & champlat = 'XLAT' & champlon = 'XLONG' & champhgt = 'HGT' & file='/home/aymeric/Science/INTERCOMP/LMD_MMM/LMD_MMM_d04_2024-05-57_06:00:00' & nomfich='nest4_MSLtemp' & charphtot = 'PHTOT'

;champ1='tempk' & champ2='tempk' & champlat = 'glat' & champlon = 'glon' & champhgt = 'topo' & file='/home/aymeric/Science/INTERCOMP/MRAMS/MRAMS_d04_2024-05-57_06:00:00' & nomfich='nest4_MSL_mramstemp' & charphtot = 'z_lyrmid_agl'


;;Decommenter pour un fichier d'API
;api='on'	;;a tester, voir REF si ca ne marche pas

;;Choix du fichier de simulation
;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_68.5/jour20/wrfout_d01_2024-03-20_06:00:00'
;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau2/wrfout_d01_2024-09-08_00:00:00'
;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau05/wrfout_d01_2024-09-08_00:00:00'
;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau5/wrfout_d01_2024-09-08_00:00:00'
;file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau1/wrfout_d01_2024-09-08_00:00:00'
;file='/tmp7/aslmd/wrfout_d01_2024-09-10_06:00:00'
;file='/tmp7/aslmd/MERIDIANI_FLUSH/wrfout_d01_2024-09-07_06:00:00'
;file='/tmp7/aslmd/MERIDIANI_FLUSH/wrfout_d01_2024-09-13_06:00:00'
;file='/tmp7/aslmd/MERIDIANI_FLUSH2/wrfout_d01_2024-09-07_06:00:00'
;file='/tmp7/aslmd/MERIDIANI_FLUSH2/wrfout_d01_2024-09-13_06:00:00'
;file='/tmp7/aslmd/MERIDIANI_BOMB/wrfout_d01_2024-09-07_06:00:00'
;file='/home/aymeric/Science/INTERCOMP/LMD_MMM/LMD_MMM_d04_2024-05-57_06:00:00'

;;Nom du fichier de sortie
	;nomfich='flush'
        ;nomfich='nest4_MSL'
	;nomfich='flush_late'
	;nomfich='flush2'
	;nomfich='flush2_late'
	;nomfich='bomb'
        nomfich= nomfich+'_vert'
        ;nomfich='LS244_tau1_vent_bas'
	;nomfich='test'
	
;; Choix des options :
;;Heure de depart du trace sur les abscisses (ne fait que decaler le trace : le mettre a l'heure de debut de la simu)
	start_time=0. ;6. ;0.
;;Intervalle horaire (ne fait qu'elargir le trace, en placant une heure de simu tous les interv)
	interv=1. 
	
;;Gestion des sites
n = 1 ;;nombre de sites
sites = make_array(n,/STRING)
coord = make_array(2,n,/FLOAT)
	sites[0] ='landing_site'
	;sites[0] = 'landing_site'
	;sites[1] = 'crater'
	;sites[2] = 'crater_rim'
	;sites[3] = 'high_area'
	;sites[4] = 'low_area'
	;sites[5] = 'plain'
	;sites[] = ''
	;; si on ne veut pas de nom avec les coordonnes rajouter une case avec rien dedans
        coord[0,0] =325.149-360.        & coord[1,0] =-26.372
	;coord[0,0] =353.87-360.	& coord[1,0] =-1.88
	;coord[0,1] =-9.93 		& coord[1,1] =5.45
	;coord[0,2] =-9.78		& coord[1,2] =6.63
	;coord[0,3] =5.77		& coord[1,3] =-6.84
	;coord[0,4] =-11.78		& coord[1,4] =7.80
	;coord[0,5] =-2.67		& coord[1,5] =-1.8
	;coord[0,] =-2.67		& coord[1,] =-1.8
	
;;sites de reference:
;'landing_site'	coord[0,] =353.87-360.		& coord[1,] =-1.88
;'crater'		coord[0,] =-9.93 		& coord[1,] =5.45
;'crater_rim'		coord[0,] =-9.78		& coord[1,] =6.63
;'high_area'		coord[0,] =5.77			& coord[1,] =-6.84
;'low_area'		coord[0,] =-11.78		& coord[1,] =7.80
;'plain'		coord[0,] =-2.67		& coord[1,] =-1.8
;'crater_slope'	coord[0,] =-9.59		& coord[1,] =5.96
;'crater_LS'		coord[0,] =-3.68		& coord[1,] =-3.65
;'crater_LS_rim'	coord[0,] =-4.70		& coord[1,] =-3.65
;'holden' 		coord[0,] =-34.851		& coord[1,] =-26.372

	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	

;;Pour regler les parametres de trace dans autre chose que param_plot qui risaue de mal finir un jour
SPAWN, '\rm param_plot.idl ; cp -f windautom.idl param_plot.idl'	
	
	
;;Recuperation de la latitude, longitude, topographie et pression (pour le calcul de la heuteur)

       	getcdf, file=file, charvar=champlat, invar=lat
	getcdf, file=file, charvar=champlon, invar=lon
	getcdf, file=file, charvar=champhgt, invar=elevation	


	;;Recuperation des variables
	if champ1 eq champ2 then begin
		getcdf, file=file, charvar=champ1, invar=c
	endif  else begin
		getcdf, file=file, charvar=champ1, invar=c1
		getcdf, file=file, charvar=champ2, invar=c2
	endelse
	
	
for i = 0,n-1 do begin

;;Site du trace
	name_site = sites[i]
	landing_site=[ coord[0,i] , coord[1,i] ]
	
;;Recuperation de [i,j] a partir de [lon,lat]
	findxy, $
		x=lon, $
		y=lat, $
       	 	point=landing_site, $
		ind=indp, $
		flag=2

	
if (n_elements(api) eq 0) then begin
	print,'----------- Standard output'
		
;;Recuperation de la hauteur
	getcdf, file=file, charvar=charphtot, invar=ph
        height=reform(ph(indp[0],indp[1],*,*))
                ;;on se debarasse de la position
        height=total(height,2)/n_elements(height(0,*))
                ;;on somme sur le temps et on moyenne: on obtient un tableau de valeurs d'altitudes moyennes
if (charphtot eq 'z_lyrmid_agl') then begin
	heightp=height/1000.
	print, heightp
endif else begin
	height=height/1000./3.72
	heightp=height(0:n_elements(height(*))-2)-elevation(indp[0],indp[1])/1000.
		;;on soustrait la topographie et on obtient une colonne des hauteurs finales
        print, heightp
endelse
;if (champ1 eq 'W') then heightp=height(0:n_elements(height(*))-1)-elevation(indp[0],indp[1])/1000.

		
	;;Utiliser la boucle ci dessous pour imprimer les hauteurs à entrer dans la création d'un API
	;print='on'
	if (n_elements(print) ne 0) then begin
		for i=0,n_elements(heightp)-1 do begin
			print, STRTRIM(heightp(i),2)+','
		endfor
	endif
	
endif else begin
	print,'----------- API output'
	
	;;créer avec API un fichier contenant Um, Vm (ou autres) et vert, dans namelist.api mettre les interp_levels voulus (par ex avec la boucle print ci dessus)
	getcdf, file=file, charvar='vert', invar=height
	heightp=height/1000.		
endelse	
		
		
;;Creation du champ a tracer

if champ1 eq champ2 then begin
	print,'Vous tracez les valeurs de '+champ1
	champi = reform(c(indp[0],indp[1],*,*)) 
	what_I_plot = champi

if (champ1 eq 'W') then what_I_plot =  what_I_plot[0:n_elements(what_I_plot(*,0))-2,*]
	
endif  else begin

if (champ1 eq 'T') and (champ2 eq 'PTOT') then begin 
	print,'Vous tracez la temperature normale'
	T = reform( c1(indp[0],indp[1],*,*) )  
	PTOT = reform( c2(indp[0],indp[1],*,*) ) 
	what_I_plot =(T+220.)*(PTOT/610.)^(192./844.)
	 
endif else begin
	print,'Vous tracez le module de '+champ1+' et '+champ2
	champi1 = reform(c1(indp[0],indp[1],*,*)) 
	champi2 = reform(c2(indp[0],indp[1],*,*)) 
	what_I_plot = sqrt(champi1^2 + champi2^2)

endelse
endelse


;print, 'mean(what_I_plot) ', mean(what_I_plot) 
;print, what_I_plot

;;Gestion du temps
;;definition des abscisses : il va tracer toutes les valeurs du fichier de simu, mais en les decalant a start time
	localtime = float(start_time) + float(interv)*findgen(n_elements(what_I_plot(0,*))) 	
	
;;Nom du fichier de sortie

PS_Start,filename=nomfich+'_'+name_site+'_lon'+string(coord[0,i],'(I0)')+'_lat'+string(coord[1,i],'(I0)')+'.ps'
!P.Charsize = 1.2

help, transpose(what_I_plot), localtime, heightp
	
map_latlon, $
        transpose(what_I_plot), $      ; 2D field
        localtime, $                   ; 1D latitude
        heightp, $                      ; 1D longitude
	overcontour=0		
	

;;Conversion en .png
PS_End, /PNG

endfor


end	
	
	

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;yeah=where(localtime lt 0) & if (yeah(0) ne -1) then localtime[yeah]=24+localtime[yeah]
;;; 24h issue
;localtime=localtime MOD 24 & yeah=where(localtime eq 0) & if (yeah(0) ne -1) then localtime[yeah]=24
;localtime=24.+localtime

;!p.charthick = 2.0
;!p.thick = 3.0
;!x.thick = 1.0
;!y.thick = 2.0

;overcontour=transpose(what_I_plot)     ; another 2D field to overplot with contour lines (=0: no)

;;ATTENTION : a ne JAMAIS, JAMAIS, JAMAIS faire	
	;field=sqrt(c1^2+c2^2) & c1=0. & c2=0.
	;what_I_plot=reform(field(indp[0],indp[1],*,*))
	;what_I_plot = field
