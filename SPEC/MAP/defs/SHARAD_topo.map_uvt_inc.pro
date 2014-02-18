;;;;;*************************************** FOLDER
        folder  = '/donnees/aslmd/MODELES/MESOSCALE/LMD_MM_MARS/mpi_64/WPS/'
        ;coord2d = 'false'       ;; for non-regular projections  ;;CHANGER windowx and windowy
	coord2d = 'polar'       ;; for polar projection with map_set ;; utiliser use_lt=99
        filename        = folder + 'geo_em.d01.nc'
        ini_utc         = 6           ;; cf. name of file
        freq            = 1           ;; cf. 1 output per ... hour
        save_ps         = 'MOLA_topo'

;;;;;*************************************** LVL & FLD & TIME
        nlevel          = 0
        use_lt          = 99               ;; cf. what user wants (99 pour tous)
        field1          = 'HGT_M'            ;; comment to trace horizontal velocity
        ;field2         = 'PTOT'           ;; contour
        ;no3d           = 'true'
        overvector_x    = 0               ;; comment out to get rid of vectors
        overvector_y    = 0               ;; comment out to get rid of vectors
        fieldscoord     = ['XLONG_M','XLAT_M','HGT_M']

;;;;;*************************************** TWEAK VAR
        what_I_plot = what_I_plot / 1000.
        ;print, max(what_I_plot), min(what_I_plot)

;;;;;*************************************** PLOT TITLES
        title_user 	= 	'MOLA Topography' 
        title_axis      =       ['East longitude','North latitude']
        ;subtitle_user   = 'LMD_MM'
        ;;subtitle_user   = subtitle_user + ' / LT = '+string(use_lt,'(I0)')+'h'
        ;;subtitle_user   = subtitle_user + ' / UTC!D0!N = 06:00am' ;; pour les series
        ;subtitle_user   = subtitle_user + ' / Ls = 120!Uo!N'
        ;subtitle_user   = subtitle_user + ' / dx = 10km [S]'

filename        = folder + 'paleo/geo_em.d01.nc'
save_ps         = 'SHARAD_topo'
title_user      = 'SHARAD Paleo-Topography'

;;;;;*************************************** COLOR TABLES
        flag_cb         =       'true'
	poscb		=	0.85
        format          =       '(F4.1)' 
        colors          =       64 
        pal             =       15 ;33             ;; GOOD: 4, 18, 22, 16, 37, 33, 39, 6, 11

;;;;;*************************************** FILL LIMITS
        minfield_init   =       -6.0
        maxfield_init   =       -2.5
        ndiv            =       7

	;;;;;*************************************** LIMIT TRICKS
        ;;;;;********************* must always follow FILL LIMITS
        	missing_value=1.e30
        	lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
        	lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
        	;lim_blank = 0.05 & w=where(abs(what_I_plot) le lim_blank) & if (w[0] ne -1) then what_I_plot[w]=missing_value

;;;;;*************************************** WINDS
        windex          =       30.             ;; DEF: 20.
        stride          =       3.              ;; DEF: 5.

;;;;;*************************************** CONTOUR
	overcontour	=	overcontour/1000.
        lev             =       -10. + 0.5*findgen(30)

;;;;;*************************************** AXIS
        isotropic       =       'true'
        intervalx       =       30.0
        intervaly       =	05.0

;;;;;*************************************** MAP LIMITS
	windowx 	= 	[-180.,180.]        
	windowy 	= 	[70.,90.]

;;;;;***************************************;;;;;
;;;;;***************************************;;;;;
;;;;;***************************************;;;;;
;;;;;***************************************;;;;;

;;;;;*************************************** SETTING TIME (do not modify)
        if (n_elements(windowx) eq 0) then windowx = 0 
        utc_to_lt       = mean(windowx) / 15.    ;; cf. longitude -- LT = UTC + utc_to_lt
        use_utc         = use_lt - utc_to_lt
        zentime         = floor(((24 + use_utc - ini_utc) MOD 24)/freq)        ;; TRUE IDL SUBSCRIPT... 
        if (use_lt ne 99) then ntime = zentime else ntime = 99    ;; ou commenter pour avoir tous les pas de temps
        print, zentime, ntime

;;;;;***************************************
;;;;;***************************************
;;;;;***************************************
;;;;;***************************************

;;;
;;; VERTICAL SECTION LIMITS
;;;
;minalt=-5.             ;; grepSEC
;maxalt=40.             ;; grepSEC
;minspace=0.0           ;; grepSEC
;maxspace=35.0          ;; grepSEC

;;;
;;; METRIC UNITS FOR VERTICAL SECTION
;;;
;factor=10.        ;; grepSEC
;space=space*60. & spacekm='true'       ;; grepSEC
;minspace=minspace*60./factor           ;; grepSEC
;maxspace=maxspace*60./factor           ;; grepSEC
;intervalx=round(intervalx*60./factor)   ;; grepSEC


;;;
;;; TRICKS
;;;

;; pour tracer juste les vecteurs sur un fond uni ou vide
;pal=0                                          ;; 1/4 grepMAP
;what_I_plot(*,*)=what_I_plot(*,*)*0.+0.2       ;; 2/4 grepMAP
;what_I_plot(0,0)=0.                            ;; 3/4 grepMAP
;flag_cb='false'                                        ;; 4/4 grepMAP

;;; truc pour tracer juste les contours et une zone grisee de topo
;pal=0                                           ;; 1/9 grepALL
;w=where(abs(what_I_plot) lt missing_value)     ;; 2/9 grepALL
;what_I_plot[w]=0.                              ;; 3/9 grepALL
;w=where(abs(what_I_plot) gt missing_value)     ;; 4/9 grepALL
;what_I_plot[w]=0.75                            ;; 5/9 grepALL
;w=where(what_I_plot eq 0.)                     ;; 6/9 grepALL
;what_I_plot[w]=missing_value                           ;; 7/9 grepALL
;what_I_plot(0,0)=1. & what_I_plot(1,0)=0.      ;; 8/9 grepALL
;flag_cb='false'                                ;; 9/9 grepALL

