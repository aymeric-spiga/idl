;;;;;*************************************** FOLDER
	folder  = '/d5/aslmd/LMD_MM_MARS_SIMUS/OM/'
        coord2d = 'regular'       ;; for non-regular projections  ;;CHANGER windowx and windowy
	filename        = folder + 'OM6_TI85/wrfout_d01_2024-06-43_06:00:00_zabg'
	ini_utc         = 6           ;; cf. name of file
	freq            = 1           ;; cf. 1 output per ... hour
	save_ps         = 'OM6_TI85_tsurf'

;;;;;*************************************** LVL & FLD & TIME
	nlevel          = 0 
	use_lt          = 99 ;12           ;; cf. what user wants (99 pour tous)
	field1          = 'TSURF'       ;; comment to trace horizontal velocity
	;field2  	= 'XLAT' 	;; contour
	no3d            = 'true'
	;overvector_x    = 0             ;; comment out to get rid of vectors
	;overvector_y    = 0             ;; comment out to get rid of vectors

;;;;;*************************************** TWEAK VAR
        ;what_I_plot = what_I_plot * 0.
        ;print, max(what_I_plot), min(what_I_plot)

;;;;;*************************************** PLOT TITLES
        title_user    = 'Horizontal wind (m s!U-1!N) 10m AGL'
	title_user    = 'Horizontal wind vectors 20m AGL'
	;title_user    = 'Horizontal wind (m s!U-1!N) 20km AMR'
	title_user	= 'Surface temperature (K)'
	subtitle_user   = 'LMD_MM'
	;subtitle_user   = subtitle_user + ' / LT = '+string(use_lt,'(I0)')+'h'
	subtitle_user   = subtitle_user + ' / UTC!D0!N = 06:00am' ;; pour les series
	subtitle_user   = subtitle_user + ' / Ls = 173!Uo!N'
	subtitle_user   = subtitle_user + ' / dx = 6km [S]'

;;;;;*************************************** COLOR TABLES
        ;flag_cb         =       'false'
        poscb           =       0.96
        format          =       '(I0)'
        ;format          =       '(F4.1)' 
        ;colors          =       128
        colors          =       32
        pal             =       22              ;; GOOD: 4, 18, 22, 16, 37, 33, 39, 6, 11, 19, 0

;;;;;*************************************** FILL LIMITS
        minfield_init   =      0;150.
        maxfield_init   =      0;180.

        ;;;;;*************************************** LIMIT TRICKS
        ;;;;;********************* must always follow FILL LIMITS
                ;missing_value=1.e30
                ;lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
                ;lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
                ;;lim_blank = 2. & w=where(abs(what_I_plot) le lim_blank) & if (w[0] ne -1) then what_I_plot[w]=missing_value

;;;;;*************************************** WINDS
        windex          =       30.            ;; DEF: 20.
        ;windex          =       20.             ;; DEF: 20.
        ;stride          =       1.             ;; DEF: 5.
        ;stride          =       2.              ;; DEF: 5.
        stride          =       5.              ;; DEF: 5.

;;;;;*************************************** CONTOUR
        overcontour     =       overcontour/1000.
        ;lev             =       50. + 50.*findgen(20)
        ;lev             =       -10. + 0.2*findgen(20)
        lev             =       -10. + 2.*findgen(20)
        ;lev             =       -10. + 1.*findgen(40)
	;lev		=	-10. + 0.5*findgen(50)

;;;;;*************************************** AXIS
        ;isotropic       =       'true'		;; DEF: 'true'
        ;intervalx       =       0.5
        intervalx       =        2.0
        ;intervaly       =      intervalx
        intervaly       =       1.0

;;;;;*************************************** MAP LIMITS
        windowx         =       [-146.,-126.] 
        windowy         =       [11.,27.] 

;;;;;***************************************;;;;;
;;;;;***************************************;;;;;
;;;;;***************************************;;;;;
;;;;;***************************************;;;;;


;;;;;*************************************** SETTING TIME (do not modify)
	utc_to_lt       = mean(windowx) / 15.    ;; cf. longitude -- LT = UTC + utc_to_lt
	use_utc 	= use_lt - utc_to_lt
	zentime   	= floor(((24 + use_utc - ini_utc) MOD 24)/freq)        ;; TRUE IDL SUBSCRIPT... ajouter un modulo...
	if (use_lt ne 99) then ntime = zentime else ntime = 99    ;; ou commenter pour avoir tous les pas de temps


;;;;;***************************************;;;;;
;;;;;***************************************;;;;;
;;;;;***************************************;;;;;
;;;;;***************************************;;;;;


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

