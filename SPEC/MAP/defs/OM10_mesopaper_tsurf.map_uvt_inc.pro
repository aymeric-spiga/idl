;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FOLDER 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FOLDER
	folder 	= '/d5/aslmd/LMD_MM_MARS_SIMUS/OM/' 
        coord2d = 'false'	;; for non-regular projections  
filename 	= folder + 'OM10_mesopaper/wrfout_d01_2024-06-43_18:00:00_zabg' 
save_ps 	= 'OM10_mesopaper_tsurf'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; TIME
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; TIME
ini_utc		= 18	;; cf. name of file
freq		= 2  	;; cf. 1 output per ... hour
utc_to_lt       = -9    ;; cf. longitude -- LT = UTC + utc_to_lt 
	use_lt	= 2	;; cf. what user wants
	use_utc	= use_lt - utc_to_lt
        ntime	= ((24 + use_utc - ini_utc) MOD 24)/freq	;; TRUE IDL SUBSCRIPT ... ajouter un modulo... 
        print, 'CHECK ... ', filename, ntime
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; LVL & FLD
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; LVL & FLD
nlevel          = 0
field1   	= 'TSURF'  	;; comment to trace horizontal velocity
no3d 		= 'true'
overvector_x	= 0		;; comment out to get rid of vectors
overvector_y	= 0		;; comment out to get rid of vectors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; TWEAK VAR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; TWEAK VAR
	;a0 = 65.052165 & a1 = 3.1228993 & a2 = 0.0053787417 ;; Ls ~ 173
	;;;a0 = 64.039300 & a1 = 3.1378104 & a2 = 0.0055369148 ;; Ls ~ 120
	;what_I_plot = - ( alog ( a1 - what_I_plot / a0 ) ) /  a2
	;print, max(what_I_plot), min(what_I_plot)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 


;;;;;*************************************** PLOT TITLES
        ;title_user 	= 	'Surface temperature (K) and Winds 10m ABG (m s!U-1!N)' 
	;title_user    	= 	'Wind-induced apparent thermal inertia (J m!U-2!N s!U-0.5!N K!U-1!N)'
	;title_user 	= 	'Winds 10m ABG (m s!U-1!N)'
        title_user      =       'Surface temperature (K)'
        ;subtitle_user   =       'LT = 03:00am  /  Ls = 173!Uo!N  /  dx = 10km [single]  /  Uniform TI = 85 J m!U-2!N s!U-0.5!N K!U-1!N'
        subtitle_user   =       'LT = 02:00am  /  Ls = 173!Uo!N  /  dx = 10km [single]  /  Uniform TI = 85 J m!U-2!N s!U-0.5!N K!U-1!N'
        title_axis      =       ['East longitude','North latitude']

;;;;;*************************************** COLOR TABLES
        flag_cb         =       'true'
	poscb		=	0.96
        format          =       '(I0)'
        ;format          =       '(F4.1)' 
        ;colors          =       128
        colors          =       32
        pal             =       22              ;; GOOD: 4, 18, 22, 16, 37, 33, 39, 6, 11
        ;pal             =       4

;;;;;*************************************** FILL LIMITS
        minfield_init   =       150.
        maxfield_init   =       180.
	;minfield_init	= 	50.
	;maxfield_init	= 	250.

	;;;;;*************************************** LIMIT TRICKS
        ;;;;;********************* must always follow FILL LIMITS
        	missing_value=1.e30
        	lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
        	lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
        	;lim_blank = 2. & w=where(abs(what_I_plot) le lim_blank) & if (w[0] ne -1) then what_I_plot[w]=missing_value

;;;;;*************************************** WINDS
        ;windex          =       10.		;; DEF: 20.
        windex          =       30.             ;; DEF: 20.
        ;stride          =       1. 		;; DEF: 5.
        ;stride          =       2.              ;; DEF: 5.
        stride          =       3.              ;; DEF: 5.

;;;;;*************************************** CONTOUR
	overcontour	=	overcontour/1000.
        ;lev             =       50. + 50.*findgen(20)
        ;lev             =       -10. + 0.2*findgen(20)
        lev             =       -10. + 2.*findgen(20)
        ;lev             =       -10. + 1.*findgen(40)

;;;;;*************************************** AXIS
        isotropic       =       'true'
        ;intervalx       =       0.5
        intervalx       =       2.0
        ;intervaly       =	intervalx
	intervaly	=	1.0

;;;;;*************************************** MAP LIMITS
        ;windowx         =       [-144.,-126.] 
        ;windowy         =       [10.,26.] 
	windowx 	= 	[-146.,-126.]        
	;windowy 	 =	[10.,28.]
	;windowx 	 = 	[-146.,-128.]
        ;windowy         =      [12.,26.]
	windowy 	= 	[11.,27.]

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

