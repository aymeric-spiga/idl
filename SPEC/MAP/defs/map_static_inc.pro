;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filename = '/d5/aslmd/LMD_MM_MARS_SIMUS/OM/OM6_TI85/geo_em.d01.nc'
save_ps  = 'olympus_slope'
	field1  = 'THERMAL_INERTIA'  ;; omettre trace pentes
	save_ps  = 'olympus_ti'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
format        = '(I0)'
colors        = 32 ;64 ;128 ;32
pal           = 6 ;39 ;4 ;33 ;4 ;33 ;4 ;22 ;33 ;39 ;6 ;11 ;0 ;11 ;6 ;0 ;6 ;39 ;19 
	pal   = 4	
title_user    = 'Slope angle (deg)'
	title_user    =  'Apparent nighttime thermal inertia (J m!U-2!N s!U-0.5!N K!U-1!N)'
title_axis    = ['East longitude','North latitude']
poscb=0.80
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;what_I_plot = what_I_plot / 1000.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FILL LIMITS				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
minfield_init= 0.      
maxfield_init= 25.    
	minfield_init= 50.
	maxfield_init= 250.
missing_value=1.e30
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max 
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; WINDS				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
windex=20.  	;; default: 20.
stride=3.       ;; default: 5.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; CONTOUR				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
overcontour=overcontour/1000.
lev=-10. + 2.*findgen(20) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; AXIS				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
intervalx=2.0            
intervaly=1.0            
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; MAP LIMITS 				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
windowx = [-146.,-126.]
windowy = [11.,27.]
