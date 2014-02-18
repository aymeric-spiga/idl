;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ntime    = 19  ;; le vrai indice IDL 
nlevel   = 0 ;1 
filename = '/donnees/aslmd/MODELES/LMD_MM_MARS/OUTILS_CONVERSION/wrfout_d03_2024-01-51_06z00z00_zabg' ;;UTC+1h
save_ps  = 'meridiani'
;coord2d  = 'true'
;field1   = 'TSURF'  ;; omettre trace vitesse horizontale
no3d = 'true'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
format        = '(I0)'
colors        = 128
pal           =  22 ;33 ;39 ;6 ;11 ;0 ;11 ;6 ;0 ;6 ;39 ;19 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; PLOT TITLES                                ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
title_user = 'Surface temperature (K) and Winds 10m ABG (m s!U-1!N)'
subtitle_user='LT = 02:00am  /  Ls = 25!Uo!N  /  dx = 2.25km [nest 3]  /  Uniform TI = 65 J m!U-2!N s!U-0.5!N K!U-1!N'        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FILL LIMITS				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
minfield_init= 150.      
maxfield_init= 180.    
;minfield_init= 155.
;maxfield_init= 175.
missing_value=1.e30
;lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max 
;lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; WINDS				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
windex=15.  	;; default: 20.
stride=2.       ;; default: 5.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; CONTOUR				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
overcontour=overcontour/1000.
lev=-10. + 0.2*findgen(100) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; AXIS				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
intervalx=0.5            
intervaly=0.5            
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; MAP LIMITS 				       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
windowx=[14.0,16.0]        
windowy=[10.0,12.0]


;
; OLYMPUS
;
ntime = 8
filename = '/donnees/aslmd/MODELES/LMD_MM_MARS/OUTILS_CONVERSION/wrfout_d01_2024-06-42_18z00z00_zabg' ;;UTC-9h ;;sorties tt 2h
save_ps = 'olympus'
subtitle_user='LT = 03:00am  /  Ls = 173!Uo!N  /  dx = 20km [single]  /  Uniform TI = 85 J m!U-2!N s!U-0.5!N K!U-1!N'
windex = 30.
stride = 2.
lev = -10. + 2.*findgen(20)
intervalx = 2.            
intervaly = 2.            
windowx = [-146.,-126.]        
windowy = [10.,28.]

filename = '/donnees/aslmd/MODELES/LMD_MM_MARS/OUTILS_CONVERSION/wrfout_d01_2024-06-43_18:00:00_zabg' ;;mm chose a 10 km
save_ps = 'olympus10'
subtitle_user='LT = 03:00am  /  Ls = 173!Uo!N  /  dx = 10km [single]  /  Uniform TI = 85 J m!U-2!N s!U-0.5!N K!U-1!N'
stride = 3.

;pal           =  4 ;33 ;4 ;22 ;33 ;39 ;6 ;11 ;0 ;11 ;6 ;0 ;6 ;39 ;19 
;title_user    = 'Wind-induced apparent thermal inertia (J m!U-2!N s!U-0.5!N K!U-1!N)'
;minfield_init= 30. ;50.
;maxfield_init= 230. ;250.
;a0 = 65.052165 & a1 = 3.1228993 & a2 = 0.0053787417 ;; Ls ~ 173
;;a0 = 64.039300 & a1 = 3.1378104 & a2 = 0.0055369148 ;; Ls ~ 120
;what_I_plot = - ( alog ( a1 - what_I_plot / a0 ) ) /  a2
;print, max(what_I_plot), min(what_I_plot)
;overvector_x=0
;overvector_y=0

;minfield_init= 0.
;maxfield_init= 30.


;;
;; % of error due to winds
;;
;pal=4;18;0;22
;minfield_init=-100.
;maxfield_init=150.
;;what_I_plot = 100. * (what_I_plot / 60. - 1.)
;;what_I_plot = 100. * (what_I_plot / 70. - 1.)
;;what_I_plot = 100. * (what_I_plot / 85. - 1.)
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max 
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min 
;;;commenter au dessus

ntime    = 6
;;filename = '/tmp7/aslmd/results60/wrfout_d01_2024-05-03_06:00:00_zabg'
;;filename = '/tmp7/aslmd/results70/wrfout_d01_2024-05-03_06:00:00_zabg'
;;filename = '/tmp7/aslmd/lowres/wrfout_d01_2024-05-01_06:00:00_zabg'
;filename = '/tmp7/aslmd/wrfout_d01_2024-05-03_06:00:00_zabg'
;filename = '/tmp7/aslmd/OM6/wrfout_d01_2024-06-43_06:00:00_zabg'


windowx = [-146.,-128.]
windowy = [12.,26.]

;windowx=[-115.,-111.]       ;; grepMAP
;windowy=[-1.,3.]       ;; grepMAP

title_user = 'Winds 10m ABG (m s!U-1!N)' 
