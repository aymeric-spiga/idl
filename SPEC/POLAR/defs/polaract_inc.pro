;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filename = 'wrfout_d01_2024-07-03_06:00:00_z'
;filename = 'wrfout_d01_2024-07-02_06:00:00_z'
;filename = 'wrfout_d01_2024-07-04_06:00:00_z'
;filename = 'wrfout_d01_2024-07-99_06:00:00_z' ;;1234
filename = 'wrfout_d01_2024-07-234_06:00:00_z' ;;234
;filename = 'wrfout_d01_2024-07-34_06:00:00_z'
;filename = 'wrfout_d01_2024-07-23_06:00:00_z'
save_ps  = 'polar'
cfield1  = 'W'
cfield2  = 'HGT'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
format        = '(F4.1)'
colors        = 32 ;128
minfield_init = -5000.
maxfield_init = 5000.
pal           =  0 ;4;16;22 ;33 ;39 ;11 ;0 39 19 ;NON: 6,11
title_user    = ''
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
missing_value=1.e30
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
;lim_blank = 1.5 & w=where(abs(what_I_plot) le lim_blank) & if (w[0] ne -1) then what_I_plot[w]=missing_value
;overcontour = overcontour / 1000.
;lev = -16. + 2.*findgen(40)
lev = [0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
;lev = [0.5, 0.8,1.2,1.6,2.0] 
;lev = [0.75, 1.,2.,3.]
lev = [0.2, 0.3, 0.5, 1., 2., 3., 5., 10.]
;lev = [0.3, 0.5, 1., 2., 3., 5., 10.]
;lev = [0.25, 0.5, 0.75, 1.]
flag_cb = 'false'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

