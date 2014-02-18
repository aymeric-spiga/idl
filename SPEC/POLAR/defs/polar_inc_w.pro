;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ntime    = 2 
nlevel   = 1 
latmin = -90 & latmax = 20 & lonmin = -180 & lonmax = 180 & hem = -1.
latmin = -40 & latmax = -10 & lonmin = -45 & lonmax = -15 & hem = -1.
filename = '/tmp1/aslmd/wrfout_d01_2024-06-17_06:00:00' & save_ps  = 'wpolarz1'
;filename = '/tmp1/aslmd/wrfout_d02_2024-06-17_06:00:00' & save_ps  = 'wpolarz2'
;filename = '/tmp1/aslmd/wrfout_d03_2024-06-17_06:00:00' & save_ps  = 'wpolarz3'
;filename = '/tmp1/aslmd/wrfout_d04_2024-06-17_06:00:00' & save_ps  = 'wpolarz4'
latmin = -30 & latmax = -24 & lonmin = -38 & lonmax = -32
filename = '/tmp1/aslmd/wrfout_d04_2024-06-17_06:00:00' & save_ps  = 'wpolarzz4'
cfield1  = 'WAVE'
cfield2  = 'HGT'
;cfield3  = 'HGT_M'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
format        = '(F4.1)'
colors        = 128
minfield_init = -8. 
maxfield_init = +2.
pal           =  22 ;33 ;39 ;6 ;11 ;0 ;11 ;6 ;0 ;6 ;39 ;19 
title_user    = ''
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
missing_value=1.e30
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
;lim_blank = 0.5 & w=where(abs(what_I_plot) le lim_blank) & if (w[0] ne -1) then what_I_plot[w]=missing_value
overcontour = overcontour / 1000.
lev = -10. + 1.*findgen(40)
;flag_cb = 'false'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;
;; T. Appere paper
;;
latmin = 76 & latmax = 90. & lonmin = -180 & lonmax = 180 & hem = 1.
ntime    = 0  ;; le vrai indice IDL 
nlevel   = 0
filename = '/tmp15/aslmd/polar_61pts/wrfout_d01_2024-03-04_06:00:00_zabg'
save_ps  = 'appere'
cfield1  = 'USTM'
minfield_init = 0.1
maxfield_init = 0.6
pal           =  22
format        = '(F5.2)'
colors        = 32
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
lev = -10. + 0.2*findgen(50)
;lev = -10. + 0.5*findgen(20)
overcontour = - overcontour
lev = 0.5*findgen(20)
lev = 0.2*findgen(50)




