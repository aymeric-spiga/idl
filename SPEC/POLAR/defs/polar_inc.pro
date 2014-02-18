;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ntime    = 10
nlevel   = 0
filename = './wrfout_d01_2024-07-03_06:00:00_z'
save_ps  = 'polar'
cfield1  = 'W'
cfield2  = 'HGT'
cfield3 = 'SWDOWNZ'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
format        = '(F4.1)'
colors        = 32 ;128
minfield_init = -2.50
maxfield_init = 2.50
pal           =  33 ;39 ;6 ;11 ;0 ;11 ;6 ;0 ;6 ;39 ;19 
title_user    = ''
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
missing_value=1.e30
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
lim_blank = 0.5 & w=where(abs(what_I_plot) le lim_blank) & if (w[0] ne -1) then what_I_plot[w]=missing_value
overcontour = overcontour / 1000.
lev = -16. + 2.*findgen(40)
;lev = [1, 100, 200, 300, 400, 500]
flag_cb = 'false'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

