;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
format        = '(I0)'
colors	      = 64 ;32	
minfield_init = -8.
maxfield_init = 12. 
pal = 33 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
missing_value=1.e30
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
flag_cb = 'true'
windowx = [0.,9.] 
windowy = [0.,12.] 
intervalx = 1. 
intervaly = 1. 
title_axis = ['Horizontal coordinate (km)','Altitude above surface (km)']
title_user = 'Vertical wind (m s!U-1!N) in Meridiani LES'
subtitle_user = 'DOD = 0.5 / BW = 10 m s!U-1!N'
poscb = 0.55
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
