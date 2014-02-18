;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
format        = '(F4.1)'
colors        = 128
;colors	      = 16	
minfield_init = -8.
maxfield_init =  12. 
pal = 33
	;
	;minfield_init = -2.
	;maxfield_init =  2.
	;
;minfield_init = 0.;0.1
;maxfield_init =  1.;1.1
;;pal           =  22;33  
;minfield_init =  0.
;maxfield_init =  10. 
;pal           =  22  
;title_user    = ''
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
missing_value=1.e30
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
flag_cb = 'true'
windowx = [0.,20.] 
windowy = [0.,20.] 
intervalx = 2. 
intervaly = 2. 
	windowx = [0.,9.]
	windowy = [0.,9.]
	intervalx = 1.
	intervaly = 1.
;;
        windowx = [0.,8.]
        windowy = [6.,14.]
        intervalx = 1.0
        intervaly = 1.0
;;
title_axis = ['x-axis (km)','y-axis (km)']
poscb = 0.70
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
windex = 20 
stride = 1
