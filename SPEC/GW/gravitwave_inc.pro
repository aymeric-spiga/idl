;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;; W 
colors        = 64
minfield_init = -10.
maxfield_init = 10. 
ndiv	      = 10
format        = '(I0)'
pal           = 17 ;4 ;16 ;18 ;4
title_user    = 'Vertical velocity (m s!U-1!N)'
lev=80.+10.*findgen(30)
;;;;;;;; tk
;minfield_init = 110.
;maxfield_init = 230.
;ndiv	      = 12
;format        = '(I0)'
;pal	      = 33
;title_user    = 'Temperature (K)'
;lev = [-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10.]
;lev=80.+10.*findgen(30)
;;;;;;;; tk - tsat
;minfield_init =  -4. ;-8.
;maxfield_init =  4. ;+8.
;ndiv	      =  8	
;format        = '(F4.1)'
;pal           =  0
;title_user    = 'T!Datm!N - T!Dsat!N (K)'
;lev = [-6.,-5.,-4.,-3.,-2.,-1.,0.]
;;;;;;;;; hr_nlte
;	minfield_init =  -6.
;	maxfield_init =  6.
;	ndiv          =  6
;	format        = '(F4.1)'
;	pal           =  33 
;	title_user    = 'NLTE heating rate (K / hour)'
lev=[-20.,-18.,-16.,-14.,-12.,-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.]
title_user    = 'NLTE cooling rate (K hour!U-1!N)'
w=where(abs(what_I_plot) le 1e5) & if (w[0] ne -1) then what_I_plot[w] = - what_I_plot[w]
minfield_init = 0.
maxfield_init = 20.
ndiv          = 10
pal	      = 1 ;4 ;22
lev=80.+10.*findgen(30)
format	      = '(I0)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
missing_value=1.e30
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
flag_cb = 'true'
windowx = [5.,55.]
windowy = [0.,100.]
windowy = [75.,85.]
windowx = [5., 35.] 
intervalx = 5.
intervaly = 5.
windowy = [0.,130.]
windowx = [10., 110.] 
;windowy = [50.,130.]
intervalx = 10.
intervaly = 10.
title_axis = ['horizontal coordinate (x5km)','altitude (km)']
poscb=0.55
;isotropic='false'
;lev = 10.*findgen(100)
;lev = [-6.,-4.,-2.,2.,4.,6.]
;;lev = [-8.,-5.,-1.,1.,5.,8.]
;;lev = [-6.,-5.,-4.,-3.,-2.,-1.,0.]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
