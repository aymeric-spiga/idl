;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ntime    = 0
nlevel   = 0
latmin = -90 & latmax = 20 & lonmin = -180 & lonmax = 180
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d01.nc' & save_ps  = 'polar1'
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d02.nc' & save_ps  = 'polar2'
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d03.nc' & save_ps  = 'polar3'
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d04.nc' & save_ps  = 'polar4'
latmin = -40 & latmax = -10 & lonmin = -45 & lonmax = -15
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d04.nc' & save_ps  = 'polarz4'
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d03.nc' & save_ps  = 'polarz3'
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d02.nc' & save_ps  = 'polarz2'
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d01.nc' & save_ps  = 'polarz1'
latmin = -30 & latmax = -24 & lonmin = -38 & lonmax = -32
filename = '/donnees/aslmd/MODELES//LMD_MM_MARS/mpi_64/WPS/geo_em.d04.nc' & save_ps  = 'polarzzz4'
cfield1  = 'HGT_M'
;cfield1  = 'THERMAL_INERTIA'
cfield2  = 'HGT_M'
cfield3  = 'HGT_M'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
format        = '(I0)'
colors        = 128
minfield_init = -8000 
maxfield_init = +12000
pal           =  33 ;39 ;6 ;11 ;0 ;11 ;6 ;0 ;6 ;39 ;19 
;minfield_init = 150.
;maxfield_init = 550.
;pal           =  3
minfield_init = -2500
maxfield_init = +1500
pal           =  33  
title_user    = ''
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
missing_value=1.e30
lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
;lim_blank = 0.5 & w=where(abs(what_I_plot) le lim_blank) & if (w[0] ne -1) then what_I_plot[w]=missing_value
overcontour = overcontour / 1000.
lev = -10. + 1.*findgen(40)
lev = -10. + 0.5*findgen(80)
flag_cb = 'false'
flag_cb = 'true'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

