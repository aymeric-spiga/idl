

SPAWN, '\rm param_plot.idl'
SPAWN, "echo 'intervaly=5.' > param_plot.idl"
SPAWN, "echo 'intervalx=5.' >> param_plot.idl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
file='gcm' & header='' & nlines_header=4 & ncol = 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nlines = FILE_LINES(file)-nlines_header & data=FLTARR(ncol,nlines)
OPENR, lun, file, /GET_LUN & READF, lun, header & READF, lun, header & READF, lun, header & READF, lun, header 
READF, lun, data & CLOSE, lun & FREE_LUN, lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
altitude        = reform(data(0,*))
temp            = reform(data(1,*)) 
data = 0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
file='meso' & header='' & nlines_header=4 & ncol = 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nlines = FILE_LINES(file)-nlines_header & data=FLTARR(ncol,nlines)
OPENR, lun, file, /GET_LUN & READF, lun, header & READF, lun, header & READF, lun, header & READF, lun, header
READF, lun, data & CLOSE, lun & FREE_LUN, lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
altitude2       = reform(data(0,*))         
temp2           = reform(data(1,*)) 
data = 0.


;
; extremes
;
temp = temp - 12.
temp2 = temp2 - 12.
print, 'ATTENTION CACA CACA CACA'

what_I_plot = temp ;- temp2
column = altitude
alt = [20.,60.]
minfield_init = 100. ;-20 ;100.
maxfield_init = 200. ;20. ;200.
overplot = temp2
overplot_column = altitude2
discrete = 0
title_user = 'GCM vs. mesoscale simulations'
title_axis = ['Temperature (K)','Altitude above MOLA (km)']
mention = ''

profile, $
        what_I_plot, $                          ; 1D vertical profile
        column, $                               ; altitudes     
        alt=alt, $                              ; altitude range [altmin, altmax]
        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
        inprofile=overplot, $                   ; another vertical profile to overplot
        incolumn=overplot_column, $             ; altitudes of the other vertical profile (in case /= column)
        discrete=discrete, $                    ; show the profile points (= type of points in !psym)
        title_plot=title_user, $                ; title of the plot ('Profile' is default)
        title_axis=title_axis, $                ; title of the [x,y] axis (['Field','Altitude'] is default)
        mention=mention                         ; add text precision within the plot window (default is nothing or '')



;
; condensation
;
latcond=5.9e5
tcond1mb=136.27
r=192.
bcond=1./tcond1mb
acond=r/latcond
press=610.*exp(-altitude/10.)
overplot=1./(bcond-acond*alog(.0095*press))
overplot_column=altitude
oplot, overplot, overplot_column, psym=5

