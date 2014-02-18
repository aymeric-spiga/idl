pro similarity

!p.multi = [0,2,1]


lim=0.30 & zelim=(-3.85/alog(lim)+0.07*alog(lim))*exp(-4.61*lim)
up=1. & & zeup= - 1.55 * up + 1.27
slope = (zelim - zeup) / (lim - up) & ord = zeup - slope * up
t1 = findgen(1001.)*lim/1000. & param1 = (-3.85/alog(t1)+0.07*alog(t1))*exp(-4.61*t1)
t2 = lim + findgen(1001.)*(1.0 - lim)/1000. & param2 = slope * t2 + ord
;;print, slope, ' slope'
;;print, ord, ' ord'
;;      -1.51562 slope      1.23562 ord


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
set_title     =  "Similarity functions: vertical eddy heat flux"
set_titlex    =  "Dimensionless vertical heat flux <w'T'>/<w'T'>!Dmax!N"
set_titley    =  'Dimensionless height z/z!Di!N'
set_subtitle  =  ''
set_xrange    =  [-0.5,1.]
set_yrange    =  [0.0,1.0]
set_tickx     =  0.25
set_ticky     =  0.1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;zesim = 1 - 1.2 * a_h
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PS_Start, FILENAME='similarity.ps'
!P.Charsize = 1.2
plot, param1, t1, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle;, /ISOTROPIC
oplot, param2, t2
oplot, [0.,0.], set_yrange, linestyle=1
oplot, 1 - 1.2 * t1, t1, linestyle=2
oplot, 1 - 1.2 * t2, t2, linestyle=2
xyouts, 0.6, 0.47, 'MARS'
xyouts, 0.25, 0.3, 'EARTH'
;PS_End, /PNG


t1 = findgen(2001.)*1./2000.
param1 = 2.05 * t1^(2./3.) * ( 1. - 0.64 * t1 )^2
param2 = 1.80 * t1^(2./3.) * ( 1. - 0.80 * t1 )^2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
set_title     =  "Similarity functions: vert. vel. variance"
set_titlex    =  "Dimensionless vert. vel. variance <w'!U2!N>/w!D*!N!U2!N"
set_titley    =  'Dimensionless height z/z!Di!N'
set_subtitle  =  ''
set_xrange    =  [0.,1.]
set_yrange    =  [0.0,1.0]
set_tickx     =  0.25
set_ticky     =  0.1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PS_Start, FILENAME='similarity.ps'
!P.Charsize = 1.2
plot, param1, t1, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle;, /ISOTROPIC
oplot, param2, t1, linestyle=2
xyouts, 0.2, 0.3, 'EARTH'
xyouts, 0.67, 0.4, 'MARS'
PS_End, /PNG



end
