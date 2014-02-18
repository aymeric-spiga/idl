

PS_Start, FILENAME='velmax.ps'
!P.Charsize = 1.2


history_interval_s = 400.
history_interval_s = 100.
smoothampl=3700/history_interval_s
;smoothampl=0.  ;; no smooth
;smoothampl=20.
smoothampl=5.

lev=1
;lev=0 qu'il faut faire ???


@report_tasi.inc




;restore, filename='REAL_50m_145_145_201_12km/case_A/addturb.dat'
restore, filename=saveplot+cases(0)+'/addturb.dat'
restore, filename=saveplot+cases(0)+'/getturb.dat'
print, h(lev)
velmaxye = reform(velmax(lev,*))
velmaxye = SMOOTH(TEMPORARY(velmaxye), [smoothampl], /EDGE_TRUNCATE)
plot, localtime, velmaxye, $
xrange=[8.,17.], xtickinterval=1., ytickinterval=2.0, $
yrange=[0.,14.], ytitle='Horizontal velocity (m s!U-1!N)', xtitle='Local time (h)'

;restore, filename='REAL_50m_145_145_201_12km/case_B/addturb.dat'
restore, filename=saveplot+cases(1)+'/addturb.dat'
restore, filename=saveplot+cases(1)+'/getturb.dat'
print, h(lev)
velmaxye = reform(velmax(lev,*))
velmaxye = SMOOTH(TEMPORARY(velmaxye), [smoothampl], /EDGE_TRUNCATE)
oplot, localtime, velmaxye, linestyle=1

;restore, filename='REAL_50m_145_145_201_12km/case_C/addturb.dat'
restore, filename=saveplot+cases(2)+'/addturb.dat'
restore, filename=saveplot+cases(2)+'/getturb.dat'
print, h(lev)
velmaxye = reform(velmax(lev,*))
velmaxye = SMOOTH(TEMPORARY(velmaxye), [smoothampl], /EDGE_TRUNCATE)
oplot, localtime, velmaxye, linestyle=2

;restore, filename='REAL_50m_145_145_201_12km/case_I/addturb.dat'
restore, filename=saveplot+cases(3)+'/addturb.dat'
restore, filename=saveplot+cases(3)+'/getturb.dat'
print, h(lev)
velmaxye = reform(velmax(lev,*))
velmaxye = SMOOTH(TEMPORARY(velmaxye), [smoothampl], /EDGE_TRUNCATE)
oplot, localtime, velmaxye, linestyle=3

PS_END, /PNG
