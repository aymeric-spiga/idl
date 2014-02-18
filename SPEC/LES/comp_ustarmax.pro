
PS_Start, FILENAME='ustarmax.ps'
!P.Charsize = 1.2



history_interval_s = 400.
history_interval_s = 100.
smoothampl=3700/history_interval_s
smoothampl=0.  ;; no smooth
;smoothampl=20.

lev=1


;restore, filename='REAL_50m_145_145_201_12km/case_A/addturb.dat'
;ustarmaxye = SMOOTH(TEMPORARY(ustarmax), [smoothampl], /EDGE_TRUNCATE)
;ustarmye = SMOOTH(TEMPORARY(ustarm), [smoothampl], /EDGE_TRUNCATE)
;plot, localtime, ustarmaxye, yrange=[0.0,1.0], xtickinterval=1., ytickinterval=0.1, ytitle='m/s', xtitle='Local time (h)'
;oplot, localtime, ustarmye

restore, filename='REAL_50m_145_145_201_12km/case_B/addturb.dat'
ustarmaxye = SMOOTH(TEMPORARY(ustarmax), [smoothampl], /EDGE_TRUNCATE)
ustarmye = SMOOTH(TEMPORARY(ustarm), [smoothampl], /EDGE_TRUNCATE)
plot, localtime, ustarmaxye, yrange=[0.0,1.2], xtickinterval=1., ytickinterval=0.1, ytitle='m/s', xtitle='Local time (h)'
oplot, localtime, ustarmye

restore, filename='REAL_50m_145_145_201_12km/case_C/addturb.dat'
ustarmaxye = SMOOTH(TEMPORARY(ustarmax), [smoothampl], /EDGE_TRUNCATE)
ustarmye = SMOOTH(TEMPORARY(ustarm), [smoothampl], /EDGE_TRUNCATE)
oplot, localtime, ustarmaxye, linestyle=2
oplot, localtime, ustarmye, linestyle=2

;restore, filename='REAL_50m_145_145_201_12km/case_I/addturb.dat'
;ustarmaxye = SMOOTH(TEMPORARY(ustarmax), [smoothampl], /EDGE_TRUNCATE)
;oplot, localtime, ustarmaxye, linestyle=3

;restore, filename='REAL_50m_145_145_201_12km/case_HIGH/addturb.dat'
;ustarmaxye = SMOOTH(TEMPORARY(ustarmax), [smoothampl], /EDGE_TRUNCATE)
;oplot, localtime, ustarmaxye, linestyle=4




PS_END, /PNG
