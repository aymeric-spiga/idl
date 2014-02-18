pro ye, r, t

;; paper
min_lt = 13.
max_lt = 16.

;min_lt = 12.
;max_lt = 16.

;min_lt = 13.
;max_lt = 15.


wherefiles='./REAL_50m_145_145_201_12km/'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_A'+'/compturb.dat'
restore, filename=wherefiles+'case_A'+'/getturb.dat'
wspe = where( (localtime ge min_lt) and (localtime le max_lt) )
nspe = n_elements(wspe)
r = TOTAL(a_wt_bot[*,wspe],2) 
t = TOTAL(a_h[*,wspe],2) 
rtot = [r]
ttot = [t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_B'+'/compturb.dat'
restore, filename=wherefiles+'case_B'+'/getturb.dat'
r = TOTAL(a_wt_bot[*,wspe],2)
t = TOTAL(a_h[*,wspe],2) 
rtot = [rtot,r]
ttot = [ttot,t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_C'+'/compturb.dat'
restore, filename=wherefiles+'case_C'+'/getturb.dat'
r = TOTAL(a_wt_bot[*,wspe],2)
t = TOTAL(a_h[*,wspe],2)
rtot = [rtot,r]
ttot = [ttot,t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_I'+'/compturb.dat'
restore, filename=wherefiles+'case_I'+'/getturb.dat'
r = TOTAL(a_wt_bot[*,wspe],2)
t = TOTAL(a_h[*,wspe],2)
rtot = [rtot,r]
ttot = [ttot,t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_HIGH'+'/compturb.dat'
restore, filename=wherefiles+'case_HIGH'+'/getturb.dat'
r = TOTAL(a_wt_bot[*,wspe],2)
t = TOTAL(a_h[*,wspe],2)
rtot = [rtot,r]
ttot = [ttot,t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



r = r / float(nspe)
t = t / float(nspe)

sssss = sort(t)
r = r[sssss]
t = t[sssss]

;;
;; OK
;;
;lim=0.30 & zelim=(-3.85/alog(lim)+0.07*alog(lim))*exp(-4.61*lim)
;up=1. & & zeup= - 1.55 * up + 1.27
;slope = (zelim - zeup) / (lim - up) & ord = zeup - slope * up
;t1 = findgen(1001.)*lim/1000. & param1 = (-3.85/alog(t1)+0.07*alog(t1))*exp(-4.61*t1)
;t2 = lim + findgen(1001.)*(1.0 - lim)/1000. & param2 = slope * t2 + ord
;
;r = [param1, param2]
;t = [t1, t2]


end
