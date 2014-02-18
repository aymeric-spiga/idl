pro ye2, r, t

;;paper
min_lt = 12.
max_lt = 14.

;min_lt = 13.
;max_lt = 16.

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
r = TOTAL(a_vel_bot[*,wspe],2)
t = TOTAL(a_h[*,wspe],2)
rtot = [r]
ttot = [t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_B'+'/compturb.dat'
restore, filename=wherefiles+'case_B'+'/getturb.dat'
r = TOTAL(a_vel_bot[*,wspe],2)
t = TOTAL(a_h[*,wspe],2) 
rtot = [rtot,r]
ttot = [ttot,t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_C'+'/compturb.dat'
restore, filename=wherefiles+'case_C'+'/getturb.dat'
r = TOTAL(a_vel_bot[*,wspe],2)
t = TOTAL(a_h[*,wspe],2)
rtot = [rtot,r]
ttot = [ttot,t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_I'+'/compturb.dat'
restore, filename=wherefiles+'case_I'+'/getturb.dat'
r = TOTAL(a_vel_bot[*,wspe],2)
t = TOTAL(a_h[*,wspe],2)
rtot = [rtot,r]
ttot = [ttot,t]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
restore, filename=wherefiles+'case_HIGH'+'/compturb.dat'
restore, filename=wherefiles+'case_HIGH'+'/getturb.dat'
r = TOTAL(a_vel_bot[*,wspe],2)
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

;t1 = findgen(2001.)*1./2000.
;param1 = 2.05 * t1^(2./3.) * ( 1. - 0.64 * t1 )^2
;r = param1
;t = t1

end
