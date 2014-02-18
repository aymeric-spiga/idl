restore, filename='case_B/addturb.dat'
caseb=ustarm
restore, filename='case_C/addturb.dat'
casec=ustarm
restore, filename='case_HIGH/addturb.dat'
casez=ustarm
plot, localtime, caseb
oplot, localtime, casec, linestyle=2
oplot, localtime, casez, linestyle=3

