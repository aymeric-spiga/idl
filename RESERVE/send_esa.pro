restore, filename='yeahyeah'
filename = 'mesoscale_profile'
tab1 = HEIGHTP

;tab2 = reform(WHAT_I_PLOT(*,0))

tab2 = reform(WHAT_I_PLOT(*,1)+WHAT_I_PLOT(*,2)+WHAT_I_PLOT(*,3)+WHAT_I_PLOT(*,4)+WHAT_I_PLOT(*,5))/5.

make_ascii, filename+'_MEAN', tab1, tab2

plot, tab1, tab2

tab2 = reform(WHAT_I_PLOT(*,1))
make_ascii, filename+'_LT7', tab1, tab2
plot, tab1, tab2

tab2 = reform(WHAT_I_PLOT(*,2))
make_ascii, filename+'_LT8', tab1, tab2
plot, tab1, tab2

tab2 = reform(WHAT_I_PLOT(*,3))
make_ascii, filename+'_LT9', tab1, tab2
plot, tab1, tab2

tab2 = reform(WHAT_I_PLOT(*,4))
make_ascii, filename+'_LT10', tab1, tab2
plot, tab1, tab2

tab2 = reform(WHAT_I_PLOT(*,5))
make_ascii, filename+'_LT11', tab1, tab2
plot, tab1, tab2, yrange=[0.,40.]




