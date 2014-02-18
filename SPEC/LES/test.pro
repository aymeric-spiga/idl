pro test


getcdf, $
;        file='./wrfout_d01_9999-01-01_05:08:21', $
        file='./wrfout_d01_9999-01-01_06:10:00', $
;        file='./wrfout_d01_9999-01-01_07:11:40', $
        charvar='W', $
        invar=w
w1 = max ( w, DIMENSION=3 )
wm1 = min ( TEMPORARY(w), DIMENSION=3 )

;getcdf, $
;        file='./wrfout_d01_9999-01-01_06:10:00', $
;        charvar='W', $
;        invar=w
;w2 = max ( w, DIMENSION=3 )
;wm2 = min ( TEMPORARY(w), DIMENSION=3 )
;
;getcdf, $
;        file='./wrfout_d01_9999-01-01_07:11:40', $
;        charvar='W', $
;        invar=w
;w3 = max ( w, DIMENSION=3 )
;wm3 = min ( TEMPORARY(w), DIMENSION=3 )
;
;getcdf, $
;        file='./wrfout_d01_9999-01-01_04:06:40', $
;        charvar='W', $
;        invar=w
;w4 = max ( w, DIMENSION=3 )
;wm4 = min ( TEMPORARY(w), DIMENSION=3 )

getcdf, $
        file='./wrfout_d01_9999-01-01_06:10:00', $
        charvar='PSFC', $
        invar=psfc

gros=50 ;20;5;10
;gros=5 ;; pour le cas hill
var1 = psfc[ where( w1  gt 10.) ]
var2 = psfc
var3 = psfc[ where( wm1  lt -5.)  ]
;var4 = psfc[ where( w1  gt 13.) ] ;; bien mais deja montre avec courbes max
r1 = ( float(histogram(floor(var1*gros))) ) ;/ n_elements(var1) * 100.
r2 = ( float(histogram(floor(var2*gros))) ) ;/ n_elements(var2) * 100. 
r3 = ( float(histogram(floor(var3*gros))) ) ;/ n_elements(var3) * 100.
r1 = 100. * r1 / max(r1)
r2 = 100. * r2 / max(r2)
r3 = 100. * r3 / max(r3)
t1 = min(var1) + findgen(n_elements(r1))/gros
t2 = min(var2) + findgen(n_elements(r2))/gros
t3 = min(var3) + findgen(n_elements(r3))/gros
PS_START, file='stats_w.ps'
  !P.Charsize = 1.2
  !p.charthick = 2.0
  !p.thick = 2.0
  !x.thick = 2.0
  !y.thick = 2.0
!p.psym = 10 
plot, t2, r2, xrange=[min(t2)+(max(t2)-min(t2))/2., max(t2)], xtickinterval=1, xtitle='Surface pressure (Pa)', ytitle='Relative quantity (ratio to max value in %)';, yrange=[0,4]
oplot, t3, r3, psym=2
oplot, t1, r1, psym=5
PS_END, /PNG

getcdf, $
        file='./wrfout_d01_9999-01-01_00:00:00', $
        charvar='PSFC', $
        invar=psfc

xx=findgen(n_elements(psfc(*,0,0)))*50.
yy=findgen(n_elements(psfc(0,*,0)))*50. 

PS_START, file='psfc.ps'
  !P.Charsize = 1.2
  !p.charthick = 2.0
  !p.thick = 2.0
  !x.thick = 2.0
  !y.thick = 2.0
levu=[685.,687.,689.,691.,693.,695.,697.,699.,701.,703.,705.,707.,709.,711.,713.,715.,717.,719.]

psfc(0,0,1)=684.
psfc(0,1,1)=718.

loadct, 4
contour, reform(psfc(*,*,1)), xx, yy, /cell_fill, nlevels=60, xtitle='x dimension (m)', ytitle='y dimension (m)', /isotropic, max_value=718., min_value=684.
contour, reform(psfc(*,*,1)), xx, yy, lev=levu, c_labels=findgen(n_elements(levu))*0+1, /overplot
PS_END, /PNG

stop


;m = max(psfc)
;psfc = -11000. * alog ( psfc / m )

PS_START, file='stats_w.ps'
  !P.Charsize = 1.2
  !p.charthick = 2.0
  !p.thick = 2.0
  !x.thick = 2.0
  !y.thick = 2.0
plot, psfc, w1, psym=3, title='Maximum of vertical wind amplitude along the vertical (between LT 11 and 15)', xtitle='Surface pressure (Pa)', ytitle='Maximum vertical wind (m s!U-1!N)', yrange=[5.,20.], xrange=[685.,720.]
oplot, psfc, w2, psym=3
oplot, psfc, w3, psym=3
oplot, psfc, w4, psym=3
PS_END, /PNG

PS_START, file='stats_w2.ps'
  !P.Charsize = 1.2
  !p.charthick = 2.0
  !p.thick = 2.0
  !x.thick = 2.0
  !y.thick = 2.0
plot, psfc, wm1, psym=3, title='Maximum of vertical wind amplitude along the vertical (between LT 11 and 15)', xtitle='Surface pressure (Pa)', ytitle='Maximum vertical wind (m s!U-1!N)', yrange=[-12.,-4.], xrange=[685.,720.]
oplot, psfc, wm2, psym=3
oplot, psfc, wm3, psym=3
oplot, psfc, wm4, psym=3

PS_END, /PNG



end
