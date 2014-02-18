pro ti_vs_tsurf

getcdf, $
        file='/d6/vblmd/MERIDIANI_EXOMARS/saves_simu_LS_244_tau05/wrfout_d01_2024-09-08_00:00:00', $
        charvar='TSURF', $
        invar=tsurf
tsurf = max ( tsurf, DIMENSION=3 )

getcdf, $
        file='./geo_em.d01.nc', $
        charvar='THERMAL_INERTIA', $
        invar=ti

getcdf, $
        file='./geo_em.d01.nc', $
        charvar='ALBEDO_GCM', $
        invar=alb

PS_START, file='ti.ps'
  !P.Charsize = 1.2
  !p.charthick = 2.0
  !p.thick = 2.0
  !x.thick = 2.0
  !y.thick = 2.0
plot, ti, tsurf, psym=3, title='Daily maximum of surface temperature at ExoMars landing season and region', xtitle='Thermal Inertia (tiu)', ytitle='Surface temperature (K)', xrange=[0,500], yrange=[270,320]
oplot, ti[where(tsurf eq max(tsurf))], tsurf[where(tsurf eq max(tsurf))], psym=4
xyouts, 380, 310, 'TI = '+string(ti[where(tsurf eq max(tsurf))],'(I0)')+' tiu'
PS_END, /PNG

PS_START, file='alb.ps'
  !P.Charsize = 1.2
  !p.charthick = 2.0
  !p.thick = 2.0
  !x.thick = 2.0
  !y.thick = 2.0
plot, alb, tsurf, psym=3, title='Daily maximum of surface temperature at ExoMars landing season and region', xtitle='Albedo', ytitle='Surface temperature (K)', xrange=[0.1,0.3], yrange=[270,320]
oplot, alb[where(tsurf eq max(tsurf))], tsurf[where(tsurf eq max(tsurf))], psym=4
xyouts, 0.25, 310, 'ALB = '+string(alb[where(tsurf eq max(tsurf))],'(F4.2)')
PS_END, /PNG


end
