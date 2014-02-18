pro gravitwaveprof
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
; faire avant api avec kappa 3.9
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	path='../TESTGW/'
	path='/tmp7/aslmd/'
path='/home/aslmd/GRAVITWAVE/'
	experiment='mcd_lt16_rcpcst'
	;experiment='mcd_lt16'
	;experiment='short_mcd_lt16_rcpcst'
	;experiment='short_mcd_lt16'
	;experiment='short_mcd_lt16_calcprho'
	;experiment='mcd_lt16_calcprho'
	;experiment='short_mcd_lt16_calcprho_tot'
	;experiment='colder_gcm_add_short'
	;experiment='colder_gcm_add'
	experiment=''
	;experiment='wind15'
	;experiment='wind20c'
	;experiment='colder_gcm_add_3.9_fine_rad_mountain'
	experiment='test'
experiment='GW_MARS_highwind_3D.157077'
experiment='GW_MARS_highwind_3D_lower_opacity.161436'
experiment='GW_MARS_highwind_3D_loweropacity_lotspoints.161468'
experiment='GW_MARS_highwind_3D_loweropacity_morepoints.161542'
;experiment='GW_MARS_highwind_3D_loweropacity_morepoints_LT15.161544'
;experiment='GW_MARS_highwind_3D_loweropacity_morepoints_LTcst.168440'
;experiment='GW_MARS_highwind_3D_loweropacity_morepoints_widehill.168468'
	file=path+experiment+'_wrfout_d01_9999-09-09_09:00:00_z'
file=path+experiment+'/wrfout_d01_9999-09-09_09:00:00_z'
;file=path+experiment+'/wrfout_d01_9999-09-09_10:00:00_z'
;file=path+experiment+'/wrfout_d01_9999-09-09_11:00:00_z'
charvar='W' & charvarc='W'
charvar='tk' & charvarc='W'
charvar='tk' & charvarc='W' & cond=1
charvar='tpot' & charvarc='tpot'
;charvar='PTOT' & charvarc='PTOT'
charvar='tk' & charvarc='tk'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nts=1  & nte=15   & nxs=0   & nxe=59
nts=3  & nte=nts  & nxs=999 & nxe=999 ;; moyenne ;; attention a la montagne !!!
nts=4  & nte=nts  & nxs=31  & nxe=nxs ;; AU CENTRE ! sinon choisir le point en traçant de 0a 59
nts=1  & nte=nts  & nxs=0   & nxe=59
nts=4  & nte=nts  & nxs=0   & nxe=59
nts=1  & nte=15   & nxs=31  & nxe=31
nts=1  & nte=15   & nxs=45  & nxe=45
nts=9  & nte=9    & nxs=0   & nxe=59
nts=1  & nte=13   & nxs=59  & nxe=59
nts=8  & nte=8    & nxs=0   & nxe=121
nts=1  & nte=9    & nxs=59  & nxe=59
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nts=1  & nte=10  & nxs=50   & nxe=50
;nts=1  & nte=1  & nxs=6   & nxe=6
nts=1  & nte=100  & nxs=50   & nxe=50
nts=1  & nte=100  & nxs=60   & nxe=60
;nts=1  & nte=100  & nxs=20   & nxe=20
nts=29 & nte=29   & nxs=60   & nxe=60
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
what_I_plot=0. & overcontour=0.
;SPAWN, '\rm param_plot.idl ; cp gravitwave_inc2.pro param_plot.idl'
;
;
;
getcdf, $
        file=file, $
        charvar=charvar, $
        invar=invar
;getcdf, $
;        file=file, $
;        charvar=charvarc, $
;        invar=invarc
getcdf, $
        file=file, $
        charvar='vert', $
        invar=vert
getcdf, $
        file=file, $
        charvar='PTOT', $
        invar=columnp
;;
;; 
;;
	s = size(invar) & middle = floor(s[2]/2) & print, 'PLOT at y subs ', middle
        nte = min(s[4]-1,nte) & print, 'nte ', nte
        mean_invar = TOTAL(invar(*,*,*,nts:nte),4) / float(nte-nts+1)
	;mean_columnp = TOTAL(columnp(*,*,*,nts:nte),4) / float(nte-nts+1)
        anomal_invar = invar & for i=0,s[4]-1 do anomal_invar(*,*,*,i) = anomal_invar(*,*,*,i) - mean_invar
;;
;;
;;
;	goto, perturb
for nt=nts,nte do begin
for nx=nxs,nxe do begin
;
zefile=experiment+'_'+charvar+'_'+string(nt+100,'(I0)')+'_'+string(nx+100,'(I0)')
;
;
;
  PS_Start, filename=zefile+'.ps'
  print, zefile+'.ps'
  !P.Charsize = 1.2
  !p.charthick = 2.0
  !p.thick = 2.0
  !x.thick = 2.0
  !y.thick = 2.0
;
;
;
;;what_I_plot = reform(invar(*,1,*,nt))
;what_I_plot = reform(invar(*,59,*,nt))
what_I_plot = reform(invar(*,middle,*,nt))

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file=path+experiment+'/input_sounding' & header='' & nlines_header=1 & ncol = 5
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   nlines = FILE_LINES(file)-nlines_header & data=FLTARR(ncol,nlines)
   OPENR, lun, file, /GET_LUN & READF, lun, header & READF, lun, data & CLOSE, lun & FREE_LUN, lun
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   mcd_column      = reform(data(0,*)) 
   calc_tpot       = reform(data(1,*))

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file=path+experiment+'/input_therm' & ncol = 5
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   nlines = FILE_LINES(file) & data=FLTARR(ncol,nlines)
   OPENR, lun, file, /GET_LUN & READF, lun, data & CLOSE, lun & FREE_LUN, lun
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   pressure        = reform(data(2,*))
   mcd_tpot        = reform(data(4,*))   

   if (nx lt 999) then begin
     what_I_plot = reform(what_I_plot(nx,*))
   endif else begin
     what_I_plot = TOTAL(what_I_plot,1) / n_elements(reform(what_I_plot(*,0)))
   endelse

     ;;;;;;;;;;;
     ;;;;;;;;;;;
     ;;;;;;;;;;;
     column = vert/1000.
     overplot_column = mcd_column/1000.
     alt = [20.,150.]
     alt = [40.,100.]
     alt = [40.,80.]
     alt = [0.,120.]
     alt = [10.,110.]
     ;;;;;;;;;;;
     ;;;;;;;;;;;
     ;mention = ''
     SPAWN, '\rm param_plot.idl'
     minfield_init = 85.  ;100.
     maxfield_init = 185. ;250.
     minfield_init = 90.
     maxfield_init = 230.
     title_axis = ['Temperature (K)','Altitude (km)']
     overplot = mcd_tpot
     ;;;;;;;;;;;
     SPAWN, "echo 'intervaly=10.' >> param_plot.idl"
     SPAWN, "echo 'intervalx=5.' >> param_plot.idl"
     SPAWN, "echo 'intervalx=10.' >> param_plot.idl"
     ;;;;;;;;;;;

	;overplot = reform(mean_invar(nx,middle,*))
	;overplot_column = column

;;;;;;;;;;;;;;;;;;;;;;;;;;;; PLOT PRESSION mettre profile en log
;column = reform(columnp(nx,middle,*,nt))
;overplot_column = pressure
;	;;overplot = mean_invar
;	;;overplot_column = mean_columnp
;alt = [1.e2,1.e-3]
;title_axis = ['Temperature (K)','Pressure (Pa)']
;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PLOT PRESSION

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; calculer la CAPE ? Tsat - Tenv / Tenv


			;what_I_plot = pressure
			;	what_I_plot = calc_tpot
			;overplot = alog(overplot) / alog(10.)
			;overplot = overplot[where(overplot lt 35.)] ;; evite les valeurs impossibles
			;what_I_plot = alog(what_I_plot) / alog(10.)
			;minfield_init = -5.
			;maxfield_init = 3.
			;	;;tpot
			;	minfield_init = 2.
			;	maxfield_init = 4.
			;what_I_plot = 100. * abs(calc_tpot - overplot) / overplot
			;minfield_init = 0.000001
			;maxfield_init = 100.

   profile, $
        what_I_plot, $                          ; 1D vertical profile
        column, $                               ; altitudes     
        alt=alt, $                              ; altitude range [altmin, altmax]
        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
        inprofile=overplot, $                   ; another vertical profile to overplot
        incolumn=overplot_column, $             ; altitudes of the other vertical profile (in case /= column)
;        discrete=discrete, $                    ; show the profile points (= type of points in !psym)
;        title_plot=title_user, $                ; title of the plot ('Profile' is default)
        title_axis=title_axis, $                ; title of the [x,y] axis (['Field','Altitude'] is default)
        mention=mention                         ; add text precision within the plot window (default is nothing or '')

	;;; sponge layer
	;w=where(FINITE(overplot) ne 0)
	;oplot, [minfield_init,maxfield_init], [max(overplot_column[w])-50.,max(overplot_column[w])-50.], linestyle=1

	;column = columnp
	;@tempcond.inc
	;yeye = reform(overplot(nx,1,*,nt)) 
	;;w = where(yeye le 0.) & yeye[w] = !VALUES.F_NAN
	;;oplot, yeye, vert/1000., linestyle=1


loadct, 4
for i=16,37 do begin  ;10-37
for j=0,40 do begin
dec = j - 20  
overoverplot = reform(invar(nx+dec,middle+dec,*,i))
oplot, overoverplot, column, color=220  ;psym=3
endfor
endfor

!p.thick = 4.0
oplot, overplot, overplot_column, linestyle=2
!p.thick = 2.0

!p.thick = 4.0
oplot, reform(invar(nx,middle,*,28)), column
!p.thick = 2.0

@tempcond.inc
overplot = reform(overplot(nx,middle,*,nt))
w = where(overplot le 0.) & overplot[w] = !VALUES.F_NAN
oplot, overplot, overplot_column, linestyle=1

	;calc_tk = calc_tpot * (pressure/610.)^(1/4.4)
	;oplot, calc_tk, column, linestyle=3

PS_End, /PNG
endfor
endfor

perturb:
PS_Start, filename=experiment+'_'+charvar+'_perturb_'+string(nxs+100,'(I0)')+'.ps'

     what_I_plot_loop = reform(anomal_invar(nxs,middle,*,*))

     what_I_plot = reform(what_I_plot_loop(*,nts))
     column = vert/1000.
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     alt = [0.,120.]
     minfield_init = -30. 
     maxfield_init = 30. 
     title_axis = ['Temperature anomaly (K)','Altitude (km)']
		;; peut pas changer	
		;SPAWN, "echo 'intervaly=10.' >> param_plot.idl"
     		;SPAWN, "echo 'intervalx=5.' >> param_plot.idl"
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     profile, $
        what_I_plot, $                          ; 1D vertical profile
        column, $                               ; altitudes     
        alt=alt, $                              ; altitude range [altmin, altmax]
        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
        inprofile=overplot, $                   ; another vertical profile to overplot
        incolumn=overplot_column, $             ; altitudes of the other vertical profile (in case /= column)
;        discrete=discrete, $                    ; show the profile points (= type of points in !psym)
;        title_plot=title_user, $                ; title of the plot ('Profile' is default)
        title_axis=title_axis, $                ; title of the [x,y] axis (['Field','Altitude'] is default)
        mention=mention                         ; add text precision within the plot window (default is nothing or '')

     for iii=nts+1, nte do begin
        overplot = reform(what_I_plot_loop(*,iii))
        overplot_column = column
        oplot, overplot, overplot_column
     endfor

	;temppp = reform(mean_invar(nxs,middle,*))
	;hache = 191.*temppp/3.72/1000.
	;ok = column / ( 2. * hache )
	;print, exp(ok)
	;yeah = exp(ok)
	;oplot, yeah, column

PS_End, /PNG

end
