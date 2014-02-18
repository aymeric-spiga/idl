pro scorerprof
;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
path='../TESTGW/'
experiment=''
	;;path='../TESTGW_save/gw/TESTGW/'
	;;experiment='mcd_lt16_calcprho'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
what_I_plot=0. & overcontour=0.
zefile=experiment+"_scorer"
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
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file=path+experiment+'/input_sounding' & header='' & nlines_header=1 & ncol = 5
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   nlines = FILE_LINES(file)-nlines_header & data=FLTARR(ncol,nlines)
   OPENR, lun, file, /GET_LUN & READF, lun, header & READF, lun, data & CLOSE, lun & FREE_LUN, lun
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   z      = reform(data(0,*))  
   theta  = reform(data(1,*))
   windu  = reform(data(3,*))   
           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           file=path+experiment+'/input_therm' & ncol = 5
           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           nlines = FILE_LINES(file) & data=FLTARR(ncol,nlines)
           OPENR, lun, file, /GET_LUN & READF, lun, data & CLOSE, lun & FREE_LUN, lun
           ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           rrr      = reform(data(0,*))
           press    = reform(data(2,*))
           rho      = reform(data(3,*))
           temp     = reform(data(4,*))

nsmooth=1
;nsmooth=3
;nsmooth=5
nsmooth=7  ;; 1 lev each km interp, 1 lev each 7 km GCM
;nsmooth=10  
nsmooth=15
nsmooth=20
nsmooth=22
z = smooth(z,nsmooth,/EDGE_TRUNCATE)
theta = smooth(theta,nsmooth,/EDGE_TRUNCATE) 
windu  = smooth(windu,nsmooth,/EDGE_TRUNCATE)
	;windu  = smooth(windu,60,/EDGE_TRUNCATE) ;; juste pour lower shear

;;; altitude above local surface
z = z - z(0)

;;; brunt vaisala frequency
g=3.72 & dlnthetadz=DERIV(z,alog(theta)) & n2=g*dlnthetadz
st2 = n2 / (windu^2) 
	;;w = where(n2 le 0.) 
  	;;n = n2 & n[w] = 0. & n = sqrt(TEMPORARY(n)) & n[w] = -1.
	;;st2[w] = 0. 

;;; shear term
dudz=DERIV(z,windu)
d2udz2=DERIV(z,dudz)
dlnrhodz = DERIV(z,alog(rho))
sh2 = - d2udz2 / windu + dlnrhodz * dudz / windu

;;;; scale height term
drhodz = DERIV(z,rho)
d2rhodz2 = DERIV(z,drhodz)
hh2 = d2rhodz2 / 2. / rho - 0.75 * dlnrhodz * dlnrhodz

	;for i=0,n_elements(windu)-1 do begin
	;  ;print, windu(i), n2 / (windu^2), - d2udz2(i) / windu(i), dlnrhodz(i) * dudz(i) / windu(i), d2rhodz2 / 2. / rho, - 0.75 * dlnrhodz * dlnrhodz
	;  print, 'YO YO ... ', st2(i), sh2(i), hh2(i)
	;endfor
	;stop

;;; scorer parameter ;; erreur dans article Tobie, Lott98 OK
scorer2 = st2 + sh2 + hh2    

;;; relative contributions 
part_buoy = 100. * abs(st2) / (abs(st2)+abs(sh2)+abs(hh2))
part_shea = 100. * abs(sh2) / (abs(st2)+abs(sh2)+abs(hh2))
part_heig = 100. * abs(hh2) / (abs(st2)+abs(sh2)+abs(hh2))

;;; scorer wavelength
scorer = sqrt(scorer2)
lscorer = scorer & w = where(FINITE(lscorer) eq 0) & if (w(0) ne -1) then lscorer[w] = !Values.F_NAN 
lscorer = 2 * !pi / lscorer

        for i=0, n_elements(z)-1, 5 do begin
         print, 'alt, buoy, shear, sc. height, scorer, scorer wvl'
         print, z(i), part_buoy(i), part_shea(i), part_heig(i), scorer2(i), lscorer(i)/1000.
        endfor

     ;;;;;;;;;;;
     what_I_plot = lscorer / 1000.
     column = z / 1000.
     overplot_column = column 
     alt = [0.,150.]
     mention = ''
     minfield_init = 0.001 
     maxfield_init = 200. 
     title_axis = ['Scorer wavelength (km)','Altitude above local surface (km)']
     ;;;;;;;;;;;

!p.multi=[0,2,1]

	;plot, temp, column, $
	;      xtickinterval=50., $
	;      ytickinterval=10., $
	;      xtitle=title_axis(0), $
	;      ytitle=title_axis(1), $
	;      xrange=[50.,200.], $
	;      yrange=alt
	;oplot, 200.+windu, column, linestyle=1

plot, what_I_plot, column, $
      xtickinterval=50., $
      ytickinterval=10., $
      xtitle=title_axis(0), $
      ytitle=title_axis(1), $
      xrange=[minfield_init,maxfield_init], $
      yrange=alt
oplot, temp, column, linestyle=2
oplot, 150.+windu, column, linestyle=1
xyouts, 40., 20., 'dotted: 150+U(z) [m s!U-1!N]', charsize=0.8
xyouts, 40., 10., 'dashed: T(z) [K]', charsize=0.8

plot, part_buoy, column, $
      xtickinterval=20., $
      ytickinterval=10., $
      xtitle='Contributions to Scorer parameter (%)', $
      ytitle=title_axis(1), $
      xrange=[minfield_init,100.], $
      yrange=alt, $
      linestyle=0
oplot, part_shea, column, linestyle=2
oplot, part_heig, column, linestyle=1
xyouts, 30., 140., 'full: buoyancy';, charsize=0.8
xyouts, 30., 130., 'dashed: shear';, charsize=0.8
xyouts, 30., 120., 'dotted: sc. height';, charsize=0.8



PS_End, /PNG
end
