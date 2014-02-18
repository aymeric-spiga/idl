        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        OPENR, 22, 'input_coord' & READF, 22, lonu & READF, 22, latu & READF, 22, lsu & READF, 22, lctu & CLOSE, 22
        OPENR, 23, 'input_more' & READF, 23, hgtu, tsurfu & CLOSE, 23
        domain='d01' & filesWRF = FindFile('wrfout_'+domain+'_????-??-??_??:??:??') & nf=n_elements(filesWRF)
        id=ncdf_open(filesWRF(0))
                NCDF_DIMINQ, id, NCDF_DIMID(id, 'west_east'    ), toto, nx & NCDF_DIMINQ, id, NCDF_DIMID(id, 'south_north'  ), toto, ny
                NCDF_DIMINQ, id, NCDF_DIMID(id, 'bottom_top'   ), toto, nz & NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), toto, nt
                NCDF_CLOSE, id
        id=ncdf_open(filesWRF(nf-1))  ;; for interrupted runs
                NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), toto, ntlast
                NCDF_CLOSE, id
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        yeye = 0 & nttot = (nf-1)*nt + ntlast & localtime = lctu + history_interval_s*findgen(nttot)/3700.
        std_tprime = fltarr(nz,nttot) 
        field = fltarr(nz,nttot)
        anomalt = 1.
        for loop  = 0, nf-1 do begin
                                                  timetime = SYSTIME(1)
          print, loop
          if (loop ne nf-1) then nloop2=nt else nloop2=ntlast
        for loop2 = 0, nloop2-1 do begin
!QUIET=1

 ;;;; TEMPPOT
 ;tprime = getget(filesWRF(loop), 'T', anomaly=anomalt, count=[0,0,0,1], offset=[0,0,0,loop2])   ;; t' = t - <t>
 ;;;; TEMP (K)
 tempk = reform ( ( t0 + getget(filesWRF(loop), 'T', count=[0,0,0,1], offset=[0,0,0,loop2]) ) * ( getget(filesWRF(loop), 'PTOT', count=[0,0,0,1], offset=[0,0,0,loop2]) / p0 )^r_cp )
 meanmean = reform(TOTAL(TOTAL(tempk,1),1) / float(nx) / float(ny))
 for i=0,n_elements(tempk(*,0,0,0))-1 do for j=0,n_elements(tempk(0,*,0,0))-1 do tempk(i,j,*) = tempk(i,j,*) - meanmean
 tprime = TEMPORARY(tempk)^2
 
 std_tprime(*, yeye) = sqrt ( reform(TOTAL(TOTAL(TEMPORARY(tprime),1),1) / float(nx) / float(ny)) )
 ;print, std_tprime(*, yeye)
 ;print, yeye

!QUIET=0
        yeye = TEMPORARY(yeye) + 1
        endfor
        endfor
        ;save, std_tprime, localtime, filename='addturb.dat' ;; pour l'instant ecrase systematiquement
        ;stop
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(std_tprime)
zex           =  localtime
zey           =  h
set_name      =  'std_tprime.ps'
set_title     =  "Temperature standard deviation (K)" ;"Maximum updraft speed (m.s!U-1!N)"
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  ''
set_xrange    =  [8.,18.]
set_yrange    =  [0.,7.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  0.
maxval        =  2. ;15.;14.;13.
nlev          =  10 ;(maxval-minval)*10.
pal           =  22
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
  ;; 2. color field
  loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
              contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zefield       =  transpose(std_tprime)
zex           =  localtime
zey           =  h
set_name      =  'std_tprime_ns.ps'
set_title     =  "Temperature standard deviation (K)" ;"Maximum updraft speed (m.s!U-1!N)"
set_titlex    =  'Local Time (h)'
set_titley    =  'Altitude above surface (km)'
set_subtitle  =  ''
set_xrange    =  [8.,18.]
set_yrange    =  [0.,1.]
set_tickx     =  1.
set_ticky     =  1.
minval        =  0.
maxval        =  2. ;15.;14.;13.
nlev          =  10 ;(maxval-minval)*10.
pal           =  22
rrr           =  'no'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveps eq 'true') then PS_Start, FILENAME=set_name
!P.Charsize = 1.2
;; 0. levels
lev = minval + (maxval-minval)*findgen(nlev+1)/float(nlev) & if (minval ne 0.) then lev = lev[where(lev ne 0.)]
;;; 1. background
loadct, 0 & contour, /NODATA, zefield, zex, zey, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0
  ;; 2. color field
  loadct, pal & if (rrr eq 'yes') then TVLCT, r, g, b, /Get & if (rrr eq 'yes') then TVLCT, Reverse(r), Reverse(g), Reverse(b)
              contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /overplot, /cell_fill
;; 3. contour field
loadct, 0 & contour, zefield, zex, zey, levels=lev, c_labels=findgen(n_elements(lev))*0.+1., /noerase, xtitle=set_titlex, xrange=set_xrange, xtickinterval=set_tickx, ytitle=set_titley, yrange=set_yrange, ytickinterval=set_ticky, title=set_title, subtitle=set_subtitle, color=0, C_LINESTYLE = (lev LT 0.0)
;; 4. choose output
if (saveps eq 'true') then PS_End, /PNG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




stop
