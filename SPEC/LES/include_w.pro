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
        for loop  = 0, nf-1 do begin
                                                  timetime = SYSTIME(1)
          if (loop ne nf-1) then nloop2=nt else nloop2=ntlast
        for loop2 = 0, nloop2-1 do begin
                 ;wprime = getget(filesWRF(loop), 'W', count=[0,0,0,1], offset=[0,0,0,loop2])

;if (localtime(yeye) eq 11.) then begin
if (localtime(yeye) eq 14.) then begin
wprime = getget(filesWRF(loop), 'W', count=[0,0,0,1], offset=[0,0,0,loop2])
print, 'save !! ', localtime(yeye)
save, wprime, h, filename='savew'
stop
uprime = getget(filesWRF(loop), 'U', count=[0,0,0,1], offset=[0,0,0,loop2])
vprime = getget(filesWRF(loop), 'V', count=[0,0,0,1], offset=[0,0,0,loop2])
help, uprime
help, vprime
veltot = sqrt(uprime^2+vprime^2)
print, 'save !! ', localtime(yeye)
save, veltot, h, filename='saveu'
stop
endif

        yeye = TEMPORARY(yeye) + 1
        endfor
        endfor
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

