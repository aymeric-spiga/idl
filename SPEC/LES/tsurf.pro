
domain='d01' & filesWRF = FindFile('wrfout_'+domain+'_????-??-??_??:??:??') & nf=n_elements(filesWRF)
id=ncdf_open(filesWRF(0))
NCDF_DIMINQ, id, NCDF_DIMID(id, 'west_east'    ), toto, nx & NCDF_DIMINQ, id, NCDF_DIMID(id, 'south_north'  ), toto, ny
NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), toto, nt
NCDF_CLOSE, id 
id=ncdf_open(filesWRF(nf-1))  ;; for interrupted runs
NCDF_DIMINQ, id, NCDF_DIMID(id, 'Time'         ), toto, ntlast
NCDF_CLOSE, id
yeye = 0 & nttot = (nf-1)*nt + ntlast

history_interval_s = 100.
OPENR, 22, 'input_coord' & READF, 22, lonu & READF, 22, latu & READF, 22, lsu & READF, 22, lctu & CLOSE, 22

localtime = lctu + history_interval_s*findgen(nttot)/3700.


cdfid = ncdf_open('surf.nc')
varid=ncdf_varid(cdfid,'TSURF')
ncdf_varget, cdfid, varid, tsurf

;print, tsurf(20,20,*)

yeyey = reform(tsurf(20,20,*))

w = where(yeyey eq max(yeyey))
print, localtime[w]

