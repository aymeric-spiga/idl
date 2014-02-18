;pro try


restore, filename='getturb.dat'

;;wt, tke, ztke, h, ht, t, p, pt, stst, localtime, xtke, ytke, wmax, wmin, depressions, psmin, filename='getturb.dat'


nz = n_elements(wt(*,0))
nt = n_elements(wt(0,*))

;
; nz: colonnes
; nt: lignes
;

openw, lun, 'LES_MERIDIANI_vert_levels.txt', /GET_LUN
printf, lun, h, format='(F12.6)'
close, lun
free_lun, lun

openw, lun, 'LES_MERIDIANI_local_times.txt', /GET_LUN
printf, lun, localtime, format='(F12.6)'
close, lun
free_lun, lun

openw, lun, 'LES_MERIDIANI_TKE.txt', /GET_LUN
printf, lun, tke, format='('+string(nz,'(I0)')+'F12.6)'
close, lun
free_lun, lun

openw, lun, 'LES_MERIDIANI_xTKE.txt', /GET_LUN
printf, lun, xtke, format='('+string(nz,'(I0)')+'F12.6)'
close, lun
free_lun, lun

openw, lun, 'LES_MERIDIANI_yTKE.txt', /GET_LUN
printf, lun, ytke, format='('+string(nz,'(I0)')+'F12.6)'
close, lun
free_lun, lun

openw, lun, 'LES_MERIDIANI_zTKE.txt', /GET_LUN
printf, lun, ztke, format='('+string(nz,'(I0)')+'F12.6)'
close, lun
free_lun, lun

openw, lun, 'LES_MERIDIANI_Wmax.txt', /GET_LUN
printf, lun, wmax, format='('+string(nz,'(I0)')+'F12.6)'
close, lun
free_lun, lun

openw, lun, 'LES_MERIDIANI_Wmin.txt', /GET_LUN
printf, lun, wmin, format='('+string(nz,'(I0)')+'F12.6)'
close, lun
free_lun, lun

openw, lun, 'LES_MERIDIANI_Tpot.txt', /GET_LUN
printf, lun, t, format='('+string(nz,'(I0)')+'F12.6)'
close, lun
free_lun, lun

restore, filename='addturb.dat'
openw, lun, 'LES_MERIDIANI_VELmax.txt', /GET_LUN
printf, lun, VELMAX, format='('+string(nz,'(I0)')+'F12.6)'
close, lun
free_lun, lun



;yeye = fltarr(nz,nt)
;openr, lun, 'LES_MERIDIANI_Tpot.txt', /GET_LUN
;readf, lun, yeye
;close, lun
;free_lun, lun
;contour, transpose(yeye), nlevels=20

yeye2 = fltarr(nz,nt)
openr, lun, 'LES_MERIDIANI_TKE.txt', /GET_LUN
readf, lun, yeye2
close, lun
free_lun, lun
contour, transpose(yeye2), nlevels=20


;end
