pro model_levels, plotlev
;; first is one (as in gw.def)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
file = './zefolder/'
file = file + 'wrfout_d01_2024-03-21_00:00:00'

file = './wrfout_d01_9999-09-09_09:00:00' 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;
;; ini
;;
;if (n_elements(plotlev) eq 0) then plotlev=7

;set_plot, 'ps'
;device, filename='model_levels.ps'

;
; constantes
;
grav = 3.72

;
; retrieve model geopotential
;
openr,unit,'heights.dat',/get_lun,error=err
;IF (err ne 0) THEN BEGIN
	getcdf, $
	        file=file, $
	        charvar='PHTOT', $
	        invar=geop
        getcdf, $
                file=file, $
                charvar='HGT', $
                invar=hgt
;	save, geop, filename='heights.dat'
;ENDIF ELSE BEGIN
;	restore, filename='heights.dat'
;ENDELSE

;
; size
;
s = size(geop)
nx = s[1] ;& print, nx 
ny = s[2] ;& print, ny
nz = s[3] ;& print, nz
nt = s[4] ;& print, nt

if (n_elements(plotlev) eq 0) then plotlev=nz-1

for lev = 1, nz-1 do begin

print, '*************************************************************************'
print, 'model level ', lev
 
;
; heights of model levels
;
;height = (geop(*,*,lev-1,*) - geop(*,*,0.,*)) / grav
height = fltarr(nx,ny,nt)
for i=0,nt-1 do height(*,*,i) = geop(*,*,lev-1,i) / grav - hgt(*,*)

;
; spatial mean (temporal variations)
;
yes = TOTAL(TOTAL(height,1),1) / float(nx) / float(ny)
print, 'spat. mean (temp. var.) ', mean(yes), max(yes), min(yes);, max(yes) - min(yes)

;
; temporal mean (spatial variations)
;
yet = TOTAL(height,3) / float(nt)  
print, 'temp. mean (spat. var.) ', mean(yet), max(yet), min(yet);, max(yet) - min(yet)

;
; complete mean
;
ye = TOTAL(TOTAL(TOTAL(height,1),1),1) / float(nx) / float(ny) / float(nt)
print, '******* mean ', ye 

if (lev eq plotlev) then begin
	!p.multi=[0,2,1]
	contour, yet, /cell_fill, nlevels=30, title=string(min(yet))+string(max(yet))
	plot, yes, yrange=[min(yes),max(yes)]
	stop
endif

endfor

;device, /close
end
