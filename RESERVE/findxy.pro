pro findxy, 		$
	x=x, 		$	;;IN
	y=y, 		$	;;IN
	point=point,	$	;;IN
	ind=ind, 	$	;;OUT		
	flag=flag		;;IN OPT



pointint=[0,0]

;;;if (point eq 'VL1') then pointint=[-48.22,22.70]
;;;if (point eq 'VL2') then pointint=[133.74,47.93]
;;;if (point eq 'PAT') then pointint=[-33.5,19.1]

if (n_elements(point) eq 2) then begin
	pointint=point
endif else begin
;;Si lon=312.050381 et lat=22.269628, alt=-3637.1396
;;Si lon=134.282101 et lat=47.668093, alt=-4497.0474
;;Si lon=326.747364 et lat=19.099641, alt=-3681.9993
if (point eq 'VL1') then pointint=[312.050381-360.,22.269628]
if (point eq 'VL2') then pointint=[134.282101,47.668093]
if (point eq 'PAT') then pointint=[326.747364-360.,19.099641]
endelse



if (n_elements(x) eq 0) then print, 'set x='
if (n_elements(y) eq 0) then print, 'set y='
if (n_elements(flag) eq 0) then flag=1


x=reform(x)
y=reform(y)

ind=intarr(2)



CASE (flag) OF
1: begin

	minim=abs(x-pointint(0))
	w=where(minim eq min(minim))
	print, 'distance is ', minim(w(0))
	ind[0]=w(0)

        minim=abs(y-pointint(1))
        w=where(minim eq min(minim))
	print, 'distance is ', minim(w(0))
        ind[1]=w(0)

end
2: begin

	distance = sqrt((x - pointint(0))^2 + (y - pointint(1))^2)
	mx = min(distance, location)
	print, 'distance 2D is ',mx
	ind = array_indices(distance, location)
	print, 'on ... ', pointint(0), pointint(1)
	print, 'found ... ', x[ind(0),ind(1)], y[ind(0),ind(1)]

end
endcase

end



