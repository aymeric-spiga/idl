pro grad, d_param, param, x, y, flag

; pour calculer le champ de gradient d'un champ 2D
;	d_param = résultat 2D
; 	param: tableau 2D
; 	x: tableau 1D 
; 	y: tableau 1D
;	flag: 	1 pour grad_x
;		2 pour grad_y


missing=1.e30
w=where(abs(param) gt missing)
if (w(0) ne -1) then param[w]=0.


imax = n_elements(param(*,0))
jmax = n_elements(param(0,*))

d_param=dblarr(imax,jmax) 
case flag of
1: begin
for j = 0, jmax-1 do begin
		d_param(*,j) = DERIV(x,param(*,j))
;               d_param(*,j) = DERIV(reform(param(*,j)),x)
endfor
end
2: begin
for i = 0, imax-1 do begin
		d_param(i,*) = DERIV(y,param(i,*))
;               d_param(i,*) = DERIV(reform(param(i,*)),y)
endfor
end
endcase
end
