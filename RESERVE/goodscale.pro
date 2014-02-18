;;Fixe l'echelle pour une tab2diable de type (x,y,z,t) ou (x,y,t) ou (z,t)
	;;Possibilite de restreindre a une zone et des heures donnees

pro goodscale, $
    tab2d, $
    mini, $
    maxi, $
    coeff=coeff_sigma
	
print,'Yeah !'

;
; input: the 2D array which is ready to be plotted
;


;if (n_elements(coeff_sigma) eq 0) then coeff_sigma=3.3 ;; vincent
if (n_elements(coeff_sigma) eq 0) then coeff_sigma=2.5 ;; aymeric

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	
;;calcul de la moyenne
	N=n_elements(tab2d)
	moy=(1/float(N))*total(tab2d) 	

;;calcul de la tab2diance
	moyarray = make_array(N, 1, /integer, value = 1)	
	varia=(1/float(N-1))*total( ( tab2d - (moyarray * moy) )^2 )
	
;;calcul de sigma
	sigma = sqrt( varia )
	
;;echelle
	mini = moy - float(coeff_sigma) * sigma & stratmin = 'sigma'
	maxi = moy + float(coeff_sigma) * sigma & stratmax = 'sigma'

;;au cas ou 
if (mini lt min(tab2d)) then begin
	mini=min(tab2d)
	stratmin = 'real min'
endif
if (maxi gt max(tab2d)) then begin
	maxi=max(tab2d)
	stratmax = 'real max'
endif
	
	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
;;affichage a l'ecran :
print,'min :'+STRTRIM(mini,2)+' '+stratmin	
print,'max :'+STRTRIM(maxi,2)+' '+stratmax


end
