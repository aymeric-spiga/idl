pro hovmuller, $
	what_I_plot, $				; 2D field
	space, $				; space coordinate
	time, $					; time coordinate
	minfield=minfield_init, $		; minimum value of plotted field (=0: calculate)
	maxfield=maxfield_init, $		; maximum value of plotted field (=0: calculate)
	minspace=minspace, $			; minimum value of space window (=0: calculate)
	maxspace=maxspace, $			; maximum value of space window (=0: calculate)
	overcontour=overcontour, $		; another 2D field to overplot with contour lines (=0: no)
	ct=pal, $				; color table (33-rainbow is default)
	colors=colors, $			; number of colors/levels (32 is default)
	title=title				; title of the plot ('' is default)

;--------------------------------------------
;
; A general routine to plot hovmuller maps
;
; USE:	- hovmuller, what_I_plot 	
;	- hovmuller, what_I_plot, space, time
;	- ... and see options above
;
;--------------------------------------------
; A. Spiga, September 2007
;--------------------------------------------




;---------------------------
; a few permanent settings
;---------------------------
;
missing_value=1.e30
@param_plot.idl
;
;------------------


;------------------
; user settings
;------------------
if (n_elements(space) eq 0) then begin
	space=findgen(n_elements(what_I_plot(*,0)))
endif
if (n_elements(time) eq 0) then begin
	time=findgen(n_elements(what_I_plot(0,*)))
endif

if (n_elements(minfield_init) eq 0) then minfield_init=0.
if (n_elements(maxfield_init) eq 0) then maxfield_init=0.
if (n_elements(minspace) eq 0) then minspace=0.
if (n_elements(maxspace) eq 0) then maxspace=0.
if (n_elements(colors) eq 0) then colors=32
if (n_elements(title) eq 0) then title=''
if (n_elements(pal) eq 0) then pal=33
if (n_elements(overcontour) eq 0) then overcontour=0.
if (n_elements(titlex) eq 0) then titlex='space'
if (n_elements(titley) eq 0) then titley='time'

;------------------
; limits
;------------------
if (minfield_init*maxfield_init eq 0) then begin   
	;---different min/max for each layer
	w=where(abs(what_I_plot) lt missing_value)
	if (w(0) ne -1) then begin
		minfield=min(what_I_plot[w])
		maxfield=max(what_I_plot[w])
	endif else begin
		what_I_plot=0.
		print, 'no valid values in this field' 
	endelse
endif else begin
	;---user-defined limits	
	minfield=minfield_init 
	maxfield=maxfield_init
	;;w=where(what_I_plot lt minfield) & if (w(0) ne -1) then what_I_plot[w]=minfield
	;;w=where(what_I_plot gt maxfield) & if (w(0) ne -1) then what_I_plot[w]=maxfield
endelse

if ((minfield+maxfield eq 0) and (abs(minfield) ne abs(maxfield))) then begin
	print, 'nothing to plot ... skipping this plot ...'
endif else begin



;-----------------------------------------
; fix the possible longitude ugliness :)
; ...get rid of the -180/180 limit
;-----------------------------------------

w=where(space eq min(space))
if (w(0) ne 0) then begin
	space[0:w(0)-1]=space[0:w(0)-1]-360
endif


;---------------
; adapt ticks
;---------------
if ((minspace eq 0) and (maxspace eq 0)) then begin
	minspace=min(space)
	maxspace=max(space)
endif
intervalx=round((maxspace-minspace)/8.)
intervaly=round((max(time)-min(time))/5.)
intervalx=60.
intervaly=4.


;------------------
; plot window
;------------------

	; to ensure the right limits
	dumb_what_I_plot=what_I_plot
	dumb_what_I_plot(*,*)=minfield & print, minfield
	dumb_what_I_plot(0,0)=maxfield & print, maxfield


loadct,0,/silent  
contour, $
	/nodata, $ 
	dumb_what_I_plot, $
	space, $
	time, $
	title=title, $
	xtitle=titlex, $
	xrange=[minspace,maxspace], $
	ytitle=titley, $
	yrange=[min(time),max(time)], $
	xtickinterval=intervalx, ytickinterval=intervaly, $
;	position=[0.15, 0.15, 0.95, 0.75], $
	color=0



;------------------
; plot field
;------------------

;lim=45.
;ecart=5.
;yeah=-lim + ecart*findgen(lim*2./5.)

loadct,pal,/silent
contour, what_I_plot, $
	space, $
	time, $
	nlevels=colors, $
;lev=yeah, $
	/cell_fill, $
	max_value=maxfield, $
	min_value=minfield, $	
	/overplot


;------------------
; colorbar
;------------------

;;format_colorbar='(F6.2)'
;;format_colorbar='(F10.3)'
;;format_colorbar='(F12.2)'
;format_colorbar='(F4.0)'  ;temperature
;colorbar, $
;;	position=[0.15, 0.85, 0.95, 0.90], $
;	/vertical, $
;	divisions=8, $
;	range=[minfield,maxfield], $
;	format=format_colorbar

;--------------------
; overplot contour
;--------------------
if (n_elements(overcontour) ne 1) then begin

	w=where(abs(overcontour) lt missing_value)
	if (w(0) ne -1) then begin
		min_contour=min(overcontour[w])
		max_contour=max(overcontour[w])
	endif

loadct,0,/silent
contour, overcontour, $
	space,time, $
	nlevels=colors/2, $
	max_value=max_contour, $
	min_value=min_contour, $	
	c_labels=findgen(colors/2), $
	color=0, $
	/noerase, $
	xtitle='space', $
	xrange=[min(space),max(space)], $
	ytitle='time', $
	yrange=[min(time),max(time)], $
	xtickinterval=intervalx, ytickinterval=intervaly, $
	position=[0.15, 0.15, 0.95, 0.75]

endif


endelse

end
