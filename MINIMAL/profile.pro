pro profile, $
	what_I_plot, $				; 1D vertical profile
	column, $				; altitudes	
	alt=alt, $				; altitude range [altmin, altmax]
	minfield=minfield_init, $		; minimum value of plotted field (=0: calculate)
	maxfield=maxfield_init, $		; maximum value of plotted field (=0: calculate)
	inprofile=overplot, $			; another vertical profile to overplot
	incolumn=overplot_column, $		; altitudes of the other vertical profile (in case /= column)
	discrete=discrete, $			; show the profile points (= type of points in !psym)
	title_plot=title_user, $		; title of the plot ('Profile' is default)
	title_axis=title_axis, $		; title of the [x,y] axis (['Field','Altitude'] is default)
	mention=mention				; add text precision within the plot window (default is nothing or '')

;--------------------------------------------
;
; A general routine to plot 1D profiles
;
; USE: see comments above
;
; EXAMPLE:
;	profile, findgen(30)^2, findgen(30)
;	profile, findgen(30)^2, findgen(30), overplot=findgen(30)^2.2, discrete=4, title_plot='An example of profile'
;	profile, findgen(30)^2, findgen(30), title_axis=['Dumb profile', 'Dumb Altitude']
;
;--------------------------------------------
; A. Spiga, September 2007
;--------------------------------------------



;---------------------------
; a few permanent settings
;---------------------------
;
missing_value=1.e30
SPAWN, 'touch param_plot.idl'
@param_plot.idl
;
;------------------


;-----------------------------
; user and default settings
;-----------------------------
if (n_elements(title_user) eq 0) then title='Profile'
if (n_elements(title_axis) ne 2) then begin
	xtitle='Field'
	ytitle='Altitude'
endif else begin
	xtitle=title_axis(0)
	ytitle=title_axis(1)
endelse	
if (n_elements(maxfield_init) eq 0) then maxfield_init=0
if (n_elements(minfield_init) eq 0) then minfield_init=0

if (n_elements(column) eq 0) then column=findgen(n_elements(what_I_plot))
if (n_elements(mention) eq 0) then mention='' 
if (mention ne '') then mention='('+mention+')'

;------------------
; missing values
;------------------
w=where(abs(what_I_plot) gt missing_value)
if (w(0) ne -1) then what_I_plot[w]=!Values.F_NAN


;------------------
; limits
;------------------
if (minfield_init*maxfield_init eq 0) then begin 
	xrange=[min(what_I_plot),max(what_I_plot)]
        print, xrange
endif else begin
	xrange=[minfield_init,maxfield_init]
        print, xrange
endelse

if (xrange(0)*xrange(1) eq 0) then begin
	print, 'nothing to plot ... skipping this plot ...'
	stop
endif

  
;---------------
; adapt ticks
;---------------
;intervalx=round((xrange(1)-xrange(0))/5.)
;intervaly=round(abs(max(column)-min(column))/5.)
;intervalx=5.
;intervaly=5.

;if (n_elements(intervalx) eq 0) then intervalx=round((xrange(1)-xrange(0))/8.)
;;if (n_elements(intervaly) eq 0) then intervaly=round((abs(max(column)-min(column))/5.))


;; ----------------------------
;; manage the altitude range
;; ----------------------------
;if (n_elements(alt) eq 2) then begin
;	altmin=alt(0)
;	altmax=alt(1)
;	w=where((column ge altmin) and (column le altmax))
;	if (n_elements(intervaly) eq 0) then intervaly=(altmax-altmin)/5
;	if (w(0) ne -1) then begin
;		column=column[w]
;		what_I_plot=what_I_plot[w]
;		minfield=min(what_I_plot)
;		maxfield=max(what_I_plot)	
;	endif
;endif



;-------------
; plot
;-------------
plot, $
	what_I_plot,column, $
	xrange=xrange, $
;	yrange=[min(column),max(column)], $
yrange=alt, $
	xtickinterval=intervalx, $
	ytickinterval=intervaly, $	
	xtitle=xtitle, $
	ytitle=ytitle,$
	title=title_user,$
;/ylog, $
	subtitle=mention 

;-------------
; overplot
;-------------
if (n_elements(discrete) ne 0) then begin
if (discrete ne 0) then begin
	!psym=discrete
	oplot, $
		what_I_plot,column
	!psym=0
endif
endif	

if (n_elements(overplot) gt 1) then begin
	column_overplot=column
	w=where(abs(overplot) gt missing_value)
	if (w(0) ne -1) then overplot[w]=!Values.F_NAN
	if (n_elements(overplot_column) gt 1) then column_overplot=overplot_column
oplot, $
	overplot, column_overplot, $
	linestyle=2
	if (n_elements(discrete) ne 0) then begin
		!psym=discrete
		oplot, $
			overplot,column_overplot
		!psym=0
	endif	
endif	



;;--------------
;; text inside
;;--------------
;if (n_elements(mention) ne 0) then begin
;	alignx=(max(what_I_plot)+min(what_I_plot))/2
;	aligny=max(column)-intervaly/2
;	xyouts, alignx, $
;		aligny, $
;		mention
;endif

end
