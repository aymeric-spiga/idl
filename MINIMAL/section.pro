pro section, $
	what_I_plot, $				; 2D field
	space, $				; horizontal coordinate
	altitude, $				; altitude coordinate
	minfield=minfield_init, $		; minimum value of plotted field (=0: calculate)
	maxfield=maxfield_init, $		; maximum value of plotted field (=0: calculate)
	minspace=minspace, $			; minimum value of space window (=0: calculate)
	maxspace=maxspace, $			; maximum value of space window (=0: calculate)
	overcontour=overcontour, $		; another 2D field to overplot with contour lines (=0: no)
	overvector_x=overvector_x, $		; wind vector - x component (=0: no)
	overvector_y=overvector_y, $		; wind vector - y component (=0: no)
	colors=colors, $			; number of colors/levels (32 is default)
	title_plot=title_user, $		; title of the plot ('Profile' is default)
	title_axis=title_axis, $		; title of the [x,y] axis (['Field','Altitude'] is default)
	ct=pal, $				; color table (33-rainbow is default)
topo=topography, $
	format=format                           ; format of colorbar annotations ('(F6.2)' is default)

	
;--------------------------------------------
;
; A general routine to plot H/V section maps
;
; USE:	- section, what_I_plot 	
;	- section, what_I_plot, space, altitude
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
flag_cb='true'
SPAWN, 'touch param_plot.idl'
space_save=space
@param_plot.idl
;
;------------------



;------------------
; user settings
;------------------
if (n_elements(title_user) eq 0) then title='Section' else title=title_user
if (n_elements(subtitle_user) eq 0) then subtitle_user=''
if (n_elements(title_axis) ne 2) then begin
	xtitle='Horizontal coordinate'
	ytitle='Altitude'
endif else begin
	xtitle=title_axis(0)
	ytitle=title_axis(1)
endelse	

if (n_elements(space) eq 0) then begin
	space=findgen(n_elements(what_I_plot(*,0)))
endif
if (n_elements(altitude) eq 0) then begin
	altitude=findgen(n_elements(what_I_plot(0,*)))
endif

if (n_elements(minfield_init) eq 0) then minfield_init=0.
if (n_elements(maxfield_init) eq 0) then maxfield_init=0.
if (n_elements(minspace) eq 0) then minspace=0.
if (n_elements(maxspace) eq 0) then maxspace=0.
if (n_elements(minalt) eq 0) then minalt=min(altitude)
if (n_elements(maxalt) eq 0) then maxalt=max(altitude)
if (n_elements(colors) eq 0) then colors=32
if (n_elements(title) eq 0) then title=''
if (n_elements(pal) eq 0) then pal=33
if (n_elements(format) eq 0) then format='(F6.2)'
if (format eq '') then format='(F6.2)'

if ((n_elements(overvector_x) eq 0) or (n_elements(overvector_y) eq 0)) then begin
	overvector_x=0.
	overvector_y=0.
endif
if (n_elements(overcontour) eq 0) then overcontour=0.


;
; dilatation factor (for winds mainly)
;
;if (n_elements(overvector_x)*n_elements(overvector_y) ne 1) then begin

if (n_elements(spacekm) eq 0) then spacekm='false'

if (n_elements(factor) ne 0.) then begin
	if (factor eq 0.) then stop
	space=space/factor
	if (spacekm ne 'false') then xtitle=xtitle+' (km)'
	xtitle=xtitle+' - divided by '+string(factor,'(I0)')
endif else begin
	factor=1.
endelse	
;endif
	
;------------------
; limits
;------------------
;if (minfield_init*maxfield_init eq 0) then begin   
if (minfield_init eq maxfield_init) then begin
	;---different min/max for each layer
	w=where(abs(what_I_plot) lt missing_value)
	if (w(0) ne -1) then begin
		minfield=min(what_I_plot[w])
		maxfield=max(what_I_plot[w])
	endif else begin
		what_I_plot=0.
		print, 'no valid values in this field' 
		stop
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

if (spacekm ne 'true') then begin
	w=where(space eq min(space))
	if (w(0) ne 0) then begin
		space[0:w(0)-1]=space[0:w(0)-1]-360
	endif
endif	


;---------------
; adapt ticks
;---------------
if ((minspace eq 0) and (maxspace eq 0)) then begin 
	minspace=min(space)
	maxspace=max(space)
endif
if (n_elements(intervalx) eq 0) then intervalx=round((maxspace-minspace)/8.)
if (n_elements(intervaly) eq 0) then intervaly=round((maxalt-minalt)/5.)


;------------------
; plot window
;------------------

	; to ensure the right limits
	dumb_what_I_plot=what_I_plot
	dumb_what_I_plot(*,*)=minfield & print, minfield
	dumb_what_I_plot(0,0)=maxfield & print, maxfield


;xticks=sindgen(21)
;xticks(*)=' '
;;xticks(0)='-120'
;;xticks(2)='-119'
;;xticks(4)='-118'
;;xticks(6)='-117'
;;xticks(8)='-116'
;;xticks(10)='-115'
;;xticks(12)='-114'
;;xticks(14)='-113'
;;xticks(16)='-112'
;;xticks(18)='-111'
;;xticks(20)='-110'
;xticks(0)='-110'
;xticks(2)='-109'
;xticks(4)='-108'
;xticks(6)='-107'
;xticks(8)='-106'
;xticks(10)='-105'
;xticks(12)='-104'
;xticks(14)='-103'
;xticks(16)='-102'
;xticks(18)='-101'
;xticks(20)='-100'
;print, xticks

loadct,0,/silent  
contour, $
	/nodata, $ 
	dumb_what_I_plot, $
	space, $
	altitude, $
	title=title, $
	xtitle=xtitle, $
	xrange=[minspace,maxspace], $
	ytitle=ytitle, $
	yrange=[minalt,maxalt], $
	xtickinterval=intervalx, $
	ytickinterval=intervaly, $
;	position=[0.15, 0.15, 0.95, 0.75], $ ;;OLD OLD
	position=[0.15, 0.35, 0.95, 0.95], $
	subtitle=subtitle_user, $
;/isotropic, $
xtickname=xticks, $
	color=0

;if (n_elements(topography) gt 1) then begin
;oplot, space, topography/1000.
;print, max(topography)
;endif

;------------------
; plot field
;------------------

; to avoid spurious blank zones
nlevels=colors
levp=minfield+(findgen(nlevels)/float(nlevels-1))*(maxfield-minfield)

loadct,pal,/silent
;reverse_ct
contour, what_I_plot, $
	space, $
	altitude, $
;	nlevels=colors, $
	levels=levp, $
	/cell_fill, $
	max_value=maxfield, $
	min_value=minfield, $	
	/overplot


;------------------
; colorbar
;------------------

;if (flag_cb eq 'true') then begin
;        format_colorbar=format
;        colorbar, $
;                position=[0.15, 0.85, 0.95, 0.90], $
;                divisions=8, $
;                range=[float(minfield),float(maxfield)], $
;                format=format_colorbar
;endif
;goto, skip_new


if (flag_cb eq 'true') then begin

format_colorbar=format
ndiv=10
;ndiv=8
pos=[0.15, 0.17, 0.95, 0.20]
if (subtitle_user eq '') then pos=[0.15, 0.20, 0.95, 0.23]

;colorbar = Obj_New("COLORBAR",$
;        Range=[minfield,maxfield],$
;        Ncolors=255,$
;        charsize=0.85,$
;        format=format_colorbar,$
;        major=ndiv,$
;        ticklen=-0.15,$
;	position=pos)
;colorbar->Draw
;Obj_Destroy, colorbar

        colorbar, $
;                /vertical, $
                position=pos, $
                range=[float(minfield),float(maxfield)], $
                divisions=ndiv, $
                ncolors=255,$
;               charsize=0.85,$
                ticklen=-0.15,$
                format=format_colorbar

endif

if (n_elements(nam2) ne 0) then thedate=STRSPLIT(nam2, /EXTRACT) & thedate=STRSPLIT(thedate(2),"'",/EXTRACT) & thedate=STRSPLIT(thedate(0),"_",/EXTRACT) & thedate=STRSPLIT(thedate(0),"-",/EXTRACT) & xyouts, 0.95, pos(1)-0.08, 'LMD Martian Mesoscale Model - d'+thedate(2)+'/m'+thedate(1), charsize=0.7, alignment=1., /NORMAL

skip_new:

;--------------------
; overplot contour
;--------------------
if (n_elements(overcontour) ne 1) then begin

	w=where(abs(overcontour) lt missing_value)
	if (w(0) ne -1) then begin
		min_contour=min(overcontour[w])
		max_contour=max(overcontour[w])
	endif


; to avoid spurious blank zones
nlevels=colors/2
if (n_elements(lev) eq 0) then lev=min_contour+(findgen(nlevels)/float(nlevels-1))*(max_contour-min_contour)
nlevels=n_elements(lev)

;lev=round(lev)
;yeah=UNIQ(lev) & lev=lev[yeah] 


loadct,0,/silent
contour, $
	overcontour, $
	space,altitude, $
;	nlevels=colors/2, $
	levels=lev, $
	C_LINESTYLE = (lev LT 0.0), $	;; dotted lines if negative
	max_value=max_contour, $
	min_value=min_contour, $	
	c_labels=findgen(nlevels), $
	color=0, $ ;255, $ ;0, $
	/noerase, $
	/overplot
;        xtickinterval=intervalx, $
;        ytickinterval=intervaly, $
;;/isotropic, $
;;xtickname=xticks, $
;	xtitle=xtitle, $
;	xrange=[minspace,maxspace], $
;	ytitle=ytitle, $
;	yrange=[minalt,maxalt], $
;position=[0.15, 0.15, 0.95, 0.75]



endif

;--------------------
; overplot wind
;--------------------
if (n_elements(overvector_x)*n_elements(overvector_y) ne 1) then begin

	w=where(abs(overvector_x) ge missing_value)  ; u and v have the same missing values ...
	if (w(0) ne -1) then begin
		overvector_x[w]=0.
		overvector_y[w]=0.
	endif
		

len=2.0
if (n_elements(windex) eq 0) then begin
	ref_mag=round(max(sqrt(overvector_x^2+overvector_y^2))/2)
	;ref_mag=round(median(sqrt(overvector_x^2+overvector_y^2)))
ref_mag=15.	;;ce qui est ci dessus deconne !!!
endif else begin
        ref_mag=round(windex)
endelse	
mytext=STRTRIM(STRING(ref_mag),2)+'ms!E-1!N'

mytext=''

if (n_elements(stride) eq 0) then begin
stride=3                ; STRIDE: Pick every nth vector for plotting. (1)
stride=2
stride=5
endif else begin
stride=round(stride)
endelse


ref_pos=[0.10, 0.90] 	
ref_pos=[0.20, 0.05]

myu=overvector_x
myv=overvector_y

longs=space
lats=altitude

wlon=where( (space gt minspace) and (space lt maxspace) )
wlat=where( (altitude gt minalt) and (altitude lt maxalt) )

if ((wlon(0) eq -1) or (wlat(0) eq -1)) then begin
        print, 'wrong window settings'
 	print, space
	print, 'lim ',minspace,maxspace
        print, altitude
        print, 'lim ',minalt,maxalt
 	stop
endif

myu=dblarr(n_elements(space),n_elements(altitude))
myv=dblarr(n_elements(space),n_elements(altitude))

myu(*,*)=!Values.F_NAN
myv(*,*)=!Values.F_NAN

for i=0,n_elements(wlon)-1 do begin
for j=0,n_elements(wlat)-1 do begin
        myu(wlon(i),wlat(j))=overvector_x(wlon(i),wlat(j))
        myv(wlon(i),wlat(j))=overvector_y(wlon(i),wlat(j))
endfor
endfor


	;; eliminate vector for colorbar clarity ?
	;nxs=1
	;nxe=n_elements(longs(*))-2
	;nys=1
	;nye=n_elements(lats(*))-2
	;myu=overvector_x(nxs:nxe,nys:nye)
	;myv=overvector_y(nxs:nxe,nys:nye)
	;longs=space(nxs:nxe)
	;lats=altitude(nys:nye)


;;---------------------------
;;
;factor=float(mean(overvector_x^2))/float(mean(overvector_y^2))
;factor=sqrt(factor)
;factor=factor/2
;
;factor=6.
;factor=12.
;factor=15.
;
;print, 'X/Z correction factor', factor
;
;;rapport ordre de grandeur entre mouvements verticaux et horizontaux
;;
;;---------------------------


myv=myv*factor


loadct,0,/silent
VECTOR, myu, myv, longs, lats, LENGTH=len, REF_MAG=ref_mag,$
        STRIDE=stride, TYPE=2, HEAD_LEN=0.3, ANGLE=30, $
;	REF_POS=ref_pos, $
;	REF_TEXT=mytext,$
        ALIGN=0.5, COLOR=1, /OVERPLOT;, /NOERASE

endif


endelse

space=space_save	;; careful with sequential call to section.pro
			;; combined with modified space (e.g. with factor)
			;; -- this save fix the issue


end
