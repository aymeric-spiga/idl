pro map_latlon, $
	what_I_plot, $				; 2D field
	lon, $					; 1D latitude OR 2D
	lat, $					; 1D longitude OR 2D
	minfield=minfield_init, $		; minimum value of plotted field (=0: calculate)
	maxfield=maxfield_init, $		; maximum value of plotted field (=0: calculate)
	overcontour=overcontour, $		; another 2D field to overplot with contour lines (=0: no)
	overvector_x=overvector_x, $		; wind vector - x component (=0: no)
	overvector_y=overvector_y, $		; wind vector - y component (=0: no)
	ct=pal, $				; color table (33-rainbow is default)
	colors=colors, $			; number of colors/levels (32 is default)
	title=title_user, $			; title of the plot ('' is default)
	format=format				; format of colorbar annotations ('(F6.2)' is default)

;--------------------------------------------
;
; A general routine to plot lat/lon maps
;
; USE:	- map_latlon, what_I_plot 	
;	- map_latlon, what_I_plot, lon, lat
;	- ... and see options above
;
;--------------------------------------------
; A. Spiga, August 2007
;--------------------------------------------




;---------------------------
; a few permanent settings
;---------------------------
;
missing_value=1.e30
flag_cb='true'
;-- optional:
;-- user defined parameters
; do not forget to type '.compile map_latlon' to update changes
SPAWN, 'touch param_plot.idl'
@param_plot.idl
;
;------------------


;------------------
; user settings
;------------------
if (n_elements(lon) eq 0) then begin
	lon=findgen(n_elements(what_I_plot(*,0)))
endif
if (n_elements(lat) eq 0) then lat=findgen(n_elements(what_I_plot(0,*)))
if (n_elements(minfield_init) eq 0) then minfield_init=0.
if (n_elements(maxfield_init) eq 0) then maxfield_init=0.
if (n_elements(overcontour) eq 0) then overcontour=0.
if ((n_elements(overvector_x) eq 0) or (n_elements(overvector_y) eq 0)) then begin
	overvector_x=0.
	overvector_y=0.
endif
if (n_elements(colors) eq 0) then colors=32
if (n_elements(title_user) eq 0) then title_user=''
if (n_elements(subtitle_user) eq 0) then subtitle_user=''
if (n_elements(pal) eq 0) then pal=33
if (n_elements(format) eq 0) then format='(F6.2)'
if (format eq '') then format='(F6.2)'
;if (n_elements(xtitleset) eq 0) then xtitleset='longitude'
;if (n_elements(ytitleset) eq 0) then ytitleset='latitude'
if (n_elements(title_axis) ne 2) then begin
        xtitleset='longitude'
        ytitleset='latitude'
endif else begin
        xtitleset=title_axis(0)
        ytitleset=title_axis(1)
endelse
if (n_elements(isotropic) eq 0) then isotropic='true'
if (n_elements(coord2d) eq 0) then coord2d='false'

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
;goodscale, $
;     what_I_plot, $
;     minfield, $
;     maxfield
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

;;;;; PB WITH 2D COORDINATES
;w=where(lon eq min(lon))
;if (w(0) ne 0) then begin
;	lon[0:w(0)-1]=lon[0:w(0)-1]-360
;endif


;---------------
; adapt ticks
;---------------
if (n_elements(windowx) eq 0) then windowx=[min(lon),max(lon)]
if (n_elements(windowy) eq 0) then windowy=[min(lat),max(lat)]
if (n_elements(intervalx) eq 0) then intervalx=round((windowx(1)-windowx(0))/8.)
if (n_elements(intervaly) eq 0) then intervaly=round((windowy(1)-windowy(0))/5.)


;------------------
; plot window
;------------------

	;lons=lon & lats=lat
	;;--------------------------------------------
	;;
	;; ANOMALY
	;;
	;w=where( (lon ge windowx(0)) and (lon le windowx(1)) ) & w2=where( (lat ge windowy(0)) and (lat le windowy(1)) )
	;yeah=fltarr(n_elements(w),n_elements(w2)) & for i=0,n_elements(w)-1 do for j=0,n_elements(w2)-1 do yeah(i,j)=what_I_plot(w(i),w2(j))
	;what_I_plot = yeah - mean(yeah)
	;lon=lon[w] & lat=lat[w2]
	;;--------------------------------------------


	; to ensure the right limits
	dumb_what_I_plot=what_I_plot
	dumb_what_I_plot(*,*)=minfield & print, minfield
	dumb_what_I_plot(0,0)=maxfield & print, maxfield

;;
;;
;subtitle_hunhun = subtitle_user
;subtitle_user = ''
;;
;;


pospos = [0.10, 0.12, 0.90, 0.92] ;; iso
;pospos = [0.15, 0.35, 0.95, 0.95] ;; non-iso

if (isotropic eq 'true') then begin
        loadct, 0, /silent
        if ( coord2d ne 'polar' ) then begin
                ;;trace non polaire
                contour, /nodata, $
                         /isotropic, $
                         position=pospos, $
                         dumb_what_I_plot,      lon,                            lat, $
                         title=title_user,      xtitle=xtitleset,               ytitle=ytitleset, $
                         xrange=windowx,                yrange=windowy, $
                         xtickinterval=intervalx,       ytickinterval=intervaly, $
                         subtitle=subtitle_user, color=0
        endif else begin
                ;;trace avec MAP_SET -- MAP_SET -- MAP_SET 
                print,'------------------- Polar plot'
                contour, /nodata, $
                         /isotropic, $
                         /overplot, $
                         dumb_what_I_plot,      lon,                            lat, $
                         title=title_user,      xtitle=xtitleset,               ytitle=ytitleset, $
                         xrange=windowx,                yrange=windowy, $
                         xtickinterval=intervalx,       ytickinterval=intervaly, $
                         subtitle=subtitle_user, color=0
        endelse
        subtitle_user = ''
endif else begin
        pospos = [0.15, 0.35, 0.95, 0.95] ;; non-iso
	loadct, 0, /silent
        if ( coord2d ne 'polar' ) then begin
                ;;trace non polaire
                contour, /nodata, $
                         position=pospos, $
                         dumb_what_I_plot,      lon,                            lat, $
                         title=title_user,      xtitle=xtitleset,               ytitle=ytitleset, $
                         xrange=windowx,                yrange=windowy, $
                         xtickinterval=intervalx,       ytickinterval=intervaly, $
                         color=0,subtitle=subtitle_user
        endif else begin
                ;;trace avec MAP_SET -- MAP_SET -- MAP_SET 
                print,'------------------- Polar plot'
                contour, /nodata, $
                         /overplot, $
                         dumb_what_I_plot,      lon,                            lat, $
                         title=title_user,      xtitle=xtitleset,               ytitle=ytitleset, $
                         xrange=windowx,                yrange=windowy, $
                         xtickinterval=intervalx,       ytickinterval=intervaly, $
                         color=0,subtitle=subtitle_user
        endelse
endelse	

;------------------
; plot field
;------------------

; to avoid spurious blank zones
nlevels=colors
levp=minfield+(findgen(nlevels)/float(nlevels-1))*(maxfield-minfield)

loadct,pal,/silent
;reverse_ct
if (pal eq 0) then reverse_ct
if (pal eq 21) then reverse_ct
contour, what_I_plot, $
	lon,lat, $
;	nlevels=colors, $
	levels=levp, $
	/cell_fill, $
	max_value=maxfield, $
	min_value=minfield, $
	/overplot


;nlat=128.
;oplot, lat*0.+nlat*100./1000., linestyle=1

;;; ANOMALY
;lon=lons
;lat=lats

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

;;; ndiv peut etre regle
if (n_elements(ndiv) eq 0) then ndiv=10

;if (n_elements(poscb) eq 0) then poscb=0.95
;if (n_elements(poscb) eq 0) then poscb=0.80
if (n_elements(poscb) eq 0) then poscb=0.70

     if (isotropic eq 'true') then begin
       ;poscb=0.85 & 
       pos=[poscb, 0.12, poscb+0.03, 0.92]
       ;poscb=0.95 & pos=[poscb, 0.12, poscb+0.03, 0.92]
;       poscb=0.96 & pos=[poscb, 0.12, poscb+0.03, 0.92]
;       poscb=0.05 & pos=[poscb, 0.12, poscb+0.03, 0.92]
               ;colorbar = Obj_New("COLORBAR",$
               ; Range=[float(minfield),float(maxfield)],$
               ; /vertical, $
               ; Ncolors=255,$
               ; charsize=0.85,$
               ; format=format_colorbar,$
               ; major=ndiv,$
               ; ticklen=-0.15,$
               ; position=pos)
               ;colorbar->Draw
               ;Obj_Destroy, colorbar
        colorbar, $
                /vertical, $
                position=pos, $
                range=[float(minfield),float(maxfield)], $
                divisions=ndiv, $
                ncolors=255,$
;               charsize=0.85,$
                ticklen=-0.15,$
                format=format_colorbar

xyouts, pos(2), 0.01, subtitle_user, charsize=0.8, alignment=1., /NORMAL

	if (n_elements(nam2) ne 0) then begin
		thedate=STRSPLIT(nam2, /EXTRACT)
		thedate=STRSPLIT(thedate(2),"'",/EXTRACT)
		thedate=STRSPLIT(thedate(0),"_",/EXTRACT) 
		thedate=STRSPLIT(thedate(0),"-",/EXTRACT) 
		xyouts, pos(2), 0.01, 'LMD Martian Mesoscale Model - d'+thedate(2)+'/m'+thedate(1), charsize=0.7, alignment=1., /NORMAL 
	endif else begin
		thedate=''
	endelse
;;
;;
;;
;xyouts, pos(2), 0.01, 'LMD Mars Mesoscale Model / '+subtitle_hunhun, charsize=0.95, alignment=1., /NORMAL
;;
;;
;;

     endif else begin

	pos=[0.15, 0.17, 0.95, 0.20]
	if (subtitle_user eq '') then pos=[0.15, 0.20, 0.95, 0.23]

               ;colorbar = Obj_New("COLORBAR",$
               ; Range=[float(minfield),float(maxfield)],$
               ; Ncolors=255,$
               ; charsize=0.85,$
               ; format=format_colorbar,$
               ; major=ndiv,$
               ; ticklen=-0.15,$
               ; position=pos)
               ;colorbar->Draw
               ;Obj_Destroy, colorbar
        colorbar, $
                position=pos, $
                range=[float(minfield),float(maxfield)], $
                divisions=ndiv, $
                ncolors=255,$
;               charsize=0.85,$
                ticklen=-0.15,$
                format=format_colorbar

	if (n_elements(nam2) ne 0) then begin
		thedate=STRSPLIT(nam2, /EXTRACT) 
		thedate=STRSPLIT(thedate(2),"'",/EXTRACT) 
		thedate=STRSPLIT(thedate(0),"_",/EXTRACT) 
		thedate=STRSPLIT(thedate(0),"-",/EXTRACT) 
		xyouts, 0.95, pos(1)-0.08, 'LMD Martian Mesoscale Model - d'+thedate(2)+'/m'+thedate(1), charsize=0.7, alignment=1., /NORMAL 
	endif else begin
		thedate=''
	endelse

      endelse	
endif

skip_new:

;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;
;;; ajouter un point
;loadct, 0
;xyouts, 353.87-360, -1.88, '+', ALIGNMENT=0.5, CHARSIZE=2, CHARTHICK=2
;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;goto, nocontour 
;print, 'contour' 
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

;lev=round(lev)
;yeah=UNIQ(lev) & lev=lev[yeah]

colcont = 150 ;75
if (pal eq 16) then colcont=0 ;noir
if (pal eq 4) then colcont=255 ;blanc

if (isotropic eq 'true') then begin
loadct,0,/silent
;reverse_ct
contour, overcontour, $
	lon,lat, $
/overplot, $ 
;	nlevels=colors/2, $
	levels=lev, $
	C_LINESTYLE = (lev LT 0.0), $   ;; dotted lines if negative
	max_value=max_contour, $
	min_value=min_contour, $	
	c_labels=findgen(n_elements(lev))*0., $
;        c_charsize=!P.charsize/2, $
	color=colcont;75;0;50;75;100tropclair;150;0;, $
;	/noerase, $
;	xtitle=xtitleset, $
;	xrange=windowx, $
;	ytitle=ytitleset, $
;	yrange=windowy, $
;       position=pospos, $
;	/isotropic, $
;	xtickinterval=intervalx, $
;	ytickinterval=intervaly 
endif else begin
loadct,0,/silent
contour, overcontour, $
        lon,lat, $
/overplot, $
;        nlevels=colors/2, $ 
        levels=lev, $
        C_LINESTYLE = (lev LT 0.0), $   ;; dotted lines if negative
        max_value=max_contour, $
        min_value=min_contour, $
        c_labels=findgen(n_elements(lev))*0., $
        color=colcont;, $
;        /noerase, $
;        xtitle=xtitleset, $ 
;        xrange=windowx, $
;        ytitle=ytitleset, $
;        yrange=windowy, $
;;       position=[0.15, 0.15, 0.95, 0.75], $
;        xtickinterval=intervalx, $
;        ytickinterval=intervaly 
endelse	


endif
nocontour:

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
if (n_elements(windex) eq 0) then windex=20. 
ref_mag=round(windex)
mytext=STRTRIM(STRING(ref_mag),2)+'ms!E-1!N'

if (n_elements(stride) eq 0) then begin
;stride=3 		; STRIDE: Pick every nth vector for plotting. (1)
;stride=2
stride=5
endif else begin
stride=round(stride)
endelse
	
ref_pos=[0.10, 0.90] 	
ref_pos=[0.15, 0.10]
ref_pos=[0.20, 0.05]
;ref_pos=[0.20, 0.02]

s=size(lon)
if (s[0] eq 1) then begin
 zeu=overvector_x
 zev=overvector_y
 longs=lon
 lats=lat
endif else begin
 minlat=min(lat) & maxlat=max(lat) & minlon=min(lon) & maxlon=max(lon)
 npoints=2*n_elements(lon(*,0))  ;; trop de points, mais au moins on ne perd rien (2 et 3 donnent results similaires)
 TRIANGULATE, lon, lat, tr
 zeu = GRIDDATA( lon, lat, overvector_x, /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
 zev = GRIDDATA( lon, lat, overvector_y, /LINEAR, triangles=tr, dimension=npoints, MISSING=!VALUES.F_NAN )
 longs =  minlon + (maxlon - minlon)*findgen(npoints)/float(npoints-1)
 lats =  minlat + (maxlat - minlat)*findgen(npoints)/float(npoints-1)
endelse

wlon=where( (longs gt windowx(0)) and (longs lt windowx(1)) )
wlat=where( (lats gt windowy(0)) and (lats lt windowy(1)) )

if ((wlon(0) eq -1) or (wlat(0) eq -1)) then begin
	print, 'wrong window settings'
        print, windowx(0), windowx(1)
	print, longs
	print, windowy(0), windowy(1)
	print, lats
	stop 
endif

myu=dblarr(n_elements(longs),n_elements(lats))
myv=dblarr(n_elements(longs),n_elements(lats))

myu(*,*)=!Values.F_NAN
myv(*,*)=!Values.F_NAN

for i=0,n_elements(wlon)-1 do begin
for j=0,n_elements(wlat)-1 do begin	
	myu(wlon(i),wlat(j))=zeu(wlon(i),wlat(j))
        myv(wlon(i),wlat(j))=zev(wlon(i),wlat(j))
endfor	
endfor	

	;; eliminate vector for colorbar clarity ?
	;nxs=1
	;nxe=n_elements(longs(*))-2
	;nys=1
	;nye=n_elements(lats(*))-2
	;myu=overvector_x(nxs:nxe,nys:nye)
	;myv=overvector_y(nxs:nxe,nys:nye)
	;longs=lon(nxs:nxe)
	;lats=lat(nys:nye)

loadct,0,/silent
VECTOR, myu, myv, longs, lats, LENGTH=len, REF_MAG=ref_mag,$
        STRIDE=stride, TYPE=2, HEAD_LEN=0.3, ANGLE=30, $
        REF_POS=ref_pos, REF_TEXT=mytext,$
        ALIGN=0.5, COLOR=1, /OVERPLOT;, /NOERASE

yeah:
endif


endelse



goto, why_contour_end
;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;

if (n_elements(overcontour) ne 1) then begin

        w=where(abs(overcontour) lt missing_value)
        if (w(0) ne -1) then begin
                min_contour=min(overcontour[w])
                max_contour=max(overcontour[w])
        endif


; to avoid spurious blank zones
nlevels=colors/2
if (n_elements(lev) eq 0) then lev=min_contour+(findgen(nlevels)/float(nlevels-1))*(max_contour-min_contour)
;lev=round(lev)
;yeah=UNIQ(lev) & lev=lev[yeah]


if (isotropic eq 'true') then begin
loadct,0,/silent
contour, overcontour, $
        lon,lat, $
;       nlevels=colors/2, $
        levels=lev, $
        C_LINESTYLE = (lev LT 0.0), $   ;; dotted lines if negative
        max_value=max_contour, $
        min_value=min_contour, $
        c_labels=findgen(nlevels), $ ;*0., $
        color=0, $
        /noerase, $
        xtitle=xtitleset, $
        xrange=windowx, $
        ytitle=ytitleset, $
        yrange=windowy, $
;        position=[0.15, 0.15, 0.95, 0.75], $
        /isotropic, $
        xtickinterval=intervalx, $
        ytickinterval=intervaly
endif else begin
loadct,0,/silent
contour, overcontour, $
        lon,lat, $
;       nlevels=colors/2, $ 
        levels=lev, $
        C_LINESTYLE = (lev LT 0.0), $   ;; dotted lines if negative
        max_value=max_contour, $
        min_value=min_contour, $
        c_labels=findgen(nlevels), $
        color=0, $
        /noerase, $
        xtitle=xtitleset, $
        xrange=windowx, $
        ytitle=ytitleset, $
        yrange=windowy, $
        position=[0.15, 0.15, 0.95, 0.75], $
        xtickinterval=intervalx, $
        ytickinterval=intervaly
endelse


endif

why_contour_end:

end
