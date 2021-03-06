pro gw, typeplot


;-----------------------------------------
; Read files generated by ARWpost & plot
;
; typeplot is s for see		(.ps are not saved, .def saved)
;	      t for trace	(.ps are saved, .def saved)
;	      r for report	(same as t + tex report)
;	      m for movie	(same as t + gif movie)
;	      w for web page    (same as t + html animation)
;             g for simple grads output
;
; default is 's'
;
; use:
;  1. set parameters in gw.def
;  2. @gwc
;  3. gw, 's' (or 't' or 'r' or 'm' or 'w' or 'g')
;     or simply: gw
; 
;-----------------------------------------
; A. Spiga, Mars 2008
;-----------------------------------------


nomovie='yes'
grads='blabla'
useidlfile='false'

if (n_elements(typeplot) eq 0) then typeplot='s'
if (typeplot eq 'g') then grads='true'

print, '---------------------'
print, 'init. please wait ...'
print, '---------------------'

;
; init.
;
field1='' & field2='' & winds=['','']
backup_data='no' & already_data='no'
datafolder='./' & plotfolder='./'
textitle='plot' & texcomments='plot'
topo=0
extract='yes' ;;a question is asked

;
; defining fields, coordinates, type of plot ... etc ...
;
@user.idl

;
; extracting/interpolating fields from WRF NETCDF files
;
if (typeplot eq 'g') then SPAWN, 'cp -f gw.def '+plotfolder+'/grads_gw.def'
if (useidlfile eq 'true') then extract='no'
if (extract eq 'yes') then begin
	print, 'EXTRACT DATA. OK ?'
        if (n_elements(dumb) eq 0) then begin  ;; give dumb in user.idl if one does not want to be asked questions !
            if (grads ne 'true') then read, dumb, prompt='type 1 to proceed, 0 to keep the same input >>> ' else dumb=999
        endif
	if (dumb ne 0) then begin
        	print, 'WE WILL EXTRACT DATA USING ARWPOST.EXE'
		if (n_elements(nest) eq 0) then nest='d01'
        	call_arwpost, nam1, nam2, nam3, nam4, nam5, tabnam, nest, datafolder, grads
	endif 
endif
if (already_data eq 'yes') then begin
	if (nam5 eq " interp_method = 0 ") then coord='model_level'
	if ((nam5 eq " interp_method = 1 ") and (tabnam(0) lt tabnam(1))) then coord='height'
        SPAWN, 'ln -sf  '+datafolder+'/'+coord+'.input.ctl input.ctl'
        SPAWN, 'ln -sf  '+datafolder+'/'+coord+'.input.dat input.dat'
endif

;-----------------------------;
; STARTING TO READ INPUT.CTL  ;
;-----------------------------;

;
; get the local time
;
interv=0. & openr, 99, 'timefil' & readf,99,interv & close, 99

;
; get preliminar info on coordinates
;
wnx=1 & wny=1 & wnz=1 & wnt=1 & txt='' & param=[0.,0.,0.,0.,0.,0.,0.,0.]

openr, 1, 'input.ctl'
while not eof(1) do begin
    readf,1,txt
    match=STRPOS(txt,'linear')
	if (match ne -1) then begin
		wnx=(STRPOS(txt,'xdef') eq -1)*wnx  ;0 if linear
		wny=(STRPOS(txt,'ydef') eq -1)*wny  ;0 if linear
                wnz=(STRPOS(txt,'zdef') eq -1)*wnz  ;0 if linear 
                wnt=(STRPOS(txt,'tdef') eq -1)*wnt  ;0 if linear
                paramstart=STRSPLIT(txt, 'r', LENGTH=paramlength)
        	paramchar=STRMID(txt,paramstart(1),paramlength(1))
		READS, paramchar, start, step
		param=param+[start*(STRPOS(txt,'xdef') ne -1),$
			step*(STRPOS(txt,'xdef') ne -1),$
			start*(STRPOS(txt,'ydef') ne -1),$
			step*(STRPOS(txt,'ydef') ne -1),$
			start*(STRPOS(txt,'zdef') ne -1),$
			step*(STRPOS(txt,'zdef') ne -1),$
			start*(STRPOS(txt,'tdef') ne -1),$
			step*(STRPOS(txt,'tdef') ne -1)]
	endif	
endwhile
close, 1
timebegin=param(6)
	

;---------------------------------
; READ .CTL FILE INFORMATIONS
;---------------------------------

info=read_ascii('input.ctl', missing_value=1e37)
infodat=info.field1

;
; second column : dimensions
;

w=where((infodat(0,*) eq 1e37) and (infodat(1,*) lt 1e37))
dimensions=reform(infodat(1,w))

nx=round(dimensions(0)) & ny=round(dimensions(1)) 
nz=round(dimensions(2)) & nt=round(dimensions(3)) 
fields=round(dimensions(4)) & vertical=intarr(fields)
for i=0,fields-1 do vertical(i)=round(dimensions(5+i))

;
; first column : coordinates
;

w=where((infodat(1,*) eq 1e37) and (infodat(0,*) lt 1e37))
if (w(0) ne -1) then coordinates=reform(infodat(0,w))

if (wnx ne 0) then begin

	wnx=nx
	x=coordinates(0:nx-1)

endif else begin

	xmin=param(0)
	xmax=param(0)+param(1)*(nx-1) MOD 360
	;inc=param(0)+param(1)*(nx-1)	;cas mercator
	;x=start+inc*findgen(nx) 
	x=xmin+(xmax-xmin)*findgen(nx)/(nx-1)
	
;param(1)=param(1)-param(0) ;;bricolo plot GCM
;xmin=param(0)
;xmax=param(0)+param(1)*(nx-1)
;x=xmin+(xmax-xmin)*findgen(nx)/(nx-1)
;;x = x - (LONG(x)/180)*360.0
;print, x

	;stop

endelse


if (wny ne 0) then begin
        wny=ny
	y=coordinates(wnx:wnx+wny-1)
endif else begin
        start=param(2)
        step=param(3)
        y=start+step*findgen(ny)
endelse
if (wnz ne 0) then begin
        wnz=nz
	z=coordinates(wnx+wny:wnx+wny+wnz-1)
endif else begin
        start=param(4)
        step=param(5)
        z=start+step*findgen(nz)
endelse
t=interv*findgen(nt)/3700. + timebegin + mean(x)/15.
lct = round(t+24) MOD 24	
	
;
; get info on fields
;
openr, 1, 'input.ctl'
toto=5+wnx+1+wny+1+wnz+2+fields
infofields=strarr(toto)
readf, 1, infofields
close, 1
desc_field=transpose(infofields(toto-fields:toto-1))

;
; vert coord
;
if (z(0) lt z(1)) then coord='height'
if (z(0) gt z(1)) then coord='pressure'
if ((z(0) eq 1.) and (z(nz-1) eq (z(nz-2)+1.))) then coord='model_level'



;---------------------------------
; READ .DAT FILE DATA
;---------------------------------
; GrADS views gridded data sets as 5-dimensional arrays 
; varying in longitude, latitude, vertical level, variable, and time.
;---------------------------------

;-----------------------------;
; STARTING TO READ INPUT.DATA ;
;-----------------------------;
print, '----------------'
print, 'reading data ...'
print, '----------------'
print, 'longitudes: ', min(x), max(x), nx, ' points'
print, 'latitudes: ', min(y), max(y), ny, ' points' 
print, 'altitudes: ', min(z), max(z), nz, ' points'
print, 'local times: ', min(lct), max(lct), nt, ' points'
print, fields, ' recorded fields :'
print, desc_field
print, '----------------'

; see http://www.dfanning.com/tips/ascii_column_data.html	
OPENR, lun, 'input.dat', /GET_LUN


;
; titles and plot name
;
denom=textitle+'_'+plot+'_'+coord+'_'+field1
if ((field2 ne '') and (topo eq 0)) then denom=denom+'_'+field2
if (topo eq 1) then denom=denom+'_HGT'
if (winds(0) ne '') then denom=denom+'_'+winds(0)+winds(1)
        case coord of
	'height': charlev='Altitude (km)'
	;'pressure': begin
	;               H=10.
	;               p0=6.1
	;               z=H*alog(p0/z)
	;               charlev='Log-Pressure altitude (km)'
	;end
	'model_level': charlev='Model level'
	endcase
																							

;
; graphics settings
;
!p.charthick = 2.0
!p.thick = 3.0
!x.thick = 2.0
!y.thick = 2.0

	; plot with topography
	if (topo eq 1) then field2='HGT'

save_ps='true' ;; hardwired
;save_ps='false'
if (save_ps eq 'false') then begin
		;if (save_ps ne 'false') then begin
		;	;set_plot, 'ps'
		;endif else begin
	set_plot, 'x'
	!P.MULTI = [0, 2, 3]
	window, 0, xsize=600, ysize=900
		;endelse
endif

no3D='false'	
;goto, yeye 
print, 'GENERATE PLOTS ',min([nt-1,num-1])+1

for l=0,min([nt-1,num-1]) do begin 	;time loop
		print, '************************************************************************'
		print, 'time loop ...', l+1, ' / ', string(nt,'(I0)')  
		T = SYSTIME(1)
	if (useidlfile eq 'true') then $
		restore,filename=plotfolder+'/'+denom+string(1000+l,'(I0)')+'.idl' $
	else begin
  for f=0, fields-1 do begin ;fields (whether 2D or 3D)	;variable loop
  nzf=vertical(f)
  print, 'scanning field ...', f+1, ' / ', string(fields,'(I0)')

	data = FLTARR(1,nzf*ny*nx)
	print, 'read data ...'
	READF, lun, data	;; the IDL way, two times faster than 3 FOR loops
	print, 'read data done'
	
	recognize=STRSPLIT(desc_field(f), /EXTRACT)
	recognize=recognize[0]
	
;	IF (STRPOS(desc_field(f),field1) ne -1) then begin
;       IF ( ( STRPOS(desc_field(f),field1) eq 0 ) and ( STRLEN(desc_field(f)) eq STRLEN(field1) ) ) then begin
	IF (recognize eq field1) then begin
		data_save=data
		print, 'found field1 ', field1
		data=REFORM(data,nx,ny,nzf,/OVERWRITE)
			w=where(data eq -9999.)	; missing values: -9999 (before, was 1e20)
			if (w(0) ne -1) then data[w]=1e37
		CASE (plot) OF
		'map': what_I_plot=reform(data(*,*,-1+min([level,nzf])))
		'meridional': what_I_plot=reform(data(-1+min([nlon,nx]),*,*))
		'zonal': what_I_plot=reform(data(*,-1+min([nlat,ny]),*))
		ENDCASE
		if (nzf eq 1) then no3D='true'
		data=0. & data=data_save & data_save=0.
	ENDIF
;	IF ((STRPOS(desc_field(f),field2) ne -1) and (field2 ne '')) then begin
;       IF ((STRPOS(desc_field(f),field2) eq 0) and (field2 ne '')) then begin
;       IF ( ( STRPOS(desc_field(f),field2) eq 0 ) and ( field2 ne '' ) and ( STRLEN(desc_field(f)) eq STRLEN(field2) ) ) then begin
        IF (recognize eq field2) then begin
		data_save=data
		print, 'found field2 ', field2
		if (topo eq 1) then begin
			data=data/100.
			data=round(data)
			data=float(data)/10.
		endif	
		data=REFORM(data,nx,ny,nzf,/OVERWRITE)
                        w=where(data eq -9999.)   ; missing values
			if (w(0) ne -1) then data[w]=1e37
		CASE (plot) OF
		'map': overcontour=reform(data(*,*,-1+min([level,nzf])))
		'meridional': overcontour=reform(data(-1+min([nlon,nx]),*,*))
		'zonal': overcontour=reform(data(*,-1+min([nlat,ny]),*))
		ENDCASE
		data=0. & data=data_save & data_save=0.
	ENDIF   
;	IF (STRPOS(desc_field(f),winds(0)) ne -1 and (winds(0) ne '')) then begin
;       IF (STRPOS(desc_field(f),winds(0)) eq 0 and (winds(0) ne '')) then begin
	IF (recognize eq winds(0)) then begin
		data_save=data
		print, 'found wind1 ', winds(0)
	        data=REFORM(data,nx,ny,nzf,/OVERWRITE)
                        w=where(data eq -9999.)   ; missing values
			if (w(0) ne -1) then data[w]=!Values.F_NAN
		CASE (plot) OF
                'map': overvector_x=reform(data(*,*,-1+min([level,nzf])))
                'meridional': overvector_x=reform(data(-1+min([nlon,nx]),*,*))
                'zonal': overvector_x=reform(data(*,-1+min([nlat,ny]),*))
		ENDCASE
		data=0. & data=data_save & data_save=0.
	ENDIF 
;	IF (STRPOS(desc_field(f),winds(1)) ne -1 and (winds(1) ne '')) then begin
;       IF (STRPOS(desc_field(f),winds(1)) eq 0 and (winds(1) ne '')) then begin
        IF (recognize eq winds(1)) then begin
		data_save=data
		print, 'found wind2 ', winds(1)
		data=REFORM(data,nx,ny,nzf,/OVERWRITE)
                        w=where(data eq -9999.)   ; missing values
	                if (w(0) ne -1) then data[w]=!Values.F_NAN
                CASE (plot) OF
                'map': overvector_y=reform(data(*,*,-1+min([level,nzf])))
                'meridional': overvector_y=reform(data(-1+min([nlon,nx]),*,*))
                'zonal': overvector_y=reform(data(*,-1+min([nlat,ny]),*))
                ENDCASE
		data=0. & data=data_save & data_save=0.
	ENDIF
	data=0.		;; free memory
	data_save=0.         ;; free memory
	
  endfor
	  save, $
	 	what_I_plot, $
	        overcontour, overvector_x, overvector_y, $
	 	x, y, z, lct, filename=plotfolder+'/'+denom+string(1000+l,'(I0)')+'.idl'
	  endelse

	if (n_elements(overcontour) eq 0) then overcontour=0.
	if (n_elements(overvector_x) eq 0) then overvector_x=0.
	if (n_elements(overvector_y) eq 0) then overvector_y=0.

	@gw_name.idl

;;--------------------------------------------------------
;; NB: you can use the winds fields to combine variables
;;--------------------------------------------------------
@wind_stress.idl	;** WIND STRESS
@sea_level.idl		;** HYDROSTATIC REDUCTION TO LEVEL OF REFERENCE
;@ertel_vorticity.idl	;** ERTEL POTENTIAL VORTICITY
@tau_cloud.idl          ;** CLOUD OPACITY
@wind_velocity.idl	;** 3D WIND VELOCITY
@hwind_velocity.idl	;** HORIZONTAL WIND VELOCITY 
;what_I_plot=what_I_plot-median(what_I_plot)
;;--------------------------------------------------------
;; NB: you can use the winds fields to combine variables
;;--------------------------------------------------------

	title=title+', LT '+string(lct(l),'(I0)')

;;--------------------------------------------------------
;;--------------------------------------------------------
;       getcdf, file='zefolder/wrfinput_d01', charvar='MARS_TI', invar=what_I_plot ;;; changer minfield et maxfield
       ;;print, size(what_I_plot)
       ;;getcdf, file='../WPS/geo_em.d01.nc', charvar='THERMAL_INERTIA', invar=yeyeye
       ;;title = 'Thermal Inertia (J m!U-2!N s!U-0.5!N K!U-1!N)' & no3D='true' & pal=3 & format='(I0)'
       ;;what_I_plot(2:153,2:153)=reform(yeyeye(*,*,0))
       ;what_I_plot=reform(what_I_plot(*,*,0)) & help, what_I_plot
;       getcdf, file='zefolder/wrfinput_d01', charvar='MARS_ALB', invar=what_I_plot ;;; changer minfield et maxfield
;       title = 'Albedo' & no3D='true' & pal=1 & format='(F5.2)'
;       what_I_plot=reform(what_I_plot(*,*,0)) & help, what_I_plot
;;--------------------------------------------------------
;;--------------------------------------------------------

	if (save_ps ne 'false') then begin
		save_ps=denom+string(1000+l,'(I0)')
			;        device, filename=save_ps+'.ps', /color, bits=16
		PS_Start, filename=save_ps+'.ps'
                !P.Charsize = 1.2
	endif

;**************************
smoothampl=2
smoothampl=0 ;; no smooth
;**************************

	CASE (plot) OF
	'map': begin
        	case coord of
		   'height': charlev=string(round(z(level-1)*1000.),'(I0)')+' m AMR'
		   'pressure': charlev=string(round(z(level-1)*100.),'(I0)')+' Pa'
		   'model_level': charlev=string(round(z(level-1)),'(I0)')
		endcase
        	if (no3D ne 'true') then title=title+', at level '+charlev
		print, title
	map_latlon, $
		smooth(what_I_plot,[smoothampl,smoothampl],/EDGE_TRUNCATE), $	; 2D field
        	x, $                                  ; 1D latitude
        	y, $                                  ; 1D longitude
        	minfield=0., $               ; minimum value of plotted field (=0: calculate)
        	maxfield=0., $               ; maximum value of plotted field (=0: calculate)
        	overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
        	overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
        	overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
        	ct=pal, $                               ; color table (33-rainbow is default)
        	colors=colors, $                        ; number of colors/levels (32 is default)
        	title=title, $                     ; title of the plot ('' is default)
        	format=format                           ; format of colorbar annotations ('(F6.2)' is default)
	end
	'zonal': begin
		displaylat=float(round(y(nlat-1)*10))/10.
		title=title+ ', section at latitude '+string(displaylat,format='(1F6.1)')+ '!Uo!N '
	        print, title
	section, $
		smooth(what_I_plot,[smoothampl,smoothampl],/EDGE_TRUNCATE), $   ; 2D field
                x, $                                ; horizontal coordinate
                z, $                             ; altitude coordinate
                minfield=0., $               ; minimum value of plotted field (=0: calculate)
                maxfield=0., $               ; maximum value of plotted field (=0: calculate)
                minspace=0., $                    ; minimum value of space window (=0: calculate)
                maxspace=0., $                    ; maximum value of space window (=0: calculate)
                overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
                overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
                overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
                colors=colors, $                        ; number of colors/levels (32 is default)
                title_plot=title, $                ; title of the plot ('Profile' is default)
                title_axis=['Longitude',charlev], $             ; title of the [x,y] axis (['Field','Altitude'] is default)
                ct=pal, $                               ; color table (33-rainbow is default)
                format=format
 	end
        'meridional': begin
                displaylon=float(round(x(nlon-1)*10))/10.
                title=title+', section at longitude '+string(displaylon,format='(1F6.1)')+ '!Uo!N '
                print, title
        section, $
                smooth(what_I_plot,[smoothampl,smoothampl],/EDGE_TRUNCATE), $   ; 2D field
                y, $                                ; horizontal coordinate
                z, $                             ; altitude coordinate
                minfield=0., $               ; minimum value of plotted field (=0: calculate)
                maxfield=0., $               ; maximum value of plotted field (=0: calculate)
                minspace=0., $                    ; minimum value of space window (=0: calculate)
                maxspace=0., $                    ; maximum value of space window (=0: calculate)
                overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
                overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
                overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
                colors=colors, $                        ; number of colors/levels (32 is default)
                title_plot=title, $                ; title of the plot ('Profile' is default)
                title_axis=['Latitude',charlev], $             ; title of the [x,y] axis (['Field','Altitude'] is default)
                ct=pal, $                               ; color table (33-rainbow is default)
                format=format
        end
 	endcase	 

PS_End, /PNG

;		 if (save_ps ne 'false') then begin
;	               device, /close
;	;               SPAWN, './ps2png -a ; \rm '+save_ps+'.ps '+save_ps+'.eps '
;	               SPAWN, './ps2png ; \rm '+save_ps+'.ps '+save_ps+'.eps '
	               if (typeplot eq 's') then begin
	                 SPAWN, 'display '+save_ps+'.png &'
	               endif 
;	         endif

;PRINT, 1000*(SYSTIME(1) - T), ' milliSeconds'
PRINT, SYSTIME(1) - T, ' seconds to process'
endfor

CLOSE, lun
print, '...done !'
help, /memory
print, '---------'


yeye:

CASE (typeplot) OF
;'g': begin
;        SPAWN, 'cp -f gw.def '+plotfolder+'/grads_'+denom+'.def' 
;end
's': begin
	;SPAWN, 'find '+denom+'????.ps -exec display {} \;'
	SPAWN, 'cp -f gw.def '+plotfolder+'/see_'+denom+'.def' ; \rm ./'+denom+'????.ps ./'+denom+'????.eps ' ;;./'+denom+'????.png'
end
't': begin
	SPAWN, 'mv -f ./'+denom+'????.ps ./'+denom+'????.png '+plotfolder+'/'
	SPAWN, 'cp -f gw.def '+plotfolder+'/'+denom+'.def'
end
'r': begin
	SPAWN, 'touch report.comments ; \rm report.comments ; echo '+textitle+' > report.comments ; echo '+texcomments+' >> report.comments'
	SPAWN, 'rm -f '+denom+'_.ps ; echo executing ... base '+denom+'????.ps < report.comments'
	SPAWN, 'base '+denom+'????.ps < report.comments ; mv -f base.ps '+denom+'_.ps'

	SPAWN, 'mv -f ./'+denom+'????.ps '+plotfolder+'/'
	SPAWN, 'mv -f ./'+denom+'_.ps '+plotfolder+'/'
	SPAWN, 'cp -f gw.def '+plotfolder+'/'+denom+'.def'
end
'm': begin
	print, 'generating movie ...'
	SPAWN, 'rm -rf temp ; mkdir temp ; mv -f '+denom+'????.ps temp ; cp bigconvert200 temp ; cd temp'
	SPAWN, 'bigconvert200 *.ps ; convert -delay 60 *.png movie.gif ; cd ..'
	if (nomovie ne 'yes') then SPAWN, 'animate temp/movie.gif &' 
	print, 'movie done'

	SPAWN, 'mv -f temp/'+denom+'????.ps '+plotfolder+'/'
	SPAWN, 'mv -f temp/movie.gif '+plotfolder+'/'+denom+'.gif'
	SPAWN, 'mv -f temp/'+denom+'????.png '+plotfolder+'/'
	SPAWN, 'cp -f gw.def '+plotfolder+'/'+denom+'.def'
end
'w': begin
;
;
netfolder='/u/aslmd/WWW/ANIMATIONS/'
;netfolder='/home/aslmd/'
;
; set also textitle
;
        ;print, 'generating frames ...'
	;;        SPAWN, 'rm -rf temp ; mkdir temp ; mv -f '+denom+'????.ps temp/ ; ls temp/'
	;;;        SPAWN, "find temp/* -exec convert -antialias -density 100x100 -crop 0x0 {} '{}.png' \;"
	;;SPAWN, "cd temp ; ../ps2png ; cd .."
	;;	SPAWN, 'ls temp/'
	;;        ;SPAWN, 'mv -f temp/'+denom+'????.ps '+plotfolder+'/'
        ;SPAWN, 'cp -f gw.def '+plotfolder+'/web_'+denom+'.def'
	;;; webpage
	;;SPAWN, 'rm -rf '+netfolder+textitle+'_'+denom+' ; mkdir '+netfolder+textitle+'_'+denom
	;SPAWN, 'mkdir '+netfolder+textitle
	;;SPAWN, 'mv -f temp/'+denom+'????.png '+netfolder+textitle+'/'
	;SPAWN, 'mv -f temp/'+denom+'????.png '+netfolder+textitle+'/'
	;SPAWN, 'chmod 777 -R '+netfolder+textitle
	;SPAWN, 'rm -rf '+denom+'.html ; more header.html > '+textitle+'_'+denom+'.html'
	;genanim, num, textitle, denom 
	;SPAWN, 'more templist >> '+textitle+'_'+denom+'.html ; rm -rf templist'
	;SPAWN, 'more body.html >> '+textitle+'_'+denom+'.html'
	;SPAWN, 'mv -f '+textitle+'_'+denom+'.html '+netfolder+'/'

   print, 'generating frames ...'
   cmd = 'echo unix'
   cmd = cmd + ' ; cp -f gw.def '+plotfolder+'/web_'+denom+'.def'
   ;;; webpage
   ;cmd = cmd + ' ; mkdir '+netfolder+textitle
   cmd = cmd + ' ; mv -f ./'+denom+'????.png '+netfolder+'/' ;+textitle+'/'
   ;cmd = cmd + ' ; chmod 777 -R '+netfolder+textitle
   cmd = cmd + ' ; rm -rf '+denom+'.html ; cat header.html > '+denom+'.html'
   SPAWN, cmd
   genanim, num, textitle, denom 
   cmd = 'echo unix'
   cmd = cmd + ' ; cat templist >> '+denom+'.html ; rm -rf templist'
   cmd = cmd + ' ; cat body.html >> '+denom+'.html'
   cmd = cmd + ' ; mv -f '+denom+'.html '+netfolder+'/'
   cmd = cmd + ' ; rm -rf ./*.epsi ./*.ps'
   SPAWN, cmd

end	
ENDCASE

if (backup_data eq 'yes') then SPAWN, 'cp -f namelist.ARWpost '+datafolder+'/'+coord+'.namelist.ARWpost'
if (backup_data eq 'yes') then SPAWN, 'cp -f input.ctl '+datafolder+'/'+coord+'.input.ctl'
if (backup_data eq 'yes') then SPAWN, 'cp -f input.dat '+datafolder+'/'+coord+'.input.dat'




end
