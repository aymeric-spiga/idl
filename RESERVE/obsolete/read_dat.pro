pro read_dat, output, x, y, z, t, desc_field, nodata=nodata


;---------------------------------------------------------------------------
;
; Read the files generated by ARWpost
;
; Please create symbolic link 'input.ctl' and 'input.dat'
; --- or rename files
;
; - output contains all fields
; - x,y,z,t contains coordinates
; - desc_field contains strings describing fields
; - if nodata is not blank, do not read the fields
;
;---------------------------------------------------------------------------
; A. Spiga, August 2007
;---------------------------------------------------------------------------


;
; info for localtime
;
interv=0.
openr, 99, 'timefil'
readf,99,interv
close, 99



; get preliminar info on coordinates
;
wnx=1
wny=1
wnz=1
wnt=1
txt=''
param=[0.,0.,0.,0.,0.,0.,0.,0.]
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

nx=round(dimensions(0)) 
ny=round(dimensions(1)) 
nz=round(dimensions(2)) 
nt=round(dimensions(3)) 
fields=round(dimensions(4)) 
vertical=intarr(fields)
for i=0,fields-1 do begin
	vertical(i)=round(dimensions(5+i))
endfor


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
	

;
; get info on fields
;
openr, 1, 'input.ctl'
toto=5+wnx+1+wny+1+wnz+2+fields
infofields=strarr(toto)
readf, 1, infofields
close, 1

;
desc_field=transpose(infofields(toto-fields:toto-1))


if (n_elements(nodata) eq 0) then begin
;---------------------------------
; READ .DAT FILE DATA
;---------------------------------
; GrADS views gridded data sets as 5-dimensional arrays 
; varying in longitude, latitude, vertical level, variable, and time.
;---------------------------------

print, '----------------'
print, 'reading data ...'
print, '----------------'
print, 'longitudes: ', min(x), max(x), nx, ' points'
print, 'latitudes: ', min(y), max(y), ny, ' points' 
print, 'altitudes: ', min(z), max(z), nz, ' points'
print, 'local times: ', min(t), max(t), nt, ' points'
print, fields, ' recorded fields :'
print, desc_field
print, '----------------'

	;;very slow method
	;;----------------	
	;data=read_ascii('input.dat', missing_value=1e37)
	;datadat=data.field1

; see http://www.dfanning.com/tips/ascii_column_data.html	
OPENR, lun, 'input.dat', /GET_LUN

;count=0L
output=fltarr(fields,nx,ny,nz,nt)
dummy=0.

for l=0,nt-1 do begin 						;time loop
print, 'time loop ...', l+1, ' / ', string(nt,'(I0)')  	
	for f=0, fields-1 do begin ;fields (whether 2D or 3D)	;variable loop
	nzf=vertical(f)
	print, 'reading field ...', f+1, ' / ', string(fields,'(I0)') 
		for k=0,nzf-1 do begin  			;vertical loop
			for j=0,ny-1 do begin			;latitude loop
			for i=0,nx-1 do begin			;longitude loop

				;;very slow method
				;output(f,i,j,k,l)=datadat(count)
                                ;count=count+1

				READF, lun, dummy
				output(f,i,j,k,l)=dummy

			endfor
			endfor
		endfor
	endfor
endfor

; missing values
w=where(output eq 1e20)
if (w(0) ne -1) then begin
	print, 'change missing values ...'
	output[w]=1e37
endif


CLOSE, lun
;SPAWN, '\rm input.dat_tmp'
print, '...done !'
help, /memory
print, '---------'

;if (count ne n_elements(datadat)) then begin
;	print, 'some data were not read !'
;	stop
;endif

;;
;; Manage missing data so as IDL recognize them
;;
;w=where(output eq 1e37)
;if (w(0) ne -1) then output[w]=!Values.F_NAN

endif

end
