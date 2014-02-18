pro out_geo_test, $
	field1=charvar, $
	field2=charvar2, $
	domain=domain, $
	path=path

;---------------------------------------------------------------;
;                                                               ;
; Use:                                                          ;
;       out_geo, 'HGT_M'	(or) out_geo, 'TOPO'		;
;	out_geo, 'THERMAL_INERTIA' (or) out_geo, 'TI'           ;
;       out_geo, 'ALBEDO_GCM' (or) out_geo, 'ALBEDO'            ;
;       out_geo, 'GHT'		                                ;
;								;
;	-- nested domain					;
;	out_geo, 'HGT_M', domain='d02'				;
;                                                               ;
; Options:							;
;		- domain: default is 'd01', 			;
;			but can be 'd02', 'd03', etc ...	;
;		- path: where the geo_em files are		; 
;			'/my_root/my_folder/'			;
;                                                               ;
; A. Spiga, April 2007 - July 2007 - September 2007             ;
;                                                               ;
;---------------------------------------------------------------;




;----------------------
; set parameters
;----------------------

colors=32
format='(F6.2)'
minfield_init=0.
maxfield_init=0.

default_path='/donnees/aslmd/big_data/WRF/WPS'
default_path='./'


if (n_elements(charvar) eq 0) then charvar='HGT_M'
if (n_elements(path) eq 0) then path=default_path
if (n_elements(domain) eq 0) then domain='d01'
file=path+'/geo_em.'+domain+'.nc'



;----------------------
; equivalent inputs
;----------------------
case charvar of
	'HGT_M':
	'TOPO': charvar='HGT_M'
	'THERMAL_INERTIA':
	'TI': charvar='THERMAL_INERTIA'
	'ALBEDO_GCM':
	'ALBEDO': charvar='ALBEDO_GCM'
	'GHT':
	else:
endcase	
if (n_elements(charvar2) ne 0) then begin
case charvar2 of
        'HGT_M':
        'TOPO': charvar2='HGT_M'
        'THERMAL_INERTIA':
        'TI': charvar2='THERMAL_INERTIA'
        'ALBEDO_GCM':
        'ALBEDO': charvar2='ALBEDO_GCM'
        'GHT':
	else:
endcase
endif							


;----------------------
; set graphics
;----------------------

set_plot, 'ps'
device, file='geo_em.'+domain+'_'+charvar+'.ps', /color
;set_plot, 'x'

case charvar of
        'HGT_M': pal=33 ;16  ;4
        'THERMAL_INERTIA': pal=3
        'ALBEDO_GCM': pal=0
        'GHT': pal=33
	else: pal=33
endcase
				


!p.charthick = 2.0
!p.thick = 3.0
!x.thick = 2.0
!y.thick = 2.0

case charvar of
	'HGT_M': title='Altimetry (km)'
	'THERMAL_INERTIA': title='Thermal inertia (!NJ!N.m!U-2!N.s!U-1/2!N.K!U-1!N)'
	'ALBEDO_GCM': title='Albedo LW (%)'
	'GHT': title='Geopotential height (m)'
	else: title=''
endcase


;-------------------------------
; open file and read variables
;-------------------------------
cdfid = ncdf_open(file)

; field
print, charvar
varid=ncdf_varid(cdfid,charvar)
        ncdf_varget, cdfid, varid, var
; lon
varid=ncdf_varid(cdfid,'XLONG_M')
        ncdf_varget, cdfid, varid, lon
	lon=reform(lon[*,0])
; lat
varid=ncdf_varid(cdfid,'XLAT_M')
        ncdf_varget, cdfid, varid, lat
	lat=reform(lat[0,*])

if (n_elements(charvar2) ne 0) then begin
	print, charvar2
	varid=ncdf_varid(cdfid,charvar2)
	ncdf_varget, cdfid, varid, overcontour
endif else begin
	overcontour=0.
endelse	
	

;-------------------------------
; calculation on variables
;-------------------------------
case charvar of
        'HGT_M': var=var/1000.
        'THERMAL_INERTIA': 
        'ALBEDO_GCM': var=var*100.
        'GHT': 
	else:
endcase

!p.multi=[0,3,3]
set_ls_tab=[90,270]
set_lct_tab=[09,12,16]
set_tau_tab=[0.3,1.5]

for yeahyeahyeah=0,1 do begin
for yeah=0,1 do begin
for yeahyeah=0,2 do begin
set_ls=set_ls_tab(yeah) & print, set_ls
set_lct=set_lct_tab(yeahyeah) & print, set_lct 
set_tau=set_tau_tab(yeahyeahyeah) & print, set_tau
if ((set_ls eq 90) and (set_tau eq 1.5)) then begin
		goto, beurk
endif
@slope.inc	

;---------
; plot !
;---------

minfield=min(var)
maxfield=max(var)
print, 'mM', minfield, maxfield

yysize=150
what_I_plot=var
lon_=lon
lat_=lat
;what_I_plot=congrid(reform(var),yysize,yysize,/INTERP)
;lon_=congrid(lon,yysize)
;lat_=congrid(lat,yysize)

map_latlon, $
        what_I_plot, $                          ; 2D field
        lon_, $                                  ; 1D latitude
        lat_, $                                  ; 1D longitude
        minfield=minfield_init, $               ; minimum value of plotted field (=0: calculate)
        maxfield=maxfield_init, $               ; maximum value of plotted field (=0: calculate)
        overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
;        overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
;        overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
        ct=pal, $                               ; color table (33-rainbow is default)
        colors=colors, $                        ; number of colors/levels (32 is default)
        title=title, $                          ; title of the plot ('' is default)
        format=format                           ; format of colorbar annotations ('(F6.2)' is default)
	
;;;;;;;;;;;;;;;;;
beurk:
endfor
endfor
endfor
;;;;;;;;;;;;;;;;;


device, /close

end
