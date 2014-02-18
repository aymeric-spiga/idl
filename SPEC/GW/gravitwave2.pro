pro gravitwave2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;; KEEP REFERENCES of FIRST PLOTS 
	;file='../EN_COURS/gw/ok_hires_wrfout_d01_9999-09-09_09:00:00_z'
	;file='../EN_COURS/gw/_wrfout_d01_9999-09-09_09:00:00_z'
	;file='../LMD_MM_MARS/TESTGW/_wrfout_d01_9999-09-09_09:00:00_z'
	;;file='../LMD_MM_MARS/TESTGW/wind15_wrfout_d01_9999-09-09_09:00:00_z'
	;file='../TESTGW/wind20c_wrfout_d01_9999-09-09_09:00:00_z'
folder='/d5/aslmd/GRAVITWAVE/GW_MARS_highwind_3D.157077/'
folder='/d5/aslmd/GRAVITWAVE/GW_MARS_highwind_narrowmountain_2D.156878/'
;folder='/d5/aslmd/GRAVITWAVE/GW_MARS_highwind_narrowmountain_3D.157025/'
folder='/d5/aslmd/GRAVITWAVE/GW_MARS_highwind_narrowmountain_2D.morepoints.160818/'
folder='/d5/aslmd/GRAVITWAVE/GW_MARS_highwind_3D_loweropacity_morepoints_LT15.161544/'
folder='/home/aslmd/GRAVITWAVE/GW_MARS_highwind_3D_loweropacity_morepoints.161542/'
;folder='/home/aslmd/GRAVITWAVE/GW_MARS_highwind_3D_loweropacity_morepoints_LTcst.168440/'
;folder='/home/aslmd/GRAVITWAVE/GW_MARS_highwind_3D_loweropacity_morepoints_widehill.168468/'
file=folder+'wrfout_d01_9999-09-09_09:00:00_z'
charvar='W' & charvarc='W'
charvar='W' & charvarc='tk' & cond=0
;charvar='tk' & charvarc='W'
;charvar='tk' & charvarc='W' & cond=1
;charvar='tk' & charvarc='tk' & cond=1
;charvar='tk' & charvarc='tk' & cond=0
;;charvar='tpot' & charvarc='tpot'
charvar='HR_NLTE' & charvarc='tk' & cond=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;folder='/home/aslmd/GRAVITWAVE/'
	;file=folder+'GW_MARS_highwind_3D_loweropacity_morepoints_diff.nc'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
what_I_plot=0. & overcontour=0.
SPAWN, '\rm param_plot.idl ; cp gravitwave_inc.pro param_plot.idl'
;
;
;
getcdf, $
        file=file, $
        charvar=charvar, $
        invar=invar
getcdf, $
        file=file, $
        charvar=charvarc, $
        invar=invarc
;getcdf, $
;        file=file, $
;        charvar='vert', $
;        invar=vert
getcdf, $
         file='/home/aslmd/GRAVITWAVE/GW_MARS_highwind_3D.157077/wrfout_d01_9999-09-09_09:00:00_z', $
         charvar='vert', $
         invar=vert
getcdf, $
        file=file, $
        charvar='PTOT', $
        invar=columnp
;
;
;
;for nt=4,10 do begin
;for nt=27,27 do begin
for nt=1,40 do begin
;for nt=9,9 do begin
;for nt=0,20 do begin
zefile='gravitwave_'+charvar+'_'+charvarc+'_'+string(nt+100,'(I0)')
if (n_elements(cond) ne 0) then if (cond ne 0) then zefile = 'cond_' + zefile
;
;
;
  PS_Start, filename=zefile+'.ps' & print, zefile+'.ps'
;
;
;
s = size(invar) & middle = floor(s[2]/2) & print, 'PLOT at y subs ', middle
what_I_plot = reform(invar(*,middle,*,nt)) & print, min(what_I_plot), max(what_I_plot)
overcontour = reform(invarc(*,middle,*,nt))
   if (n_elements(cond) ne 0) then begin
   if (cond ne 0) then begin
     column = columnp
     ;column = reform(columnp(*,1,*,1))
     @tempcond.inc
     yeye = reform(overplot(*,middle,*,nt))
     w = where(yeye le 0.) & yeye[w] = 0. 
     if (charvar eq 'tk') then what_I_plot = what_I_plot - yeye 
     if (charvarc eq 'tk') then overcontour = overcontour - yeye 
   endif
   endif
xx = findgen(n_elements(what_I_plot(*,0)))
zz = vert / 1000.
;
	;;w = where(what_I_plot lt 2.)
	;;w = where(what_I_plot lt 0.)
	;w = where(what_I_plot lt -1.)
	;;w = where(what_I_plot lt -3.)
	;if (w(0) ne -1) then begin
	;for i=0,n_elements(w)-1 do begin
	; ind = ARRAY_INDICES(what_I_plot, w(i))
	; print, 'x '+string(ind(0),'(I0)')+' z '+string(zz(ind(1)),'(I0)')+' t '+string(nt*925./3700.,'(F4.1)'), what_I_plot(ind[0],ind[1])
	;endfor
	;endif
;
;
;
if (charvar eq 'HR_NLTE') then begin
	what_I_plot = what_I_plot * 3700. ;; en K/hour
	;w = where((zz le 60.) and (zz ge 50.)) ;; pour enlever la ligne d'altitude ou le NLTE est active
        ;if (w(0) ne -1) then what_I_plot(*,w) = 0.
endif
;
;
;
map_latlon, $
        what_I_plot, $                          ; 2D field
        xx, $                                   ; 1D latitude
        zz, $                                   ; 1D longitude
;        minfield_init, $               ; minimum value of plotted field (=0: calculate)
;        maxfield_init, $               ; maximum value of plotted field (=0: calculate)
        overcontour=overcontour, $              ; another 2D field to overplot with contour lines (=0: no)
;        overvector_x=overvector_x, $            ; wind vector - x component (=0: no)
;        overvector_y=overvector_y, $            ; wind vector - y component (=0: no)
;        ct=pal, $                               ; color table (33-rainbow is default)
;        colors=colors, $                        ; number of colors/levels (32 is default)
;        title=title_user, $                     ; title of the plot ('' is default)
        format=format                           ; format of colorbar annotations ('(F6.2)' is default)
;
;
;
PS_End, /PNG
SPAWN, 'mv '+zefile+'.png '+folder+' ; \rm '+zefile+'.*ps*'
endfor
end
