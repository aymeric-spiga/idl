pro call_arwpost, nam1, nam2, nam3, nam4, nam5, tabnam, nest, datafolder, grads

;
; YOU MUST HAVE COMPILED ARWPOST.EXE !
;

;
; USE:
;  @plot.def
;  then: call_arwpost, nam1, nam2, nam3, nam4, nam5, tabnam
;
; but you'd rather use the trace.idl script  



;
; prepare user_defined namelist part
;
param_namelist=[ nam1, "/", " ", "&datetime", nam2, nam3, nam4, "/", " ", "&interp", nam5 ]

if (tabnam(0) ne 999.) then begin
  ntot = ( 1+tabnam(1)-tabnam(0) ) / float(tabnam(2))
  ;ntot = ( tabnam(1)-tabnam(0) ) / float(tabnam(2))
  ntot = ceil(ntot) 
  temparr =  tabnam(0)+findgen(ntot)*tabnam(2)
endif else begin
  temparr = tabnam[1:n_elements(tabnam)-1] 
endelse


        if (n_elements(temparr) le 5.) then begin
		levelstr=[ " interp_levels = ", '   '+string(temparr,'(F7.3)'), "/"]
	endif else begin
		nnnn = floor(n_elements(temparr)/3.)
		ye = strarr(nnnn)
	        for i=0,nnnn-2 do begin
		  ye[i] = string(temparr[3*i],'(F7.3)')+','+string(temparr[3*i+1],'(F7.3)')+','+string(temparr[3*i+2],'(F7.3)')+','
		endfor
		ye[nnnn-1] = string(temparr[3*(nnnn-1)],'(F7.3)')+','+string(temparr[3*(nnnn-1)+1],'(F7.3)')+','+string(temparr[3*(nnnn-1)+2],'(F7.3)')
		levelstr=[ "interp_levels = ", ye, "/"]
	endelse

param_namelist=[param_namelist,levelstr]
param_namelist=transpose(param_namelist)

;
; generate file (not available in demo mode)
;
SPAWN, 'touch toto ; \rm toto'
OPENW, 1, 'toto'
PRINTF, 1, param_namelist
CLOSE, 1
	;;
	;; fix for demo mode
	;;	
	;SPAWN, 'touch toto ; \rm toto ; touch toto'
	;for i=0, n_elements(param_namelist)-1 do begin
	;	SPAWN, 'echo "'+param_namelist[i]+'" >> toto'
	;	print, param_namelist[i]
	;endfor	
instru = 'sed s/zenest/'+nest+'/g namelist.ARWpost_ref > namelist.ARWpost_ref_tmp'
instru = instru + ' ; ' + 'rm -rf zefolder ; ln -sf '+datafolder+' zefolder'
instru = instru + ' ; ' + 'touch namelist.ARWpost ; \rm namelist.ARWpost ; mv namelist.ARWpost_ref_tmp namelist.ARWpost'
if (grads eq 'true') then instru = instru + ' ; ' + "sed s/'idl'/'grads'/g namelist.ARWpost > yeye ; \mv yeye namelist.ARWpost"
instru = instru + ' ; ' + 'cat toto >> namelist.ARWpost ; \rm toto'

;
; extract data from netcdf file
;
instru = instru + ' ; ' + './ARWpost.exe'

;
; patch for correct time in plots
;
instru = instru + " ; " + "touch input_tmp.ctl ; \rm input_tmp.ctl ; touch timefil ; \rm timefil"
	;;; this line induces problems with variables with Z in the name
if (grads ne 'true') then instru = instru + " ; " + "sed s/'Z'/' 0 Z'/g input.ctl > input_tmp.ctl ; mv -f input_tmp.ctl input.ctl"
instru = instru + " ; " + "grep interval_seconds namelist.ARWpost | tail -n 1 | awk '{print $3}' > timefil ; cat timefil"
SPAWN, instru

;;old ARWpost
;SPAWN, "sed s/'*'/'9'/g input.dat > yeah ; mv -f yeah input.dat"


;; fix debile
SPAWN, "sed s/'SWDOWN 0 Z'/'SWDOWNZ'/g input.ctl > yeah ; mv -f yeah input.ctl"

;
; ok
;
print, 'ok...'
print, 'if the extraction was successful, '
print, 'you can use the input.dat and input.ctl files' 

if (grads eq 'true') then stop

END
