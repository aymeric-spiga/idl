pro genanim, num, folder, name

last=num-1

OPENW, 1, 'templist'

FOR i=0,last do begin
	;; on ne met plus de dossier car on nomme les images
	;printf, 1, 'modImages['+STRTRIM(STRING(i),2)+'] = "'+folder+'/',$
	;;	+STRTRIM(STRING(i+1000),2)+'.ps.png";'
printf, 1, 'modImages['+STRTRIM(STRING(i),2)+'] = "',$
       +name+STRTRIM(STRING(i+1000),2)+'.png";'
ENDFOR

printf, 1, 'first_image = '+string(0,'(I0)')+';'
printf, 1, 'last_image = '+string(last,'(I0)')+';'

CLOSE, 1

end
