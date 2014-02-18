; vector.pro
; vector plotting routine 
; see http://cgam.nerc.ac.uk/graphics/idl.php for more details

; Calling procedure:
; VECTOR, U, V, X, Y
; U and V are the vector components as two dimensional arrays.
; X and Y are the vectors locations as one dimensional arrays.
;
; Optional parameters:
; LENGTH: Length of the vector, default=1 which is 300 pixels.
; OVERPLOT: Set this to overplot the current graphics plot.
; STRIDE: Pick every nth vector for plotting.
; TYPE:   Type of arrow head. 0-4, default=0
; ANGLE:  Angle of vector head in degrees, default=22.5.
; HEAD_LEN: Length of the arrow head relative to shaft, default=0.3
; ALIGN:  Alignment of arrow.  Default is 0.5, range=0.0 to 1.0.
; REF_MAG: Arrow magnitude.
; REF_POS: Position for reference vector in normalised coordinates [x,y].
; REF_TEXT: Text to place to the left of the reference vector.
; COLOR: Color index of vectors
		
PRO VECTOR, u,v,x,y, LENGTH = length,$
        Color=color, STRIDE=stride, ALIGN=align, $
	REF_MAG=ref_mag, TYPE=type, ANGLE=angle, HEAD_LEN=head_len,$
	REF_POS=ref_pos, REF_TEXT=ref_text, OVERPLOT=overplot, _EXTRA=extra

;Check mapping is cylindrical or not set - reject if otherwise.
IF ((!MAP.PROJECTION NE 0) AND (!MAP.PROJECTION NE 8)) THEN BEGIN
  PRINT, 'Mapping must be cylindrical or unset to use the vector routine'
  RETURN
ENDIF

;Basic checks of the input data
a=SIZE(u)
b=SIZE(v)
c=SIZE(x)
d=SIZE(y)
IF ((a[0] NE 2) OR (b[0] NE 2)) THEN BEGIN
  PRINT, 'u and v must be two dimensional'
  RETURN
ENDIF
IF (TOTAL(ABS(b[0:2]-a[0:2])) NE 0) THEN BEGIN
  PRINT, 'u and v must have the same size'
  RETURN
ENDIF
IF ((c[0] NE 1) AND (d[0] NE 1)) THEN BEGIN
  PRINT, 'x and y must be one dimensional'
  RETURN
ENDIF
IF ((c[1]-a[1]) NE 0) THEN BEGIN
  PRINT, 'u and x have mismatched sizes'
  RETURN
ENDIF
IF ((d[1]-b[2]) NE 0) THEN BEGIN
  PRINT, 'v and y have mismatched sizes'
  RETURN
ENDIF

;Initialise parameters if undefined
IF (N_ELEMENTS(STRIDE) EQ 0) THEN stride=0
IF N_ELEMENTS(LENGTH) EQ 0 THEN length=1.0
IF N_ELEMENTS(COLOR) EQ 0 THEN color = !P.COLOR
IF (N_ELEMENTS(ANGLE) EQ 0) THEN angle=22.5
IF (N_ELEMENTS(HEAD_LEN) EQ 0) THEN head_len=0.3
IF (N_ELEMENTS(TYPE) EQ 0) THEN TYPE=0
IF (N_ELEMENTS(ALIGN) EQ 0) THEN align=0.5	
IF (N_ELEMENTS(REF_TEXT) EQ 0) THEN ref_text=' '
IF (N_ELEMENTS(REF_MAG) EQ 0) THEN BEGIN
  maxmag=MAX(ABS(SQRT(u^2.+v^2.)))
ENDIF ELSE BEGIN
  maxmag=ref_mag
ENDELSE

;Setup the plot area if undefined
IF (NOT KEYWORD_SET(overplot)) THEN BEGIN
  xs=x[0]-(x(1)-x(0))
  xf=x[N_ELEMENTS(x)-1]+(x(1)-x(0))
  ys=y[0]-(y(1)-y(0))
  yf=y[N_ELEMENTS(y)-1]+(y(1)-y(0))  
  PLOT,[xs,xf],[ys,yf], XSTYLE=1, YSTYLE=1, /NODATA,$
       COLOR=color, _EXTRA = extra
ENDIF

;Do stride data reduction if needed	
IF (stride GT 1) THEN BEGIN
  mypts=FLTARR(a[1], a[2])
  mypts[*,*]=0.0	   
  FOR iy=0,a[2]-1,stride DO BEGIN
  FOR ix=0,a[1]-1,stride DO BEGIN
    IF ( ((ix/stride) EQ FIX(ix/stride)) AND $
         ((iy/stride) EQ FIX(iy/stride)) ) THEN mypts[ix,iy]=1.0
  ENDFOR
  ENDFOR
  pts=WHERE(mypts LT 1.0)
  u[pts]=0.0
  v[pts]=0.0
ENDIF 

;Main vector plotting loop
FOR ix=0, N_ELEMENTS(x)-1 DO BEGIN
FOR iy=0, N_ELEMENTS(y)-1 DO BEGIN
  PLOT_VECTOR, u(ix,iy), v(ix,iy), x(ix), y(iy), type=type, $
               angle=angle, head_len=head_len,$
               maxmag=maxmag, align=align, length=length,$
	       color=color, cstyle=0
ENDFOR
ENDFOR

;Plot a reference vector if requested
IF (N_ELEMENTS(REF_POS) NE 0) THEN BEGIN
  PLOT_VECTOR, maxmag, 0.0, ref_pos[0], ref_pos[1], type=type, $
               angle=angle, ref_text=ref_text, head_len=head_len,$
               maxmag=maxmag, align=align, length=length,$
	       color=color, cstyle=1
ENDIF

END

