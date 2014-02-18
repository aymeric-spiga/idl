PRO plot_vector, u, v, x, y, type=type, angle=angle, head_len=head_len,$
                 maxmag=maxmag, align=align, length=length,$
		 ref_text=ref_text, color=color, cstyle=cstyle
;Procedure to calculate and plot a vector
 
mylen=300.0 ; length of default arrow in pixels
rev=1.0    
x0=0.0
y0=0.0
x1=u/maxmag*mylen*length
y1=v/maxmag*mylen*length
dx=x1-x0
IF (dx LT 0.0) THEN rev=-1.0
dy=y1-y0
r=SQRT(dx^2+dy^2)
theta=ATAN(dy/dx)     
phi=angle*!dtor
rfrac=head_len 
x2=x1-r*rfrac*rev*COS(theta-phi)
y2=y1-r*rfrac*rev*SIN(theta-phi)
x3=x1-r*rfrac*rev*COS(theta+phi)
y3=y1-r*rfrac*rev*SIN(theta+phi)
x4=x1-rfrac/2*r*rev*COS(theta)
y4=y1-rfrac/2*r*rev*SIN(theta)
;Calculate intersection of vector shaft and head points either
;side of the shaft - see
;http://astronomy.swin.edu.au/~pbourke/geometry/lineline2d
;for more details
ua=((x3-x2)*(y0-y2)-(y3-y2)*(x0-x2))/$
   ((y3-y2)*(x1-x0)-(x3-x2)*(y1-y0))
x5=x0+ua*(x1-x0)
y5=y0+ua*(y1-y0) 
               
;Plot vectors in data space - cstyle=0
;Position in device coordinates and then convert to data coordinates
IF (cstyle EQ 0) THEN BEGIN
 pt1=CONVERT_COORD(x,y, /DATA, /TO_DEVICE)
 xpts=[x0,x1,x2,x3,x4,x5]+pt1[0]-align*dx
 ypts=[y0,y1,y2,y3,y4,y5]+pt1[1]-align*dy
 pts=CONVERT_COORD(xpts,ypts, /DEVICE, /TO_DATA)
 xpts=pts[0,*]
 ypts=pts[1,*]
 x0=xpts[0]
 x1=xpts[1]
 x2=xpts[2]
 x3=xpts[3]
 x4=xpts[4]
 x5=xpts[5]
 y0=ypts[0]
 y1=ypts[1]
 y2=ypts[2]
 y3=ypts[3]
 y4=ypts[4]
 y5=ypts[5]   
   
	       
; Plot the vectors omiting any vectors with NaNs
 z=[xpts, ypts] 	 
 IF (TOTAL(FINITE(z)) EQ 12) THEN BEGIN
  IF (N_ELEMENTS(TYPE) EQ 0) THEN TYPE=0
  IF (TYPE EQ 0) THEN BEGIN	 
    PLOTS, [x0,x1,x2,x1,x3], [y0,y1,y2,y1,y3], COLOR=color
  ENDIF
  IF (TYPE EQ 1) THEN BEGIN	 
    PLOTS, [x0,x5,x3,x1,x2,x5], [y0,y5,y3,y1,y2,y5], COLOR=color
  ENDIF  
  IF (TYPE EQ 2) THEN BEGIN	 
    PLOTS, [x0,x4,x2,x1,x3,x4], [y0,y4,y2,y1,y3,y4], COLOR=color
  ENDIF
  IF (TYPE EQ 3) THEN BEGIN	 
    PLOTS, [x0,x1], [y0,y1], COLOR=color
    POLYFILL, [x1,x2,x3,x1], [y1,y2,y3,y1], COLOR=color
  ENDIF
  IF (TYPE EQ 4) THEN BEGIN	 
    PLOTS, [x0,x4], [y0,y4], COLOR=color
   POLYFILL, [x1,x2,x4,x3,x1], [y1,y2,y4,y3,y1], COLOR=color
  ENDIF
 ENDIF
ENDIF	


;Plot reference vector - cstyle=1
;Position in device coordinates and then convert to data coordinates
IF (cstyle EQ 1) THEN BEGIN
 pt1=CONVERT_COORD(x,y, /NORMAL, /TO_DEVICE)
 xpts=[x0,x1,x2,x3,x4,x5]+pt1[0]
 ypts=[y0,y1,y2,y3,y4,y5]+pt1[1]
 x0=xpts[0]
 x1=xpts[1]
 x2=xpts[2]
 x3=xpts[3]
 x4=xpts[4]
 x5=xpts[5]
 y0=ypts[0]
 y1=ypts[1]
 y2=ypts[2]
 y3=ypts[3]
 y4=ypts[4]
 y5=ypts[5]     
	       
; Plot the vectors omiting any vectors with NaNs
 z=[xpts, ypts] 	 
 IF (TOTAL(FINITE(z)) EQ 12) THEN BEGIN
  IF (N_ELEMENTS(TYPE) EQ 0) THEN TYPE=0
  IF (TYPE EQ 0) THEN BEGIN	 
    PLOTS, [x0,x1,x2,x1,x3], [y0,y1,y2,y1,y3], COLOR=color, /DEVICE
  ENDIF
  IF (TYPE EQ 1) THEN BEGIN	 
    PLOTS, [x0,x5,x3,x1,x2,x5], [y0,y5,y3,y1,y2,y5], COLOR=color, /DEVICE
  ENDIF  
  IF (TYPE EQ 2) THEN BEGIN	 
    PLOTS, [x0,x4,x2,x1,x3,x4], [y0,y4,y2,y1,y3,y4], COLOR=color, /DEVICE
  ENDIF
  IF (TYPE EQ 3) THEN BEGIN	 
    PLOTS, [x0,x1], [y0,y1], COLOR=color, /DEVICE
    POLYFILL, [x1,x2,x3,x1], [y1,y2,y3,y1], COLOR=color, /DEVICE
  ENDIF
  IF (TYPE EQ 4) THEN BEGIN	 
    PLOTS, [x0,x4], [y0,y4], COLOR=color, /DEVICE
   POLYFILL, [x1,x2,x4,x3,x1], [y1,y2,y4,y3,y1], COLOR=color, /DEVICE
  ENDIF
 ENDIF
 
 ;Add the reference vector text
 IF (N_ELEMENTS(REF_TEXT) EQ 0) THEN ref_text=' '
 XYOUTS, x0-100, y0-100, ref_text, ALIGNMENT=1.0, /DEVICE
ENDIF
	
END		
