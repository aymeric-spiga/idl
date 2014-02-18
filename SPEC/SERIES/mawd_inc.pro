;;;;;*************************************** PLOT TITLES
        title_user 	= 	'Water vapor column [precip. microns]' 
        title_axis      =       ['L!DS!N','North Latitude']

;;;;;*************************************** COLOR TABLES
        flag_cb         =       'true'
	poscb		=	0.75
        format          =       '(I0)' 
        colors          =       32
        pal             =       33;22;4;33              ;; GOOD: 4, 18, 22, 16, 37, 33, 39, 6, 11

;;;;;*************************************** FILL LIMITS
        minfield_init   =       0.
        maxfield_init   =       100. ;120. ;65.
        ndiv            =       10 ;12 ;13

	;;;;;*************************************** LIMIT TRICKS
        ;;;;;********************* must always follow FILL LIMITS
        	missing_value=1.e30
        	lim_max = maxfield_init & w=where((what_I_plot ge lim_max) and (what_I_plot le 1e9)) & if (w[0] ne -1) then what_I_plot[w]=lim_max
        	lim_min = minfield_init & w=where(what_I_plot le lim_min) & if (w[0] ne -1) then what_I_plot[w]=lim_min
        	lim_blank = 1.0 & w=where(abs(what_I_plot) le lim_blank) & if (w[0] ne -1) then what_I_plot[w]=missing_value

;;;;;*************************************** CONTOUR
	overcontour	=	overcontour
        lev             =       minfield_init + (maxfield_init - minfield_init)*findgen(ndiv+1)/float(ndiv)

;;;;;*************************************** AXIS
        isotropic       =       'true'
        intervalx       =       15.0
        intervaly       =	15.0

;;;;;*************************************** MAP LIMITS
	windowx 	= 	[30.,150.]
        windowx         =       [30.,1500.]
	windowy 	= 	[-15.,90.]

;;;;; MULTIYEAR
	isotropic       =       'false'
        windowx         =       [0.,2880.] ;[0.,3600.] ;[1500.,4000.]
        intervalx       =       360.0
        windowy         =       [-90.,90.]
        intervaly       =       30.0
        lev		=       [600.] ;[100.]
