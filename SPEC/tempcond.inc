;
; input: log-pressure altitude
;

latcond=5.9e5
tcond1mb=136.27
r=192.

bcond=1./tcond1mb
acond=r/latcond

;press=610.*exp(-column/10.)
press=columnp

overplot=1./(bcond-acond*alog(.0095*press))
overplot_column=column
