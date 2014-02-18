;;
;; ATLAS EXOMARS STUDY
;;

;;;zelevel='1'
;;;firstextract = 'true'
;;;zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-web_map_model_level_W_HGT.def'   ;; event. modifier si on ne veut pas recalculer les .idl


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
if (firstextract eq 'true') then SPAWN, "echo 'dumb=999' >> gw.def" else SPAWN, "echo 'dumb=0' >> gw.def"  ;;  juste pour le premier
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SPAWN, "echo 'dumb=0' >> gw.def ; echo 'intervalx=1.0 ;; grepALL' >> gw.def ; echo 'intervaly=1.0 ;; grepALL' >> gw.def ; echo 'windowx=[-12.,00.] ;; grepMAP' >> gw.def ; echo 'windowy=[-08.,04.] ;; grepMAP' >> gw.def ; echo 'stride=3. ;; grepALL' >> gw.def ; sed s:'atlasMERIDIANI':'atlasMERIDIANI_closer':g gw.def | sed s:'false':'true':g > ye.def ; mv -f ye.def gw.def"  ;; attention dernier pas tres robuste, pour useidlfile
;@gwc
;gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SPAWN, "echo 'intervalx=0.5 ;; grepALL' >> gw.def ; echo 'intervaly=0.5 ;; grepALL' >> gw.def ; echo 'windowx=[-08.,-04.] ;; grepMAP' >> gw.def ; echo 'windowy=[-04.,00.] ;; grepMAP' >> gw.def ; echo 'stride=1. ;; grepALL' >> gw.def ; sed s:'atlasMERIDIANI_closer':'atlasMERIDIANI_ellipse':g gw.def > ye.def ; mv -f ye.def gw.def"
;@gwc
;gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





;;;;
;;;SPAWN, 'mkdir zefolder/idl_save_'+zelevel+' ; mv zefolder/*.idl zefolder/idl_save_'+zelevel+'/'

