;;
;; ATLAS EXOMARS STUDY
;;

zelevel='2'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-web_map_model_level_U_HGT_UV.def'   ;; event. modifier si on ne veut pas recalculer les .idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='2'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-closer_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='2'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-ellipse_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
SPAWN, 'mkdir zefolder/idl_save_2 ; mv zefolder/*.idl zefolder/idl_save_2/'
;

zelevel='4'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-web_map_model_level_U_HGT_UV.def'   ;; event. modifier si on ne veut pas recalculer les .idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='4'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-closer_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='4'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-ellipse_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
SPAWN, 'mkdir zefolder/idl_save_4 ; mv zefolder/*.idl zefolder/idl_save_4/'
;

zelevel='6'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-web_map_model_level_U_HGT_UV.def'   ;; event. modifier si on ne veut pas recalculer les .idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='6'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-closer_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='6'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-ellipse_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
SPAWN, 'mkdir zefolder/idl_save_6 ; mv zefolder/*.idl zefolder/idl_save_6/'
;

zelevel='7'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-web_map_model_level_U_HGT_UV.def'   ;; event. modifier si on ne veut pas recalculer les .idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='7'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-closer_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='7'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-ellipse_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
SPAWN, 'mkdir zefolder/idl_save_7 ; mv zefolder/*.idl zefolder/idl_save_7/'
;

zelevel='8'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-web_map_model_level_U_HGT_UV.def'   ;; event. modifier si on ne veut pas recalculer les .idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='8'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-closer_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='8'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-ellipse_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
SPAWN, 'mkdir zefolder/idl_save_8 ; mv zefolder/*.idl zefolder/idl_save_8/'
;

zelevel='18'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-web_map_model_level_U_HGT_UV.def'   ;; event. modifier si on ne veut pas recalculer les .idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='18'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-closer_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

zelevel='18'
zeatlas='/home/aslmd/MERIDIANI_EXOMARS/atlas-ellipse_web_map_model_level_U_HGT_UV.def'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
SPAWN, "ln -sf "+zeatlas+" atlas.def ; \rm gw.def ; sed s:'zeremplace':'"+zelevel+"':g atlas.def > gw.def"
@gwc
gw, 'w'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
SPAWN, 'mkdir zefolder/idl_save_18 ; mv zefolder/*.idl zefolder/idl_save_18/'
;

exit
