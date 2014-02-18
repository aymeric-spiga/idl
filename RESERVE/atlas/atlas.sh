#! /bin/bash

\rm atlas.log
echo "print, '****************************yy'
zelevel='1'
firstextract = 'true'
zeatlas='/tmp7/aslmd/atlas-web_map_model_level_U_HGT_UV.def' 
@atlas.pro
print, '****************************yy'
zelevel='2'
firstextract = 'true'
zeatlas='/tmp7/aslmd/atlas-web_map_model_level_U_HGT_UV.def' 
@atlas.pro
print, '****************************yy'
zelevel='3'
firstextract = 'true'
zeatlas='/tmp7/aslmd/atlas-web_map_model_level_U_HGT_UV.def' 
@atlas.pro
print, '****************************yy'
zelevel='4'
firstextract = 'true'
zeatlas='/tmp7/aslmd/atlas-web_map_model_level_U_HGT_UV.def' 
@atlas.pro
print, '****************************yy'
zelevel='5'
firstextract = 'true'
zeatlas='/tmp7/aslmd/atlas-web_map_model_level_U_HGT_UV.def' 
@atlas.pro
print, '****************************yy'
zelevel='6'
firstextract = 'true'
zeatlas='/tmp7/aslmd/atlas-web_map_model_level_U_HGT_UV.def' 
@atlas.pro
exit" | /opt/idl-6.4/idl/bin/idl >> atlas.log 2>& 1 


#############################
#
# have you modified atlas.pro ?
# sed s/'U_HGT_UV'/'W_HGT'/g atlas.pro > ye.pro
#
#echo atlas.pro should have been modified
#echo you should have had generated input files
#

