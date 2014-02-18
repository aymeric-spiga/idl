;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
saveplot        = './COMP_REAL_50m_145_145_201_12km/'
setting         = 'REAL_50m_145_145_201_12km/'
;cases          = ['case_A', 'case_B', 'case_C', 'case_I']
;cases          = ['case_A', 'case_B', 'case_C', 'case_HIGH']
cases          = ['case_A', 'case_B', 'case_C', 'case_I', 'case_HIGH']
;cases          = ['case_A', 'case_B', 'case_C', 'case_I']
paths           = setting + cases
tttitle         = ['[case a]  T!Ds,max!N=284 K   h= -3.9  km' $
                  ,'[case b]  T!Ds,max!N=280 K   h= -3.5  km' $
                  ,'[case c]  T!Ds,max!N=279 K   h=+2.4 km' $
                  ,'[case i]   T!Ds,max!N=279 K   h= -0.5  km' $
                  ,'[case z]  T!Ds,max!N=268 K   h=+8.4 km' ]
tttitles        = ['[a]  T!Ds,m!N=284 K   h= -3.9  km' $
                  ,'[b]  T!Ds,m!N=280 K   h= -3.5  km' $
                  ,'[c]  T!Ds,m!N=279 K   h=+2.4 km' $
                  ,'[i]   T!Ds,m!N=279 K   h= -0.5  km' $
                  ,'[z]  T!Ds,m!N=268 K   h=+8.4 km' ]
shorttt         = ['a','b','c','i','z' ]


        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;saveplot        = './COMP_case_C/'
        ;setting        = [       $
        ;;                'REAL_100m_145_145_201_15km/', $
        ;;                'REAL_100m_101_101_201_15km/', $
        ;;                'REAL_100m_101_101_151_15km/', $
        ;;                'REAL_100m_73_73_201_15km/', $
        ;;                'REAL_100m_73_73_151_15km/', $
        ;;                      'REAL_100m_45_45_71_10km/', $
        ;;                'REAL_100m_145_145_201_20km/', $
        ;;                'REAL_100m_45_45_71_10km/', $
        ;;                'REAL_100m_45_45_71_10km/', $
        ;;;                'REAL_50m_201_201_201_12km/', $
        ;;;                'REAL_50m_241_241_201_12km/', $
        ;;'REAL_50m_181_181_201_12km/', $
        ;;      'REAL_50m_145_145_251_10km/', $
        ;;              'REAL_100m_145_145_101_12km/', $
        ;;              'REAL_100m_241_241_151_15km/', $
        ;'REAL_100m_181_181_201_15km/', $
        ;                'REAL_50m_145_145_201_12km/'      $       ;;ref 
        ;               ]
        ;cases          = 'case_C'
        ;paths           = setting + cases
        ;;paths = [ paths, '/home/spiga/data/current3/']
        ;tttitle = paths
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; pour les calculs des lois de similitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;saveplot        = './'
;setting         = 'REAL_50m_145_145_201_12km/'
;cases           = ['case_A', 'case_B', 'case_C', 'case_I']
;;;cases           = ['case_C', 'case_B']
;;cases           = ['case_B', 'case_C', 'case_HIGH']
;cases           = ['case_A', 'case_B', 'case_C', 'case_I','case_HIGH']
;paths           = setting + cases
;;;paths = [ paths, './GUSEV/']
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;saveplot        = './'
;paths = [ 'REAL_50m_145_145_201_12km/case_A', '/home/spiga/data/current3/']

;;tttitle = paths

