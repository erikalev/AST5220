  �   >   k820309    `          17.0        �ąZ                                                                                                           
       bs_mod.f90 BS_MOD                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                           16                                                                                                  16                                                                                                  8#         @                                                      #Y 	   #DYDX 
   #XS    #HTOT    #NSTEP    #YOUT    #DERIVS           0  
 @                              	                   
              &                                                     
                                 
                   
              &                                                     
                                      
                
                                      
                
                                                      D @                                                 
               &                                           #         @                                       	               #X    #Y    #DYDX                                  
                                     
                
                                                   
              &                                                                                                       
               &                                           #         @                                                      #A    #B              
D@                                                 
               &                                                     
D                                                   
               &                                           #         @                                                    	   #Y    #DYDX    #X    #HTRY    #EPS    #YSCAL    #HDID    #HNEXT    #DERIVS           0  
D@                                                 
 
              &                                                     
  @                                                 
              &                                                     
D @                                   
                 
                                      
                
                                      
                
                                                    
              &                                                     D                                     
                 D                                     
       #         @    @                                  	               #X     #Y !   #DYDX "                                   
                                      
                
                                !                   
              &                                                                                    "                   
 	              &                                           (        `�                              #                                         #ARR $   #SEED %   p          H r &     7SO p        j            j                                      H r &     7SO p        j            j                                                           0  
 @                              $                                 &                                                     
 @                              %           (        `                                '                    "                    #N (   #M )     p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                                                    (                                                      )            (        `                               *                    &                
    #A +   #B ,     p        H r -     7
S
O p        j            j                                  p          H r -     7
S
O p        j            j                                    H r -     7
S
O p        j            j                                      H r -     7
S
O p        j            j                                    H r -     7
S
O p        j            j                                                           0  
 @                              +                   
 #             &                                                  0  
 @                              ,                   
 $             &                                           (        `                               .                    *                
    #A /   #B 0     p        H r 1     7
S
O p        j            j                                  p          H r 1     7
S
O p        j            j                                    H r 1     7
S
O p        j            j                                      H r 1     7
S
O p        j            j                                    H r 1     7
S
O p        j            j                                                           0  
 @                              /                   
 '             &                                                  0  
 @                              0                   
 (             &                                           (        `                               2                    ,                
    #FIRST 3   #INCREMENT 4   #N 5   p          5 O p            5 O p                                    
                                 3     
                
                                 4     
                
  @                              5           #         @                                  6                    #IEST 7   #XEST 8   #YEST 9   #YZ :   #DY ;             
                                 7                     
                                 8     
                
                                 9                   
              &                                                  0  D@                              :                   
               &                                                     D                                ;                   
               &                                           %         @                               <                           #ARR =             
                                 =                   
 -             &                                                                                     1     SIZE                                           -     SIZE                                           &     SIZE    �         fn#fn    �   @   J   HEALPIX_TYPES "   �   p       I4B+HEALPIX_TYPES !   j  p       DP+HEALPIX_TYPES "   �  p       LGT+HEALPIX_TYPES    J  r       NPAR_CUMSUM    �  r       NPAR_ARTH    .  q       NPAR2_ARTH    �  �       MMID    +  �   a   MMID%Y    �  �   a   MMID%DYDX    C  @   a   MMID%XS    �  @   a   MMID%HTOT    �  @   a   MMID%NSTEP      �   a   MMID%YOUT    �  t      MMID%DERIVS      @   a   MMID%DERIVS%X    C  �   a   MMID%DERIVS%Y !   �  �   a   MMID%DERIVS%DYDX    [  V       SWAP    �  �   a   SWAP%A    =	  �   a   SWAP%B    �	  �       BSSTEP    h
  �   a   BSSTEP%Y    �
  �   a   BSSTEP%DYDX    �  @   a   BSSTEP%X    �  @   a   BSSTEP%HTRY       @   a   BSSTEP%EPS    @  �   a   BSSTEP%YSCAL    �  @   a   BSSTEP%HDID      @   a   BSSTEP%HNEXT    L  v      BSSTEP%DERIVS     �  @   a   BSSTEP%DERIVS%X       �   a   BSSTEP%DERIVS%Y #   �  �   a   BSSTEP%DERIVS%DYDX      w      CUMSUM    �  �   a   CUMSUM%ARR      @   a   CUMSUM%SEED    ]        UPPER_TRIANGLE !   w  @   a   UPPER_TRIANGLE%N !   �  @   a   UPPER_TRIANGLE%M    �  �      OUTERDIFF    �  �   a   OUTERDIFF%A    U  �   a   OUTERDIFF%B    �  �      OUTERPROD    �  �   a   OUTERPROD%A    ?  �   a   OUTERPROD%B    �  �       ARTH    �  @   a   ARTH%FIRST    �  @   a   ARTH%INCREMENT       @   a   ARTH%N    `  v       PZEXTR    �  @   a   PZEXTR%IEST      @   a   PZEXTR%XEST    V  �   a   PZEXTR%YEST    �  �   a   PZEXTR%YZ    n  �   a   PZEXTR%DY    �  Y       IMINLOC    S  �   a   IMINLOC%ARR    �  =      OUTERPROD%SIZE       =      OUTERDIFF%SIZE    Y   =      CUMSUM%SIZE 