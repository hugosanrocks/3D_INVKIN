  "  O   k820309    �          2021.7.1    ���c                                                                                                          
       ../delaunay/delaunay_mod.f90 DELAUNAY #         @                                                       #NODE_NUM    #NODE_XY                                                                   D                                                    	     p          p          5 � p        r        p          5 � p        r                      #         @                                                       #TRIANGLE_NUM    #TRIANGLE_NODE                                                                                                                            p          p          5 � p        r        p          5 � p        r                      %         @                                                          #X0    #Y0 	   #X1 
   #Y1    #X2    #Y2    #X3    #Y3                                                   	                                                 	     	                                                 
     	                                                      	                                                      	                                                      	                                                      	                                                      	       #         @                                                       #POINT_NUM    #POINT_XY    #TRI_NUM    #TRI_VERT    #TRI_NABE              D @                                                   D @                                                  	     p          p          5 � p        r        p          5 � p        r                                D @                                                   D @                                                       p          p           5 � p        r    n                                       2    p           5 � p        r    n                                      2                                    D @                                                       p          p           5 � p        r    n                                       2    p           5 � p        r    n                                      2                           #         @                                                      #N    #A    #INDX                                                                                                                       	     p          p          5 � p        r        p          5 � p        r                               D                                                         p          5 � p        r        5 � p        r                      #         @                                                      #N    #P    #A              D @                                                   D @                                                       p          5 � p        r        5 � p        r                               D                                                    	 
    p          p          5 � p        r        p          5 � p        r                      %         @                                                          #XU    #YU     #XV1 !   #YV1 "   #XV2 #   #YV2 $   #DV %                                                  	                                                       	                                                 !     	                                                 "     	                                                 #     	                                                 $     	                  @                              %     	       #         @                                  &                    #X '   #Y (   #POINT_NUM )   #POINT_XY *   #TRI_NUM +   #TRI_VERT ,   #TRI_NABE -   #LTRI .   #LEDG /   #RTRI 0   #REDG 1             D @                              '     	                 D @                              (     	                                                 )                     D @                              *                    	     p          p          5 � p        r )       p          5 � p        r )                                                               +                                                     ,                         p          p          5 � p        r +       p          5 � p        r +                                                              -                         p          p          5 � p        r +       p          5 � p        r +                               D                                .                      D                                /                      D                                0                      D                                1            #         @                                  2                    #I 3   #TOP 4   #BTRI 5   #BEDG 6   #POINT_NUM 7   #POINT_XY 8   #TRI_NUM 9   #TRI_VERT :   #TRI_NABE ;   #STACK <   #IERR =                                             3                      D                                4                      D                                5                      D                                6                                                      7                     D @                              8                    	     p          p          5 � p        r 7       p          5 � p        r 7                                                               9                     D                                :                         p          p          5 � p        r 9       p          5 � p        r 9                              D                                ;                         p          p          5 � p        r 9       p          5 � p        r 9                              D                                <                         p          5 � p        r 7       5 � p        r 7                               D                                =            #         @                                  >                    #N ?   #P @             D @                              ?                     D @                              @                     	    p          5 � p        r ?       5 � p        r ?                     %         @                               A                           #I B   #J C              @                              B                       @                              C            %         @                               D                           #X E                                             E            %         @                               F                           #IVAL G   #ILO H   #IHI I                                             G                       @                              H                       @                              I            #         @                                  J                    #N K   #P L   #BASE M   #IERROR N                                             K                                                     L                         p          5 � p        r K       5 � p        r K                                                               M                      D                                N               �   .      fn#fn    �   c       READ_NODES $   1  @   a   READ_NODES%NODE_NUM #   q  �   a   READ_NODES%NODE_XY     E  m       WRITE_TRIANGLES -   �  @   a   WRITE_TRIANGLES%TRIANGLE_NUM .   �  �   a   WRITE_TRIANGLES%TRIANGLE_NODE    �  �       DIAEDG    V  @   a   DIAEDG%X0    �  @   a   DIAEDG%Y0    �  @   a   DIAEDG%X1      @   a   DIAEDG%Y1    V  @   a   DIAEDG%X2    �  @   a   DIAEDG%Y2    �  @   a   DIAEDG%X3      @   a   DIAEDG%Y3    V  �       DTRIS2 !   �  @   a   DTRIS2%POINT_NUM     $  �   a   DTRIS2%POINT_XY    �  @   a   DTRIS2%TRI_NUM     8  F  a   DTRIS2%TRI_VERT     ~	  F  a   DTRIS2%TRI_NABE )   �
  `       R82VEC_SORT_HEAP_INDEX_A +   $  @   a   R82VEC_SORT_HEAP_INDEX_A%N +   d  �   a   R82VEC_SORT_HEAP_INDEX_A%A .   8  �   a   R82VEC_SORT_HEAP_INDEX_A%INDX    �  ]       R82VEC_PERMUTE !   I  @   a   R82VEC_PERMUTE%N !   �  �   a   R82VEC_PERMUTE%P !   =  �   a   R82VEC_PERMUTE%A      �       LRLINE    �  @   a   LRLINE%XU    �  @   a   LRLINE%YU      @   a   LRLINE%XV1    ]  @   a   LRLINE%YV1    �  @   a   LRLINE%XV2    �  @   a   LRLINE%YV2      @   a   LRLINE%DV    ]  �       VBEDG    !  @   a   VBEDG%X    a  @   a   VBEDG%Y     �  @   a   VBEDG%POINT_NUM    �  �   a   VBEDG%POINT_XY    �  @   a   VBEDG%TRI_NUM    �  �   a   VBEDG%TRI_VERT    �  �   a   VBEDG%TRI_NABE    �  @   a   VBEDG%LTRI    �  @   a   VBEDG%LEDG      @   a   VBEDG%RTRI    ]  @   a   VBEDG%REDG    �  �       SWAPEC    d  @   a   SWAPEC%I    �  @   a   SWAPEC%TOP    �  @   a   SWAPEC%BTRI    $  @   a   SWAPEC%BEDG !   d  @   a   SWAPEC%POINT_NUM     �  �   a   SWAPEC%POINT_XY    x  @   a   SWAPEC%TRI_NUM     �  �   a   SWAPEC%TRI_VERT     �  �   a   SWAPEC%TRI_NABE    `  �   a   SWAPEC%STACK      @   a   SWAPEC%IERR    T  V       PERM_INVERSE    �  @   a   PERM_INVERSE%N    �  �   a   PERM_INVERSE%P    �  ^       I4_MODP    �  @   a   I4_MODP%I    <  @   a   I4_MODP%J    |  W       I4_SIGN    �  @   a   I4_SIGN%X      l       I4_WRAP      @   a   I4_WRAP%IVAL    �  @   a   I4_WRAP%ILO    �  @   a   I4_WRAP%IHI    ?   l       PERM_CHECK    �   @   a   PERM_CHECK%N    �   �   a   PERM_CHECK%P     �!  @   a   PERM_CHECK%BASE "   �!  @   a   PERM_CHECK%IERROR 