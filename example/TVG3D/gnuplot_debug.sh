#!/opt/homebrew/bin/gnuplot
    #========================header========================================
    set termoption dashed
    #set term pngcairo dashed
    set   autoscale                        
    unset log                              
    unset label                          
    set grid
    set xtic auto                         
    set ytic auto                         
    set terminal png 
    set term png size 900,600

    unset colorbox
    set style line 1 lt 1 lw 2.0 pt  1 ps 0.5 lc rgb "#000000"    
    set style line 2 lt 2 lw 2.0 pt  2 ps 0.5 lc rgb "#ff0000"  
    set style line 3 lt 3 lw 2.0 pt  8 ps 0.5 lc rgb "#0000ff"   
    set style line 4 lt 4 lw 2.0 pt  4 ps 0.5 lc rgb "#08000" 
    set style line 5 lt 5 lw 2.0 pt  3 ps 0.5 lc rgb "#00ffff" 
    set style line 6 lt 6 lw 2.0 pt  6 ps 0.5 lc rgb "#ff00ff" 
    set style line 7 lt 7 lw 2.0 pt 12 ps 0.5 lc rgb "#ffa500" 


    #================shared variables===================================
    Get_y_midp_C2P_3D="debug_Get_y_midp_C2P_3D.dat"



    #set ylabel "yc"
    #set xrange[100.0:]
    
    #=======================debug_init====================================
    set ylabel "y index"
    set output 'debug_Get_y_midp_C2P_3D.png'
    plot  Get_y_midp_C2P_3D using 2:1 with p ls 1 title "origin", \
          Get_y_midp_C2P_3D using 3:($1-0.5) with p ls 2 title "interpolation"
