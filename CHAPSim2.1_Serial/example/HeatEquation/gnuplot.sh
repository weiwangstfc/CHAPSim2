#!/usr/bin/gnuplot
    #========================header========================================
    set termoption dashed
    #set term pngcairo dashed
    set bmargin 5
    set   autoscale                        
    unset log                              
    unset label                          
    set grid
    set xtic auto                         
    set ytic auto                         
    set terminal png 
    set term png size 800,600 enhanced font "Helvetica,14"
    
    set border lw 2.0
    set tics scale 1.8
    set tics   font ", 18"
    set key    font ", 15"
    set xlabel font ", 18"
    set ylabel font ", 18"

    unset colorbox
 
 #========================plot========================================
    set output 'Temporal_development.png'
    set key right top
    set ylabel "u(x)"
    set xlabel "x"
    plot for [i=0:30000:1000] 'Plot_Burgers_profile'.i.'.dat' using 1:2 with lp title 'iter='.i
 
 
 
 
 
 
 
 
 
  
