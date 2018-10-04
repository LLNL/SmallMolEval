#!/usr/local/bin/gnuplot -c 
#set terminal png transparent size 240,180 enhanced font 'Verdana,10'
set terminal png transparent size ARG2,ARG3 enhanced font 'Verdana,10'
o = sprintf("%s%s", ARG1,".png")
set output o
set key left top
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic 0,1,1                          # set ytics automatically
set title ARG1 offset 0,-1
#set xlabel "t" offset 0,.75
set xr [0.0:5.0]
set yr [0.0:1.0]
i = sprintf("%s%s", ARG1, ".out.txt")
plot i using 1:2 title 'G(t)' with lines lw 3, i using 1:3 title 'F(t)' with lines lw 3
