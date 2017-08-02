set terminal png size 2048, 2048;
set output "Sample.heatmap.png";
set palette rgbformula 34, 35, 36 negative;
unset colorbox;
### Customize this using commented values at end of input data file
set xtics rotate ("ID1" 0, "ID2" NUMBER, etc.);
### Customize this
set ytics ("ID1" 0, "ID2" NUMBER, etc.);
set grid front;
plot [0:600][0:600] 'dat/BO121.heatmap.dat' using 1:2:(log($3)) notitle with image;
