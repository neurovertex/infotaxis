if (!exists("logfile")) logfile='log.txt'
if (!exists("outfile")) outfile='gnuplot_path.png'
resolution=system('head '. logfile) + 0
boundary=100/resolution
sourcex=90
sourcey=90
set terminal pngcairo size 800,600 enhanced font 'Roboto,12'
set output outfile
set xlabel "X(m)"
set ylabel "Y(m)"
set xrange [0:(100/resolution)]
set yrange [0:(100/resolution)]

plot 	logfile using ($1/resolution):($2/resolution) w l title "Path"