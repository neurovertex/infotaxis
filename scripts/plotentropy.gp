if (!exists("logfile")) logfile='log.txt'
resolution=system('head '. logfile) + 0
sourcex=90
sourcey=90
set terminal pngcairo size 800,600 enhanced font 'Roboto,12'
set output "gnuplot.png"
set xlabel "iteration"
set ylabel sprintf("Distance (m * %.2f)", 1/resolution)
set y2label "Entropy | Detections"
set ytics nomirror
set y2tics 1
set grid xtics y2tics
plot 	logfile using (sqrt(($1-sourcex)*($1-sourcex) + ($2-sourcey)*($2-sourcey))) with lines title "Distance", \
	 	logfile using 3 with lines title "Entropy" axes x1y2, \
	 	logfile using 4 with histeps title "Detections" axes x1y2