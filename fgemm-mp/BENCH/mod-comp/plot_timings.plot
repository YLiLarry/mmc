set terminal postscript eps color enhanced font ", 18"
#
set style line 1  lt 1 lc 2 lw 3 pt 4 ps 2
set style line 2  lt 1 lc 3 lw 3 pt 5 ps 2 
#
set output 'power-proj.eps'
set xlabel "Polynomial degree"
set ylabel "Time in seconds"
plot 'power-proj.dat' using 1:2 w lp ls 1 title 'NTL',\
     'power-proj.dat' using 1:3 w lp ls 2 title 'FFLAS'
#---------------------------------------------------
set output 'mod-comp.eps'
set xlabel "Polynomial degree"
set ylabel "Time in seconds"
plot 'mod-comp.dat' using 1:2 w lp ls 1 title 'NTL',\
     'mod-comp.dat' using 1:3 w lp ls 2 title 'FFLAS'
#---------------------------------------------------
set output 'min-poly.eps'
set xlabel "Polynomial degree"
set ylabel "Time in seconds"
plot 'min-poly.dat' using 1:2 w lp ls 1 title 'NTL',\
     'min-poly.dat' using 1:3 w lp ls 2 title 'FFLAS'
#---------------------------------------------------
set output 'factoring.eps'
set xlabel "Polynomial degree"
set ylabel "Time in seconds"
plot 'factoring.dat' using 1:2 w lp ls 1 title 'NTL',\
     'factoring.dat' using 1:3 w lp ls 2 title 'FFLAS'
	
