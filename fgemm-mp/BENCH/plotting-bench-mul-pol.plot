set encoding iso_8859_1 

set xlabel "Polynomial degree" 
set ylabel "Time in seconds" offset 2
set key left Left top reverse



set logscale x 2
set format x "2^{%L}"
set logscale y 2
set format y "%.2f"

#set ytics (0,001,0.01,0.1,1,10,100,1000,10000,100000)
set grid ytics mytics


#---------------------------------------------*
set style line 1  lt 1 lc 1 lw 3 pt 2 ps 2
set style line 2  lt 1 lc 2 lw 3 pt 4 ps 2
set style line 3  lt 1 lc 3 lw 3 pt 5 ps 2 
set style line 4  lt 2 lc 4 lw 3 pt 1 ps 2 
#---------------------------------------------*



set terminal postscript enhanced landscape color 18
 set datafile separator '|'

dim="16 32 64 256 512"

do for [i=0:4] {
outfile=sprintf('bench-mul-pol-%d.eps',i);
set output outfile
n=word(dim,i+1);
outtitle=sprintf("Polynomial matrix multiplication over F_p \n (matrix dim = %s)",n);
set title outtitle	  
plot[:262144] \
	'matpol-result3.timing'  index i using 4:($5) w lp ls 1 title "Flint",\
       	'matpol-result3.timing'  index i using 4:($6) w lp ls 2 title "Mathemagix",\
       	'matpol-result3.timing'  index i using 4:($7) w lp ls 3 title "LinBox/FFLAS"
 }	

