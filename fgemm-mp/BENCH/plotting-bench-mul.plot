set encoding iso_8859_1 

set xlabel "Entry bitsize" 
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

dim="32 64 256 512 1024"

do for [i=0:4] {
outfile=sprintf('bench-mul-%d.eps',i);
set output outfile
n=word(dim,i+1);
outtitle=sprintf("(matrix dim = %s)",n);
set title outtitle	  
plot[:262144] \
	'fgemm-result.timing'  index i using 3:($5) w lp ls 1 title "Flint",\
       	'fgemm-result.timing'  index i using 3:($6) w lp ls 2 title "Mathemagix",\
       	'fgemm-result.timing'  index i using 3:($4) w lp ls 3 title "FFLAS"
 }	

do for [i=0:3] {
outfile=sprintf('bench-mul-full-%d.eps',i);
set output outfile
n=word(dim,i+1);
outtitle=sprintf("Integer matrix multiplication \n (matrix dim = %s)",n);
set title outtitle	  
plot[:262144] \
	'fgemm-result.timing'  index i using 3:($5) w lp ls 1 title "Flint",\
       	'fgemm-result.timing'  index i using 3:($6) w lp ls 2 title "Mathemagix",\
       	'fgemm-result.timing'  index i using 3:($4) w lp ls 3 title "FFLAS",\
	'fgemm-result.timing'  index i using 3:($7) w lp ls 4 title "LinBox"
 }	



# dim="2048 8192 32768"
# do for [i=0:2] {

# outfile=sprintf('bench-mul-by-bitsize-%d.eps',i);
# set output outfile
# n=word(dim,i+1);
# outtitle=sprintf("Integer matrix multiplication \n (bitsize = %s)",n);
# set title outtitle	  
# plot \
# 	'fgemm-result-by-bitsize.timing'  index i using 2:($5) w lp ls 1 title "Flint",\
#        	'fgemm-result-by-bitsize.timing'  index i using 2:($6) w lp ls 2 title "Mathemagix",\
#        	'fgemm-result-by-bitsize.timing'  index i using 2:($4) w lp ls 3 title "FFLAS"
#  }	

