latex:
	pdflatex article.tex
	bibtex   article
	pdflatex article.tex
	pdflatex article.tex
	evince   article.pdf &

plot:
	cd BENCH && \
	gnuplot plotting-bench-mul.plot && \
	for f in *.eps; \
	do \
		ps2pdf $$f; \
		pdfcrop $${f%%.*}.pdf $${f%%.*}.pdf; \
	done && \
	cd mod-comp && \
	gnuplot plot_timings.plot && \
	for f in *.eps; \
	do \
		ps2pdf $$f; \
		pdfcrop $${f%%.*}.pdf $${f%%.*}.pdf; \
	done
