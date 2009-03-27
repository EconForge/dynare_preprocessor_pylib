macroprocessor.pdf: macroprocessor.tex old-design.pdf new-design.pdf
	pdflatex macroprocessor
	pdflatex macroprocessor

clean:
	rm -f macroprocessor.pdf *.toc *.aux *.log *.nav *.snm *.vrb *.out *~
