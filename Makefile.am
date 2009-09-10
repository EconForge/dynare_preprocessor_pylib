if HAVE_PDFLATEX
pdf-local: macroprocessor.pdf
endif

SRC = macroprocessor.tex old-design.pdf new-design.pdf

EXTRA_DIST = $(SRC)

macroprocessor.pdf: $(SRC)
	$(PDFLATEX) macroprocessor
	$(PDFLATEX) macroprocessor

clean-local:
	rm -f macroprocessor.pdf *.toc *.aux *.log *.nav *.snm *.vrb *.out *~
