if ENABLE_DOC
pdf-local: macroprocessor.pdf
endif

SRC = macroprocessor.tex new-design.pdf

EXTRA_DIST = $(SRC)

macroprocessor.pdf: $(SRC)
	$(PDFLATEX) macroprocessor
	$(PDFLATEX) macroprocessor

clean-local:
	rm -f macroprocessor.pdf *.toc *.aux *.log *.nav *.snm *.vrb *.out *~
