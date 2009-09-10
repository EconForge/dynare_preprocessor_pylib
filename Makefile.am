if HAVE_PDFLATEX
all-local: preprocessor.pdf
endif

SRC = preprocessor.tex expr.png expr-sharing.png matrices.png overview.png

EXTRA_DIST = $(SRC)

preprocessor.pdf: $(SRC)
	$(PDFLATEX) preprocessor
	$(PDFLATEX) preprocessor

clean-local:
	rm -f *.pdf *.toc *.aux *.log *.nav *.snm *.vrb *.out *~
