preprocessor.pdf: preprocessor.tex expr.png expr-sharing.png matrices.png overview.png
	pdflatex preprocessor
	pdflatex preprocessor

clean:
	rm -f preprocessor.{pdf,toc,aux,log,nav,snm,vrb,out} *~
