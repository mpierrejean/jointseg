doc=PSSeg

pdf: ${doc}.Rnw
	$(R_HOME)/bin/Rscript -e "library(knitr); knit('${doc}.Rnw')"
	pdflatex ${doc}
	bibtex ${doc}
	pdflatex ${doc}
	pdflatex ${doc}
