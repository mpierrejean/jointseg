doc=PSSeg

pdf: ${doc}.Rnw
	Rscript -e "library(knitr); knit('${doc}.Rnw')"
	pdflatex ${doc}
	bibtex ${doc}
	pdflatex ${doc}
	pdflatex ${doc}
