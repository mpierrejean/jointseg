doc=PSSeg

pdf: ${doc}.Rnw
	Rscript -e "library(knitr); knit('${doc}.Rnw')"
	pdflatex ${doc}
	bibtex ${doc}
	pdflatex ${doc}
	pdflatex ${doc}
	Rscript -e "tools::compactPDF('${doc}.pdf', gs_quality='ebook')"

clean: 
	rm -f PSSeg.aux
	rm -f PSSeg.bbl
	rm -f PSSeg.blg
	rm -f PSSeg.log
	rm -f PSSeg.out
	rm -f PSSeg.toc
