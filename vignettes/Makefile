all=PSSeg.pdf

%.tex: %.Rnw
	$(R_HOME)/bin/Rscript -e "library(knitr); knit('$<')"

%.pdf: %.tex jointSeg.bib
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<
	$(R_HOME)/bin/Rscript -e "tools::compactPDF('$@', gs_quality='ebook')"

clean: 
	rm -f $<.aux
	rm -f $<.bbl
	rm -f $<.blg
	rm -f $<.log
	rm -f $<.out
	rm -f $<.toc
