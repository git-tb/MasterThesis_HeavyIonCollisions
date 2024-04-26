MAIN = "main"
TEX = "pdflatex"
BIBLIO = "bibtex"

build:
	$(TEX) $(MAIN)
	$(BIBLIO) $(MAIN)
	$(TEX) $(MAIN)
	$(TEX) $(MAIN)

clean:
	rm -f *.bbl *.blg *.aux *.dvi *.idx *.ilg *.ind *.log *.nav *.out *.snm *.xdv *.toc