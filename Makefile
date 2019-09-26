all: thesis.pdf thesis.bbl

clean:
	rm -f *.dvi *.aux *.log *.bbl *.pdf 
	rm -f presentation/*.dvi presentation/*.aux *.pdf presentation/*.log

thesis.bbl:
	xelatex thesis.tex
	bibtex thesis

thesis.pdf: thesis.tex thesis.bbl thesis.bib graphics/*.eps unrthesis.cls *.c
	xelatex thesis.tex
	bibtex thesis
	xelatex thesis.tex

presentation.pdf:
	xelatex -output-directory=presentation presentation/presentation.tex 
