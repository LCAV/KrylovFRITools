TEX = pdflatex -shell-escape -interaction=nonstopmode -file-line-error

all: web doc

web: fri.html fri.webtex.html fri.wiki fri.gh.md

doc: fri.pdf fri.docx fri.odt

view: 
	open fri.html

fri.html: fri.md fri.bib
	pandoc fri.md -Ss --mathjax --filter pandoc-citeproc  -o fri.html

fri.webtex.html: fri.md fri.bib
	pandoc fri.md -Ss --webtex --filter pandoc-citeproc  -o fri.webtex.html

fri.docx: fri.md fri.bib
	pandoc fri.md -Ss --filter pandoc-citeproc  -o fri.docx

fri.odt: fri.md fri.bib
	pandoc fri.md -Ss --filter pandoc-citeproc  -o fri.odt

fri.wiki: fri.md fri.bib
	pandoc fri.md -Ss --filter pandoc-citeproc -t mediawiki -o fri.wiki

fri.gh.md: fri.md fri.bib
	pandoc fri.md -Ss --filter pandoc-citeproc -t markdown_github -o fri.gh.md

fri.pdf: fri.tex 
	$(TEX) fri.tex ; rm fri.aux fri.log fri.out

fri.tex : fri.md fri.bib
	pandoc fri.md -Ss --filter pandoc-citeproc  -o fri.tex 

clean:
	rm 	fri.html fri.wiki fri.docx fri.pdf	fri.gh.md fri.tex fri.odt fri.webtex.html