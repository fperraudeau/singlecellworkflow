RMDFILE = $(wildcard *.Rmd)
HTML = $(patsubst %.Rmd, %.html, $(RMDFILE))
SUPPORT = createData.R
all: $(HTML)

%.html: %.Rmd $(SUPPORT)
	Rscript -e "rmarkdown::render('$<')"

%.R: %.Rmd $(SUPPORT)
	Rscript -e "require(knitr); purl('$<')"

clean:
	rm -rf *.md
