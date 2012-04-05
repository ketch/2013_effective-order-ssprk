BN=Effective_order_SSPRK_methods

all:
	latex $(BN).tex
	bibtex ${BN}
	latex ${BN}.tex
	latex ${BN}.tex
	dvips ${BN}.dvi -o
	ps2pdf ${BN}.ps

quick:
	latex ${BN}.tex
	dvips ${BN}.dvi -o
	ps2pdf ${BN}.ps

clean:
	rm -f ${BN}.toc
	rm -f ${BN}.dvi
	rm -f ${BN}.lof
	rm -f ${BN}.lot
	rm -f ${BN}.aux
	rm -f ${BN}.bbl
	rm -f ${BN}.blg
	rm -f ${BN}.idx
	rm -f ${BN}.ilg
	rm -f ${BN}.log
	rm -f ${BN}.ind
	rm -f ${BN}.out
	rm -f ${BN}.ps

