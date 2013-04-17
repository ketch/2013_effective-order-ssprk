BN=Effective_order_SSPRK_methods

## SINUM version ##
FV=ExSSPRK_effective_order-siam

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

sinum:

	latex $(FV).tex
	bibtex ${FV}
	latex ${FV}.tex
	latex ${FV}.tex
	dvips ${FV}.dvi -o
	ps2pdf ${FV}.ps

sinum-quick:

	latex ${FV}.tex
	dvips ${FV}.dvi -o
	ps2pdf ${FV}.ps

sinum-clean:

	rm -f ${FV}.toc
	rm -f ${FV}.dvi
	rm -f ${FV}.lof
	rm -f ${FV}.lot
	rm -f ${FV}.aux
	rm -f ${FV}.bbl
	rm -f ${FV}.blg
	rm -f ${FV}.idx
	rm -f ${FV}.ilg
	rm -f ${FV}.log
	rm -f ${FV}.ind
	rm -f ${FV}.out
	rm -f ${FV}.ps




