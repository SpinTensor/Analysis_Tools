include Make.inc

all: fourier.x \
	autocorr.x \
	crosscorr.x \
	histogram.x \
	histogram2d.x \
	histogram3d.x

ALLOCDIR = ./external/Fortran_Reallocate/
LIBALLOCA = $(ALLOCDIR)/libredealloc.a
LIBALLOC = -L$(ALLOCDIR)/ -lredealloc
INCALLOC = -I$(ALLOCDIR)/

CMDOPTDIR = ./external/cmd_option_parser/
LIBCMDOPTA = $(CMDOPTDIR)/libcmdopt.a
LIBCMDOPT = -L$(CMDOPTDIR)/ -lcmdopt
INCCMDOPT = -I$(CMDOPTDIR)/


$(LIBALLOCA):
	(cp Make.inc $(ALLOCDIR))
	(cd $(ALLOCDIR); make)

$(LIBCMDOPTA):
	(cp Make.inc $(CMDOPTDIR))
	(cd $(CMDOPTDIR); make)

kinds.o: kinds.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) -c $<

fourier.x: fourier.o kinds.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(LIBALLOC) $(LIBCMDOPT) $(LIBFFTW)
fourier.o: fourier.F90 kinds.o $(LIBALLOCA) $(LIBCMDOPTA)
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(INCALLOC) $(INCCMDOPT) $(INCFFTW) -c $<

autocorr.x: autocorr.o kinds.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(LIBALLOC) $(LIBCMDOPT) $(LIBFFTW)
autocorr.o: autocorr.F90 kinds.o $(LIBALLOCA) $(LIBCMDOPTA)
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(INCALLOC) $(INCCMDOPT) $(INCFFTW) -c $<

crosscorr.x: crosscorr.o kinds.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(LIBALLOC) $(LIBCMDOPT) $(LIBFFTW)
crosscorr.o: crosscorr.F90 kinds.o $(LIBALLOCA) $(LIBCMDOPTA)
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(INCALLOC) $(INCCMDOPT) $(INCFFTW) -c $<

histogram.x: histogram.o kinds.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(LIBALLOC) $(LIBCMDOPT)
histogram.o: histogram.F90 kinds.o $(LIBALLOCA) $(LIBCMDOPTA)
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(INCALLOC) $(INCCMDOPT) -c $<

histogram2d.x: histogram2d.o kinds.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(LIBALLOC) $(LIBCMDOPT)
histogram2d.o: histogram2d.F90 kinds.o $(LIBALLOCA) $(LIBCMDOPTA)
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(INCALLOC) $(INCCMDOPT) -c $<

histogram3d.x: histogram3d.o kinds.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(LIBALLOC) $(LIBCMDOPT)
histogram3d.o: histogram3d.F90 kinds.o $(LIBALLOCA) $(LIBCMDOPTA)
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(INCALLOC) $(INCCMDOPT) -c $<

.PHONY: clean distclean

clean:
	rm -f *.o *.mod 
	(cd $(ALLOCDIR); make clean)
	(cd $(CMDOPTDIR); make clean)

distclean: clean
	rm -f fourier.x
	rm -f histogram.x
	rm -f histogram2d.x
	rm -f histogram3d.x
	(cd $(ALLOCDIR); make distclean)
	(cd $(CMDOPTDIR); make distclean)
