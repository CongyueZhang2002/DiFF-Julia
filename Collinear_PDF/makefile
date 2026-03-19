# one needs to compile this file to a shared library

# makefile in "Collinear PDF/"
PREFIX ?= $(HOME)/UCLA_QCD_fitter/installation
FFLAGS  = -fPIC -I$(PREFIX)/include
LDFLAGS = -shared -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib
LIBS    = -lLHAPDF

pdfs.so: pdfs.f90
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

clean:
	rm -f *.so
