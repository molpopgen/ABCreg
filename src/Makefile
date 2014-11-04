CC=cc
CXX=c++
#DEBUG=-g
#DEBUG=-DNDEBUG
DEBUG=
CFLAGS=-O3 -Wall -W $(DEBUG) -I.
CXXFLAGS=-std=c++11 -O3 -Wall -W $(DEBUG) -I.
DEBUGFLAGS=-g -Wall -W -ansi -pedantic
LIBS=-lgsl -lgslcblas -lz

TARGETS=reg.o params.o process_options.o read_prior.o transformations.o abc_algebra.o generate_posterior.o
TARGETS_DEBUG=	reg_debug.o params_debug.o process_options_debug.o read_prior_debug.o transformations_debug.o abc_algebra_debug.o generate_posterior_debug.o

REGDEP=abc_algebra.hpp generate_posterior.hpp params.hpp process_options.hpp read_prior.hpp transformations.hpp
PROCESSOPTIONSDEP=process_options.hpp params.hpp
READPRIORDEP=read_prior.hpp transformations.hpp params.hpp
TRANSFORMATIONSDEP=transformations.hpp params.hpp
ABCALGEBRADEP=abc_algebra.hpp 
GENERATEPOSTERIORDEP=abc_algebra.hpp generate_posterior.hpp params.hpp
PARAMSDEP=params.hpp
all: $(TARGETS)
	$(CXX) $(CXXFLAGS) -o reg $(TARGETS) $(LIBS)

clean:
	rm -f *.o reg

debug: $(TARGETS_DEBUG)
	$(CXX) -o reg $(TARGETS_DEBUG) $(LIBS)

dist:
	mkdir ABCreg
#	cp *.cc *.hpp *.c *.h Makefile ABCreg
	cp *.cc *.hpp Makefile getmap.R exampleABC.sh abcregdoc.tex literature.bib ABCreg
	tar czf ABCreg.tar.gz ABCreg
	rm -rf ABCreg

reg.o: reg.cc $(REGDEP)

params.o: params.cc $(PARAMSDEP)

process_options.o: process_options.cc $(PROCESSOPTIONSDEP)

read_prior.o: read_prior.cc  $(PROCESSOPTIONSDEP)

transformations.o: transformations.cc $(TRANSFORMATIONSDEP)

abc_algebra.o: abc_algebra.cc $(ABCALGEBRADEP)

generate_posterior.o: generate_posterior.cc $(GENERATEPOSTERIORDEP)

reg_debug.o: reg.cc $(REGDEP)
	$(CXX) -I. -c reg.cc -o reg_debug.o $(DEBUGFLAGS)

params_debug.o: params.cc $(PARAMSDEP)
	$(CXX) -I. -c params.cc -o params_debug.o $(DEBUGFLAGS)

process_options_debug.o: process_options.cc $(PROCESSOPTIONSDEP)
	$(CXX) -I. -c process_options.cc -o process_options_debug.o $(DEBUGFLAGS)

read_prior_debug.o: read_prior.cc $(PROCESSOPTIONSDEP)
	$(CXX) -I. -c read_prior.cc -o read_prior_debug.o $(DEBUGFLAGS)

transformations_debug.o: transformations.cc $(TRANSFORMATIONSDEP)
	$(CXX) -I. -c transformations.cc -o transformations_debug.o $(DEBUGFLAGS)

abc_algebra_debug.o: abc_algebra.cc $(ABCALGEBRADEP)
	$(CXX) -I. -c abc_algebra.cc -o abc_algebra_debug.o $(DEBUGFLAGS)

generate_posterior_debug.o:  generate_posterior.cc $(GENERATEPOSTERIORDEP)
	$(CXX) -I. -c generate_posterior.cc -o generate_posterior_debug.o $(DEBUGFLAGS)
