## Set which compiler to use by defining CCCOM:
##GNU GCC compiler
#CCCOM=g++ -m64 -std=c++11 
##Clang compiler (good to use on Mac OS)
#CCCOM=clang++ -std=c++1y
##Intel C++ compiler (good to use with Intel MKL if available)
CCCOM=mpic++ -std=c++11 -g 
#########


## Flags to give the compiler for "release mode"



#LIBFLAGS = -larmadillo
LIBSPECTRA = -I/home/keaideni/WORK/Spectral/spectra/include/ -I/home/keaideni/WORK/Spectral/eigen-git-mirror/ -fopenmp





obj=main.o SingleSub.o 
main:$(obj)
	$(CCCOM) -o main $(obj)  $(LIBSPECTRA)
main.o:main.cpp SuperEnergy.h Wave.h 
	$(CCCOM) -c main.cpp -O2 $(LIBSPECTRA)
SingleSub.o:SingleSub.cpp SingleSub.h Parameter.h
	$(CCCOM) -c SingleSub.cpp -O2 $(LIBSPECTRA)
.PHONY:clean
clean:
	rm -f main $(obj)















