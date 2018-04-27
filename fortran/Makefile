# Please select your favorite Fortran compiler here
FC=gfortran

starsylv: example.f90 starsylv.f90
	$(FC) -o starsylv example.f90 starsylv.f90 -O3 -lblas

clean:
	rm -f starsylv
