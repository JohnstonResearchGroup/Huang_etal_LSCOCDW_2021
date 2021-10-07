LAPACK = -llapack
F90 = mpif90

RIXS_OBJ = bandstructure.o phonon_dispersion.o

SPEC_OBJ = bandstructure.o phonon_dispersion.o

bandstructure.o:
	$(F90) -c bandstructure.f90

phonon_dispersion.o:
	$(F90) -c phonon_dispersion.f90

compute_spectral_functions: $(SPEC_OBJ)
	$(F90)  $(OPT) -o compute_spec compute_spectral_functions.f90 $(SPEC_OBJ) $(LAPACK)

compute_rixs:  $(RIXS_OBJ)
	$(F90)  $(OPT) -o compute_rixs compute_rixs.f90 $(RIXS_OBJ) $(LAPACK)

clean:
	rm *.mod 
	rm *.o
	rm compute_rixs
	rm compute_spec
