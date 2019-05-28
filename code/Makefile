P=programe_name
OBJECTS=
CFLAGS += -g -Wall -std=gnu11 
LDLIBS = -lhdf5 -lhdf5_hl -lgsl -lgslcblas -lm 

$(P): $(OBJECTS)

all : libabaa.a adu_test adu_vector_test dye_test illumination_test \
ccd_test aba_test ratio_test ratio_test2

test : adu_test_run dye_test_run illumination_test_run ccd_test_run \
aba_test_run adu_vector_test_run ratio_test2_run

libabaa.a : adu.o adu_vector.o illumination.o dye.o ccd.o aba.o ratio.o \
fura.o
	ar cr libabaa.a adu.o adu_vector.o illumination.o dye.o ccd.o aba.o ratio.o \
fura.o

adu_test_objects=adu_test.o
adu_test : $(adu_test_objects)
	cc $(adu_test_objects) $(LDLIBS) -o adu_test

adu_test_run : adu_test
	./adu_test

adu_vector_test_objects=adu_vector_test.o
adu_vector_test : $(adu_vector_test_objects)
	cc $(adu_vector_test_objects) $(LDLIBS) -o adu_vector_test

adu_vector_test_run : adu_vector_test
	cd .. && ./code/adu_vector_test > toto && gnuplot -persist code/adu_vector_test_plot_loading_curve.gp

dye_test_objects=dye_test.o
dye_test : $(dye_test_objects)
	cc $(dye_test_objects) $(LDLIBS) -o dye_test

dye_test_run : dye_test
	cd .. && ./code/dye_test

illumination_test_objects=illumination_test.o
illumination_test : $(illumination_test_objects)
	cc $(illumination_test_objects) $(LDLIBS) -o illumination_test

illumination_test_run : illumination_test
	cd .. && ./code/illumination_test

ccd_test_objects=ccd_test.o
ccd_test : $(ccd_test_objects)
	cc $(ccd_test_objects) $(LDLIBS) -o ccd_test

ccd_test_run : ccd_test
	cd .. && ./code/ccd_test

aba_test_objects=aba_test.o
aba_test : $(aba_test_objects)
	cc $(aba_test_objects) $(LDLIBS) -o aba_test

aba_test_run : aba_test
	cd .. && ./code/aba_test | head -n 30

ratio_test_objects=ratio_test.o
ratio_test : $(ratio_test_objects)
	cc $(ratio_test_objects) $(LDLIBS) -o ratio_test

ratio_test2.o : abaa.h

ratio_test2_objects=ratio_test2.o
ratio_test2 : $(ratio_test2_objects) libabaa.a
	cc $(ratio_test2_objects) libabaa.a $(LDLIBS) -o ratio_test2

ratio_test2_run : ratio_test2
	cd .. && ./code/ratio_test2 > titi && gnuplot -persist code/ratio_test_fig.gp

fit_ratiometric : fit_ratiometric.o libabaa.a
	cc fit_ratiometric.o libabaa.a $(LDLIBS) -o fit_ratiometric

fura_concentration : fura_concentration.o libabaa.a
	cc fura_concentration.o libabaa.a $(LDLIBS) -o fura_concentration

fura_during_stim : fura_during_stim.o libabaa.a
	cc fura_during_stim.o libabaa.a $(LDLIBS) -o fura_during_stim

aba_ratio : aba_ratio.o libabaa.a
	cc aba_ratio.o libabaa.a $(LDLIBS) -o aba_ratio

.PHONY : clean
clean :
	rm -f libabaa.a adu.o adu_vector.o illumination.o dye.o ccd.o aba.o ratio.o \
	fura.o \
	$(adu_test_objects) adu_test \
	$(adu_vector_test_objects) adu_vector_test \
	$(dye_test_objects) dye_test \
	$(illumination_test_objects) illumination_test \
	$(ccd_test_objects) ccd_test \
	$(aba_test_objects) aba_test \
	$(ratio_test_objects) ratio_test \
	$(ratio_test2_objects) ratio_test2 \
	fit_ratiometric.o fit_ratiometric \
	fura_concentration.o fura_concentration \
	fura_during_stim.o fura_during_stim \
	aba_ratio.o aba_ratio
