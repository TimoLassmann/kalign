# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@


ACLOCAL_AMFLAGS = -I m4




AM_CFLAGS =  @SIMD_FLAGS@ @OPENMP_CFLAGS@ 
CLEANFILES = ihmm_seq_itest_read_test.fa test_model_file.txt test.lfa test_dna.lfa libhdf5gluetest_exp.h5 libhdf5gluetest.h5 model_at_*h5 test1.h5 test2.h5

# top_builddir = @top_builddir@


AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef



MYLIBDIRS = ../tldevel/libtldevel.a


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS) -lpthread -lm 


#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)

AM_CPPFLAGS=${MYINCDIRS} ${HDF5_CPPFLAGS}



#LOGO_SOURCES = motif_logo.c motif_logo.h

bin_PROGRAMS =  kalign timescorealn 


kalign_SOURCES = run_kalign.c \
global.h \
msa.h \
parameters.h \
parameters.c \
aln_data.h \
aln_data.c \
alphabet.h \
alphabet.c \
alignment_parameters.h \
alignment_parameters.c \
pick_anchor.h \
pick_anchor.c \
bisectingKmeans.h \
bisectingKmeans.c \
sequence_distance.h \
sequence_distance.c \
euclidean_dist.h \
euclidean_dist.c \
alignment.h \
alignment.c \
weave_alignment.h \
weave_alignment.c \
rwalign.c \
misc.h \
misc.c \
bpm.h \
bpm.c


check_PROGRAMS =  phmm get_param edist alphabet heap_test misc_test bpm_test
TESTS =  phmm get_param edist alphabet heap_test misc_test bpm_test
TESTS_ENVIRONMENT = $(VALGRIND)

timescorealn_SOURCES = \
time_and_score_aln.c \
msa.h \
parameters.h \
parameters.c \
alphabet.h \
alphabet.c \
palign_io.h \
rwalign.c \
misc.h \
misc.c

bpm_test_SOURCES = \
bpm.h \
bpm.c \
alphabet.h \
alphabet.c
bpm_test_CPPFLAGS = $(AM_CPPFLAGS) -DBPM_UTEST

misc_test_SOURCES = \
misc.h \
misc.c \
alphabet.h \
alphabet.c

misc_test_CPPFLAGS = $(AM_CPPFLAGS) -DITEST_MISC

heap_test_SOURCES = \
heap.h \
heap.c
heap_test_CPPFLAGS = $(AM_CPPFLAGS) -DITEST

alphabet_SOURCES = \
alphabet.h \
alphabet.c
alphabet_CPPFLAGS = $(AM_CPPFLAGS) -DITEST

edist_SOURCES = \
euclidean_dist.h \
euclidean_dist.c

edist_CPPFLAGS = $(AM_CPPFLAGS) -DITEST


get_param_SOURCES = run_get_param.c \
msa.h \
rwalign.c \
parameters.h \
parameters.c \
global.h \
misc.h \
misc.c \
alphabet.h \
alphabet.c \
euclidean_dist.h \
euclidean_dist.c \
kmeans.h \
kmeans.c

phmm_SOURCES = phmm.h phmm.c
phmm_LDADD = $(MYLIBDIRS)
phmm_CPPFLAGS = $(AM_CPPFLAGS) -DITEST
