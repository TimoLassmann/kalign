# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@


CLEANFILES = ihmm_seq_itest_read_test.fa test_model_file.txt test.lfa test_dna.lfa libhdf5gluetest_exp.h5 libhdf5gluetest.h5 model_at_*h5 test1.h5 test2.h5

# top_builddir = @top_builddir@


AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef

MYLIBDIRS = ../tldevel/libtldevel.a


MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS) -lpthread


#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_LIBS)

AM_CPPFLAGS=${MYINCDIRS} ${HDF5_CPPFLAGS}



#LOGO_SOURCES = motif_logo.c motif_logo.h

bin_PROGRAMS =  kalign


kalign_SOURCES = run_kalign.c \
global.h \
parameters.h \
parameters.c \
align_io.h \
align_io.c \
alignment_parameters.h \
alignment_parameters.c \
tree_building.h \
tree_building.c \
alignment.h \
alignment.c \
weave_alignment.h \
weave_alignment.c \
misc.h \
misc.c




# TESTS =  kalign_ITEST
# TESTS_ENVIRONMENT = \
# 	$(VALGRIND)

# check_PROGRAMS = kalign_ITEST
# kalign_ITEST_SOURCES = $(KALIGNSOURCE)
# kalign_ITEST_CPPFLAGS = $(AM_CPPFLAGS) -DITEST
# kalign_ITEST_LDADD = $(MYLIBDIRS)



