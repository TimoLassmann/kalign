# Build the library in the hand subdirectory first.

ACLOCAL_AMFLAGS = -I m4

@VALGRIND_CHECK_RULES@

AM_CFLAGS =  @SIMD_FLAGS@ @OPENMP_CFLAGS@
CLEANFILES = ihmm_seq_itest_read_test.fa test_model_file.txt test.lfa test_dna.lfa libhdf5gluetest_exp.h5 libhdf5gluetest.h5 model_at_*h5 test1.h5 test2.h5

# top_builddir = @top_builddir@
AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef

bin_PROGRAMS = timescorealn kalign


kalign_SOURCES = run_kalign.c \
tldevel.h \
tldevel.c \
esl_stopwatch.h \
esl_stopwatch.c \
rng.h \
rng.c \
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


check_PROGRAMS = edist alphabet misc_test bpm_test rwaln rstdin_test
TESTS = edist alphabet  misc_test bpm_test rstdin_test

# TESTS_ENVIRONMENT = $(VALGRIND)

rwaln_SOURCES = \
tldevel.h \
tldevel.c \
esl_stopwatch.h \
esl_stopwatch.c \
rwalign.c \
alphabet.h \
alphabet.c \
msa.h
rwaln_CPPFLAGS = $(AM_CPPFLAGS) -DRWALIGN_TEST

timescorealn_SOURCES = \
tldevel.h \
tldevel.c \
esl_stopwatch.h \
esl_stopwatch.c \
time_and_score_aln.c \
msa.h \
parameters.h \
parameters.c \
alphabet.h \
alphabet.c \
align_io.h \
rwalign.c \
misc.h \
misc.c

bpm_test_SOURCES = \
tldevel.h \
tldevel.c \
esl_stopwatch.h \
esl_stopwatch.c \
rng.h \
rng.c \
bpm.h \
bpm.c \
alphabet.h \
alphabet.c
bpm_test_CPPFLAGS = $(AM_CPPFLAGS) -DBPM_UTEST

misc_test_SOURCES = \
tldevel.h \
tldevel.c \
rng.h \
rng.c \
misc.h \
misc.c \
alphabet.h \
alphabet.c

misc_test_CPPFLAGS = $(AM_CPPFLAGS) -DITEST_MISC

alphabet_SOURCES = \
tldevel.h \
tldevel.c \
alphabet.h \
alphabet.c
alphabet_CPPFLAGS = $(AM_CPPFLAGS) -DITEST_ALPHABET

edist_SOURCES = \
tldevel.h \
tldevel.c \
esl_stopwatch.h \
esl_stopwatch.c \
rng.h \
rng.c \
euclidean_dist.h \
euclidean_dist.c

edist_CPPFLAGS = $(AM_CPPFLAGS) -DITEST_EDIST

rstdin_test_SOURCES = \
read_from_stdin_test.c \
tldevel.h \
tldevel.c
