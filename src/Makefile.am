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


RANDOMKIT_FILES = distributions.h \
distributions.c \
randomkit.h \
randomkit.c


ADJUSTEDRANDINDEXSOURCE = adjusted_rand_index.c adjusted_rand_index.h

MOTIFREFINE = motif_refinement.h motif_refinement.c

ISEQSOURCE = ihmm_seq.h ihmm_seq.c

FASTHMMSOURCE = fast_hmm_param.c fast_hmm_param.h  \
fast_hmm_param_test_functions.h \
fast_hmm_param_test_functions.c

KALIGNSOURCE = kalign.h kalign.c

MODELSOURCE = model.h model.c model_io.c global.h

BEAMSOURCE = beam_sample.h beam_sample.c

THREADSOURCE = thread_data.c thread_data.h

HDF5GLUESOURCE = hdf5_glue.c hdf5_glue.h #save_random_state.c

SCORESOURCE = run_score.h run_score.c

EMITSOURCE = run_emit.c emit_random.h

CALIBRATESOURCE = calibrate.c calibrate.h

CONVERSION = hmm_conversion.h hmm_conversion.c

FINITEHMM = finite_hmm.c finite_hmm.h

DIJSTRA_SOURCES = dijkstra.c dijkstra.h

LOGO_SOURCES = motif_logo.c motif_logo.h

bin_PROGRAMS =  wims_kmer sim_seq wims_model wims_plot wims_score wims_lcs wims_iclu wims_emit cal tfbs_bench composition_bench wims_select wims_ari

if ENABLE_CAIRO
bin_PROGRAMS +=  cairo_test wims_pplot wims_seqmotif
AM_CPPFLAGS += $(CAIRO_CFLAGS)
cairo_test_LDADD = $(CAIRO_LIBS) $(MYLIBDIRS)
wims_pplot_LDADD = $(CAIRO_LIBS) $(MYLIBDIRS)
wims_seqmotif_LDADD = $(CAIRO_LIBS) $(MYLIBDIRS)
endif


wims_ari_SOURCES = model_ari_comparison.c $(MODELSOURCE) $(ISEQSOURCE) $(ADJUSTEDRANDINDEXSOURCE) $(HDF5GLUESOURCE) $(FINITEHMM) $(RANDOMKIT_FILES)

wims_seqmotif_SOURCES = draw_sequence_logo.c $(LOGO_SOURCES) $(SCORESOURCE) $(MODELSOURCE) $(THREADSOURCE) $(ISEQSOURCE) $(HDF5GLUESOURCE)  $(FINITEHMM) $(RANDOMKIT_FILES)

wims_pplot_SOURCES = pplot.c   $(LOGO_SOURCES)  $(DIJSTRA_SOURCES) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM)  $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES)  $(HDF5GLUESOURCE)

cairo_test_SOURCES = cairo_test.c

cal_SOURCES = $(CALIBRATESOURCE) $(ISEQSOURCE) $(FINITEHMM) $(EMITSOURCE) $(SCORESOURCE) $(HDF5GLUESOURCE) $(RANDOMKIT_FILES) $(THREADSOURCE)
cal_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)


wims_plot_SOURCES = \
plot_model.c  $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM)  $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES)  $(HDF5GLUESOURCE)
#wims_plot_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)

wims_model_SOURCES = \
build_model.c init_seq_label.h init_seq_label.c  $(BEAMSOURCE) $(SCORESOURCE) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE) $(THREADSOURCE)
#wims_model_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)

wims_select_SOURCES = \
select_model.c  $(SCORESOURCE) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(FASTHMMSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE) $(THREADSOURCE)

wims_kmer_SOURCES = \
kmer_enrichment.c $(KALIGNSOURCE) $(MOTIFREFINE) $(SCORESOURCE) $(EMITSOURCE) $(CONVERSION) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE) $(THREADSOURCE)


wims_score_SOURCES = \
score_sequences.c $(CONVERSION)  $(SCORESOURCE) $(ISEQSOURCE) $(FINITEHMM) $(MODELSOURCE) $(FASTHMMSOURCE)  $(RANDOMKIT_FILES) $(HDF5GLUESOURCE) $(THREADSOURCE)
#wims_score_LDADD = ${MYLIBDIRS} $(HDF5_LDFLAGS)  $(HDF5_LIBS)

wims_lcs_SOURCES = $(ISEQSOURCE) $(HDF5GLUESOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) lcs.c lcs.h
#wims_lcs_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


wims_iclu_SOURCES = $(SCORESOURCE) $(MODELSOURCE) $(MOTIFREFINE) $(THREADSOURCE) $(ISEQSOURCE) $(HDF5GLUESOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) label_suffix_array.h label_suffix_array.c  infoclust.c infoclust.h
#wims_iclu_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


wims_emit_SOURCES = $(ISEQSOURCE) $(MODELSOURCE) $(HDF5GLUESOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) $(EMITSOURCE) emit_random.c

#wims_emit_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


sim_seq_SOURCES = \
outdir.h \
outdir.c \
matrix_io.h \
matrix_io.c \
sim_seq.c \
make_dot_file.h \
make_dot_file.c

sim_seq_LDADD =  ${MYLIBDIRS}

tfbs_bench_SOURCES = $(RANDOMKIT_FILES) \
outdir.h \
outdir.c \
matrix_io.h \
matrix_io.c \
benchmark_seq.h \
benchmark_seq.c \
create_TFBS_benchmark_datasets.c

tfbs_bench_LDADD =  ${MYLIBDIRS}



composition_bench_SOURCES = $(RANDOMKIT_FILES) \
outdir.h \
outdir.c \
matrix_io.h \
matrix_io.c \
benchmark_seq.h \
benchmark_seq.c \
create_composition_benchmark_dataset.c

composition_bench_LDADD =  ${MYLIBDIRS}
#build_ihmm_SOURCES = \
model_nuc_seq.c \
matrix_io.h \
matrix_io.c \
ihmm_io.h \
ihmm_io.c \
outdir.h \
outdir.c \
make_dot_file.h \
make_dot_file.c

#build_ihmm_LDADD =  ${MYLIBDIRS}  libihmm.a



#noinst_LIBRARIES = libihmm.a

# AUTOMAKE_OPTIONS = subdir-objects
#libihmm_a_SOURCES = $(RANDOMKIT_FILES) ihmm.h ihmm.c

#libihmm_a_LIBADD  =  ${MYLIBDIRS}


TESTS =  kalign_ITEST dijkstra_ITEST fast_hmm_param_ITEST model_ITEST 	iseq_ITEST beam_ITEST dtest_ITEST   hdf5glue_ITEST ari_ITEST

TESTS_ENVIRONMENT = \
	$(VALGRIND)

check_PROGRAMS = kalign_ITEST dijkstra_ITEST pst_ITEST fast_hmm_param_ITEST model_ITEST iseq_ITEST beam_ITEST dtest_ITEST  hdf5glue_ITEST ari_ITEST

#libihmm_ITEST_SOURCES = $(libihmm_a_SOURCES)
#libihmm_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
#libihmm_ITEST_LDADD = $(libihmm_a_LIBADD)


dijkstra_ITEST_SOURCES = $(DIJSTRA_SOURCES)
dijkstra_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITEST
dijkstra_ITEST_LDADD = $(libihmm_a_LIBADD)


iseq_ITEST_SOURCES = $(ISEQSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE)
iseq_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTSEQ
iseq_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

dtest_ITEST_SOURCES  = dtest.c $(RANDOMKIT_FILES)
dtest_ITEST_CPPFLAGS = $(AM_CPPFLAGS)
dtest_ITEST_LDADD = $(MYLIBDIRS)


fast_hmm_param_ITEST_SOURCES = $(FASTHMMSOURCE) $(RANDOMKIT_FILES)
fast_hmm_param_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTFASTPARAM
fast_hmm_param_ITEST_LDADD = $(MYLIBDIRS)


hdf5glue_ITEST_SOURCES = $(HDF5GLUESOURCE)
hdf5glue_ITEST_CPPFLAGS  = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS)  -DITEST
hdf5glue_ITEST_LDADD  = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


model_ITEST_SOURCES =  $(MODELSOURCE) $(ISEQSOURCE) $(RANDOMKIT_FILES) $(HDF5GLUESOURCE)
model_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTMODEL
model_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


beam_ITEST_SOURCES = $(BEAMSOURCE) $(THREADSOURCE) $(ISEQSOURCE) $(MODELSOURCE)  $(HDF5GLUESOURCE)  $(FASTHMMSOURCE) $(RANDOMKIT_FILES)  $(CONVERSION) $(FINITEHMM)
beam_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTBEAM
beam_ITEST_LDADD = $(libihmm_a_LIBADD) $(HDF5_LDFLAGS)  $(HDF5_LIBS)

pst_ITEST_SOURCES = $(ISEQSOURCE) $(HDF5GLUESOURCE) $(RANDOMKIT_FILES)  $(FINITEHMM) lcs.c lcs.h
pst_ITEST_CPPFLAGS = $(AM_CPPFLAGS)   -DITESTLST
pst_ITEST_LDADD = $(MYLIBDIRS) $(HDF5_LDFLAGS)  $(HDF5_LIBS)


ari_ITEST_SOURCES = $(ADJUSTEDRANDINDEXSOURCE)
ari_ITEST_CPPFLAGS = $(AM_CPPFLAGS) -DITEST
ari_ITEST_LDADD = $(MYLIBDIRS)

kalign_ITEST_SOURCES = $(KALIGNSOURCE)
kalign_ITEST_CPPFLAGS = $(AM_CPPFLAGS) -DITEST
kalign_ITEST_LDADD = $(MYLIBDIRS)



