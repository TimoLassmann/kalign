AUTOMAKE_OPTIONS = serial-tests
AM_COLOR_TESTS=always
EXTRA_DIST = run_io_test.sh \
data/a2m.good.1 \
data/a2m.good.2 \
data/afa.good.1 \
data/afa.good.2 \
data/afa.good.3 \
data/clustal.good.1 \
data/clustal.good.2


TESTS_ENVIRONMENT = \
      myprog_abspath='$(top_srcdir)'/src/; \
         export myprog_abspath;\
	      testdatafiledir='$(top_srcdir)'/dev/data/; \
         export testdatafiledir;
	$(VALGRIND)


TESTS = run_io_test.sh



CLEANFILES =  rwtest.clu  rwtest.fasta  rwtest.msf
