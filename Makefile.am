EXTRA_DIST = version.sh autogen.sh VERSION

ACLOCAL_AMFLAGS = -I m4

# Build the library in the hand subdirectory first
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@


# top_builddir = @top_builddir@

SUBDIRS = tldevel src


count:
	line=$$(find . \( -name "*.[c,h]" -o -name "*.am" -o -name "*.org" -o -name "*.sh" \)  -type f -print0  | xargs -0 wc -l | grep total); echo	$$line
