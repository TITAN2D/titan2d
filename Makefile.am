#
# automake template for the top-most source directory
#
# $Id: Makefile.am 235 2012-03-29 07:04:35Z dkumar $
#

SUBDIRS = src scripts gui doc tests share

install-data-local:
	@if [ "$(TITAN2D_PORTABLE)" = "yes" ]; then \
		echo "Moving dependencies"; \
		rsync -a $(TITAN2D_DEP) $(libdir); \
		echo "Clean-up dependencies"; \
		rm /home/mikola/titan_wsp/titan2d_bld/gccportinstall/lib/titan2d_dep/lib/*.a; \
	fi
