#
# automake template for src/adapt directory
#
# $Id: Makefile.am,v 1.1.1.1 2003/08/13 19:26:11 sorokine Exp $
#

bin_PROGRAMS = lhsbed lhsvol lhstitanstats gentitanteststats

lhsbed_SOURCES = lhsbed.C lhslib.C lhslib.h

lhsvol_SOURCES = lhsvol.C lhslib.C lhslib.h

lhstitanstats_SOURCES = lhstitanstats.C lhslib.C lhslib.h

gentitanteststats_SOURCES = gentitanteststats.C lhslib.C lhslib.h

lhsbed_LDADD = ../useful/libuseful.a

lhsvol_LDADD = ../useful/libuseful.a

lhstitanstats_LDADD = ../useful/libuseful.a

gentitanteststats_LDADD = ../useful/libuseful.a

EXTRA_DIST = pbslhs dist-stats.pl  plot_stats.pl lhs.readme
