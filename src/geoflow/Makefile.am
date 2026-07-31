#
# automake template for src/geoflow directory
#
# $Id: Makefile.am 218 2009-10-22 14:38:24Z dkumar $
#

noinst_LIBRARIES = libgeoflow.a libgeoflow_solver.a
#predict.f eigen.f getcoef.f correct.f drag_model.f sgn.f 
libgeoflow_a_SOURCES = element_weight.C \
   rnr.h slopes.C move_data.C setup_geoflow.C stats.C flux_srcs.C

libgeoflow_solver_a_SOURCES = edge_states.C get_coef_and_eigen.C integrators.C

libgeoflow_solver_a_CXXFLAGS = @OPENMP_CXXFLAGS@

#integrators.$(OBJEXT):CXXFLAGS +=-I/opt/intel/advisor_xe/include -qopt-report=5 -qopt-report-phase=vec



