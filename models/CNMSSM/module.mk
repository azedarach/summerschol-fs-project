DIR          := models/CNMSSM
MODNAME      := CNMSSM
SARAH_MODEL  := NMSSM

CNMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

CNMSSM_MK     := \
		$(DIR)/module.mk

CNMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

CNMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

CNMSSM_TWO_SCALE_MK := \
		$(CNMSSM_TWO_SCALE_SUSY_MK) \
		$(CNMSSM_TWO_SCALE_SOFT_MK)

CNMSSM_SLHA_INPUT := \


CNMSSM_GNUPLOT := \
		$(DIR)/CNMSSM_plot_rgflow.gnuplot \
		$(DIR)/CNMSSM_plot_spectrum.gnuplot

CNMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBCNMSSM_SRC :=
EXECNMSSM_SRC :=

LIBCNMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBCNMSSM_SRC += \
		$(DIR)/CNMSSM_info.cpp \
		$(DIR)/CNMSSM_slha_io.cpp \
		$(DIR)/CNMSSM_physical.cpp \
		$(DIR)/CNMSSM_utilities.cpp \
		$(DIR)/CNMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/CNMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/CNMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/CNMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/CNMSSM_two_scale_model.cpp \
		$(DIR)/CNMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/CNMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/CNMSSM_two_scale_susy_scale_constraint.cpp
EXECNMSSM_SRC += \
		$(DIR)/run_CNMSSM.cpp \
		$(DIR)/scan_CNMSSM.cpp
LIBCNMSSM_HDR += \
		$(DIR)/CNMSSM_convergence_tester.hpp \
		$(DIR)/CNMSSM_high_scale_constraint.hpp \
		$(DIR)/CNMSSM_info.hpp \
		$(DIR)/CNMSSM_initial_guesser.hpp \
		$(DIR)/CNMSSM_input_parameters.hpp \
		$(DIR)/CNMSSM_low_scale_constraint.hpp \
		$(DIR)/CNMSSM_model.hpp \
		$(DIR)/CNMSSM_physical.hpp \
		$(DIR)/CNMSSM_slha_io.hpp \
		$(DIR)/CNMSSM_spectrum_generator.hpp \
		$(DIR)/CNMSSM_susy_scale_constraint.hpp \
		$(DIR)/CNMSSM_utilities.hpp \
		$(DIR)/CNMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/CNMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/CNMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/CNMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/CNMSSM_two_scale_model.hpp \
		$(DIR)/CNMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/CNMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/CNMSSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(CNMSSM_TWO_SCALE_SUSY_MK)
-include $(CNMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CNMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CNMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

endif

# remove duplicates in case all algorithms are used
LIBCNMSSM_SRC := $(sort $(LIBCNMSSM_SRC))
EXECNMSSM_SRC := $(sort $(EXECNMSSM_SRC))

LIBCNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCNMSSM_SRC)))

EXECNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECNMSSM_SRC)))

LIBCNMSSM_DEP := \
		$(LIBCNMSSM_OBJ:.o=.d)

EXECNMSSM_DEP := \
		$(EXECNMSSM_OBJ:.o=.d)

LIBCNMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_CNMSSM_OBJ := $(DIR)/run_CNMSSM.o
RUN_CNMSSM_EXE := $(DIR)/run_CNMSSM.x

SCAN_CNMSSM_OBJ := $(DIR)/scan_CNMSSM.o
SCAN_CNMSSM_EXE := $(DIR)/scan_CNMSSM.x

METACODE_STAMP_CNMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CNMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCNMSSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNMSSM_SRC) $(CNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCNMSSM_HDR) $(CNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXECNMSSM_SRC) $(CNMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(CNMSSM_MK) $(CNMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CNMSSM_TWO_SCALE_MK) $(CNMSSM_INSTALL_DIR)
ifneq ($(CNMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CNMSSM_SLHA_INPUT) $(CNMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CNMSSM_GNUPLOT) $(CNMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBCNMSSM_DEP)
		-rm -f $(EXECNMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCNMSSM_OBJ)
		-rm -f $(EXECNMSSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBCNMSSM)
		-rm -f $(RUN_CNMSSM_EXE)
		-rm -f $(SCAN_CNMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(CNMSSM_TARBALL) \
		$(LIBCNMSSM_SRC) $(LIBCNMSSM_HDR) \
		$(EXECNMSSM_SRC) \
		$(CNMSSM_MK) $(CNMSSM_TWO_SCALE_MK) \
		$(CNMSSM_SLHA_INPUT) $(CNMSSM_GNUPLOT)

$(LIBCNMSSM_SRC) $(LIBCNMSSM_HDR) $(EXECNMSSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CNMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CNMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CNMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_CNMSSM)"
		@echo "Note: to regenerate CNMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CNMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CNMSSM):
		@true
endif

$(LIBCNMSSM_DEP) $(EXECNMSSM_DEP) $(LIBCNMSSM_OBJ) $(EXECNMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCNMSSM_DEP) $(EXECNMSSM_DEP) $(LIBCNMSSM_OBJ) $(EXECNMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBCNMSSM): $(LIBCNMSSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_CNMSSM_EXE): $(RUN_CNMSSM_OBJ) $(LIBCNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(SCAN_CNMSSM_EXE): $(SCAN_CNMSSM_OBJ) $(LIBCNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBCNMSSM_DEP) $(EXECNMSSM_DEP)
ALLSRC += $(LIBCNMSSM_SRC) $(EXECNMSSM_SRC)
ALLLIB += $(LIBCNMSSM)
ALLEXE += $(RUN_CNMSSM_EXE) $(SCAN_CNMSSM_EXE)
