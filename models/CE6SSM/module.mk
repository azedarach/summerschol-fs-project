DIR          := models/CE6SSM
MODNAME      := CE6SSM
SARAH_MODEL  := E6SSM

CE6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

CE6SSM_MK     := \
		$(DIR)/module.mk

CE6SSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

CE6SSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

CE6SSM_TWO_SCALE_MK := \
		$(CE6SSM_TWO_SCALE_SUSY_MK) \
		$(CE6SSM_TWO_SCALE_SOFT_MK)

CE6SSM_SLHA_INPUT := \


CE6SSM_GNUPLOT := \
		$(DIR)/CE6SSM_plot_rgflow.gnuplot \
		$(DIR)/CE6SSM_plot_spectrum.gnuplot

CE6SSM_TARBALL := \
		$(MODNAME).tar.gz

LIBCE6SSM_SRC :=
EXECE6SSM_SRC :=

LIBCE6SSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBCE6SSM_SRC += \
		$(DIR)/CE6SSM_info.cpp \
		$(DIR)/CE6SSM_slha_io.cpp \
		$(DIR)/CE6SSM_physical.cpp \
		$(DIR)/CE6SSM_utilities.cpp \
		$(DIR)/CE6SSM_two_scale_convergence_tester.cpp \
		$(DIR)/CE6SSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/CE6SSM_two_scale_initial_guesser.cpp \
		$(DIR)/CE6SSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/CE6SSM_two_scale_model.cpp \
		$(DIR)/CE6SSM_two_scale_susy_parameters.cpp \
		$(DIR)/CE6SSM_two_scale_soft_parameters.cpp \
		$(DIR)/CE6SSM_two_scale_susy_scale_constraint.cpp
EXECE6SSM_SRC += \
		$(DIR)/run_CE6SSM.cpp \
		$(DIR)/scan_CE6SSM.cpp
LIBCE6SSM_HDR += \
		$(DIR)/CE6SSM_convergence_tester.hpp \
		$(DIR)/CE6SSM_high_scale_constraint.hpp \
		$(DIR)/CE6SSM_info.hpp \
		$(DIR)/CE6SSM_initial_guesser.hpp \
		$(DIR)/CE6SSM_input_parameters.hpp \
		$(DIR)/CE6SSM_low_scale_constraint.hpp \
		$(DIR)/CE6SSM_model.hpp \
		$(DIR)/CE6SSM_physical.hpp \
		$(DIR)/CE6SSM_slha_io.hpp \
		$(DIR)/CE6SSM_spectrum_generator.hpp \
		$(DIR)/CE6SSM_susy_scale_constraint.hpp \
		$(DIR)/CE6SSM_utilities.hpp \
		$(DIR)/CE6SSM_two_scale_convergence_tester.hpp \
		$(DIR)/CE6SSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/CE6SSM_two_scale_initial_guesser.hpp \
		$(DIR)/CE6SSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/CE6SSM_two_scale_model.hpp \
		$(DIR)/CE6SSM_two_scale_soft_parameters.hpp \
		$(DIR)/CE6SSM_two_scale_susy_parameters.hpp \
		$(DIR)/CE6SSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(CE6SSM_TWO_SCALE_SUSY_MK)
-include $(CE6SSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CE6SSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CE6SSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBCE6SSM_SRC := $(sort $(LIBCE6SSM_SRC))
EXECE6SSM_SRC := $(sort $(EXECE6SSM_SRC))

LIBCE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCE6SSM_SRC)))

EXECE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECE6SSM_SRC)))

LIBCE6SSM_DEP := \
		$(LIBCE6SSM_OBJ:.o=.d)

EXECE6SSM_DEP := \
		$(EXECE6SSM_OBJ:.o=.d)

LIBCE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_CE6SSM_OBJ := $(DIR)/run_CE6SSM.o
RUN_CE6SSM_EXE := $(DIR)/run_CE6SSM.x

SCAN_CE6SSM_OBJ := $(DIR)/scan_CE6SSM.o
SCAN_CE6SSM_EXE := $(DIR)/scan_CE6SSM.x

METACODE_STAMP_CE6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CE6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCE6SSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCE6SSM_SRC) $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCE6SSM_HDR) $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXECE6SSM_SRC) $(CE6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(CE6SSM_MK) $(CE6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CE6SSM_TWO_SCALE_MK) $(CE6SSM_INSTALL_DIR)
ifneq ($(CE6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CE6SSM_SLHA_INPUT) $(CE6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CE6SSM_GNUPLOT) $(CE6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBCE6SSM_DEP)
		-rm -f $(EXECE6SSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCE6SSM_OBJ)
		-rm -f $(EXECE6SSM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBCE6SSM)
		-rm -f $(RUN_CE6SSM_EXE)
		-rm -f $(SCAN_CE6SSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(CE6SSM_TARBALL) \
		$(LIBCE6SSM_SRC) $(LIBCE6SSM_HDR) \
		$(EXECE6SSM_SRC) \
		$(CE6SSM_MK) $(CE6SSM_TWO_SCALE_MK) \
		$(CE6SSM_SLHA_INPUT) $(CE6SSM_GNUPLOT)

$(LIBCE6SSM_SRC) $(LIBCE6SSM_HDR) $(EXECE6SSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CE6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CE6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CE6SSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_CE6SSM)"
		@echo "Note: to regenerate CE6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CE6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CE6SSM):
		@true
endif

$(LIBCE6SSM_DEP) $(EXECE6SSM_DEP) $(LIBCE6SSM_OBJ) $(EXECE6SSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCE6SSM_DEP) $(EXECE6SSM_DEP) $(LIBCE6SSM_OBJ) $(EXECE6SSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBCE6SSM): $(LIBCE6SSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_CE6SSM_EXE): $(RUN_CE6SSM_OBJ) $(LIBCE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

$(SCAN_CE6SSM_EXE): $(SCAN_CE6SSM_OBJ) $(LIBCE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(THREADLIBS)

ALLDEP += $(LIBCE6SSM_DEP) $(EXECE6SSM_DEP)
ALLSRC += $(LIBCE6SSM_SRC) $(EXECE6SSM_SRC)
ALLLIB += $(LIBCE6SSM)
ALLEXE += $(RUN_CE6SSM_EXE) $(SCAN_CE6SSM_EXE)
