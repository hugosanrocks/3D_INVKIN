#======================#
#     include file     #
#======================#

HOST=$(shell hostname)

# matches SEISCOPE
ifeq ($(findstring SEISCOPE,$(HOST)), SEISCOPE)
include make.SEISCOPE.inc
endif

#=====================#
#     directories     #
#=====================#

BASE_DIR              = $(DGCRACK2D)
# COMMON FILES
COMMON_INC_DIR        = $(BASE_DIR)/COMMON/include
COMMON_UTIL_DIR       = $(BASE_DIR)/COMMON/util
# BUTTERWORTH FILTER FILE
FILTER_DIR            = ./../src
# TEST FILES
TEST_BASE_DIR         = ./../test
TEST_INC_DIR          = ./../include
TEST_BIN_DIR          = ./../bin

VPATH = $(COMMON_UTIL_DIR):$(FILTER_DIR):$(TEST_BASE_DIR)

FFLAGS    = -I. -I $(TEST_INC_DIR) -I $(COMMON_INC_DIR) $(FOPTION)

#=====================#
#     objects         #
#=====================#

FILTER_OBJ = init_butterworth.o filtfilt_butterworth.o

COMMON_UTIL_OBJ = flush_c.o memory_stat.o tilt.o

BASE_OBJ = main.o

ALL_OBJ = $(COMMON_UTIL_OBJ) $(FILTER_OBJ) $(BASE_OBJ)

.SUFFIXES: .o .f90
.f90.o:
	$(FC) -c $(FFLAGS) $<

Test_bin: $(ALL_OBJ)
	$(FC) $(FFLAGS) -o $(TEST_BIN_DIR)/$@ $(ALL_OBJ) $(LIBS)

objs:	$(ALL_OBJ)

$(ALL_OBJ): $(COMMON_INC_DIR)/paral.h $(COMMON_INC_DIR)/param.h \
$(COMMON_INC_DIR)/global.h $(COMMON_INC_DIR)/precis.h $(TEST_INC_DIR)/type_def.h

clean:
	rm -f *.o *.mod *~ ../bin/Test_bin
