
MPICXX = mpicxx

TEMP_DIR_DEFAULT = /tmp
BASE_DIR_DEFAULT = ${PWD}

check-%:
	@ if test -z ${${*}}; then \
		echo "Environment variable $* not set"; \
		exit 1; \
	fi

test-%:
	@ if test -z ${${*}}; then \
		echo "Environment variable $* not set. Using default ${*}=${${*}_DEFAULT}"; \
		export ${*}=${${*}_DEFAULT}; \
	fi

BASE_DIR ?= $(BASE_DIR_DEFAULT)
TEMP_DIR ?= $(TEMP_DIR_DEFAULT)

NETCDF = ${NETCDF_DIR}
COMBBLAS = ${COMBBLAS_DIR}
PATOHBIN = ${PATOH_BIN}
METISBIN = ${METIS_BIN}
BASEDIR = ${BASE_DIR}
TEMPDIR = ${TEMP_DIR}

SRC_DIR = src
BIN_DIR = bin
OPT_FLAGS = -O3
DEBUG_FLAGS = -g
CXX_FLAGS = -std=c++11 -D_64BITOFFSET

INCL1_FLAGS = -I$(NETCDF)/include
INCL2_FLAGS = -I$(COMBBLAS)
CXX_INCL = $(INCL1_FLAGS) $(INCL2_FLAGS)

LIBS1_FLAGS = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf
LIBS2_FLAGS = -L$(COMBBLAS) -lCommGridlib -lMPITypelib -lMemoryPoollib

OBJ1 = $(SRC_DIR)/netcdf_utils.o $(SRC_DIR)/extract_max_level_cell.o
OBJ2 = $(SRC_DIR)/APowers.o
OBJ3 = $(SRC_DIR)/netcdf_utils.o					\
			 $(SRC_DIR)/mpas_ordering.o					\
			 $(SRC_DIR)/mpas_order.o						\
			 $(SRC_DIR)/mpas_ordering_utils.o		\
			 $(SRC_DIR)/mpas_ordering_xyzsort.o	\
			 $(SRC_DIR)/mpas_ordering_random.o	\
			 $(SRC_DIR)/mpas_ordering_morton.o	\
			 $(SRC_DIR)/mpas_ordering_hilbert.o	\
			 $(SRC_DIR)/mpas_ordering_peano.o
DEPS_HPP = $(SRC_DIR)/netcdf_utils.h

BIN1 = depths
BIN2 = apowers
BIN3 = mpasorder
MAIN = hampton
MAIN_METIS = hampton_metis

REQUIRED = check-NETCDF_DIR check-COMBBLAS_DIR check-PATOH_BIN
REQUIRED_METIS = check-NETCDF_DIR check-METIS_BIN
OPTIONAL = test-BASE_DIR test-TEMP_DIR
TESTALL = $(REQUIRED) $(OPTIONAL)
TESTALL_METIS = $(REQUIRED_METIS) $(OPTIONAL)

default: CXX_FLAGS+=$(OPT_FLAGS)
default: $(TESTALL_METIS) $(BIN1) $(BIN3) $(MAIN_METIS)

all: CXX_FLAGS+=$(OPT_FLAGS)
all: $(TESTALL) $(BIN1) $(BIN2) $(BIN3) $(MAIN)

debug: CXX_FLAGS+=$(DEBUG_FLAGS)
debug: $(TESTALL) $(BIN1) $(BIN2) $(BIN3) $(MAIN)

clean:
	rm -rf $(OBJ1) $(OBJ2) $(OBJ3) $(BIN_DIR)

$(BIN1): $(OBJ1)
	mkdir -p $(BIN_DIR)
	$(MPICXX) -o $(BIN_DIR)/$@ $^ $(CXX_FLAGS) $(LIBS1_FLAGS)

$(BIN2): $(OBJ2)
	mkdir -p $(BIN_DIR)
	$(MPICXX) -o $(BIN_DIR)/$@ $^ $(CXX_FLAGS) $(LIBS2_FLAGS)

$(BIN3): $(OBJ3)
	mkdir -p $(BIN_DIR)
	$(MPICXX) -o $(BIN_DIR)/$@ $^ $(CXX_FLAGS) $(LIBS1_FLAGS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS_HPP)
	$(MPICXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)


$(MAIN_METIS): $(SRC_DIR)/workflow.metis.base
	mkdir -p $(BIN_DIR)
	cat $< | sed s,ENV_METIS_BIN,$(METISBIN),g \
		| sed s,ENV_BASE_DIR,$(BASEDIR),g \
		| sed s,ENV_TEMP_DIR,$(TEMPDIR),g > $(BIN_DIR)/$(MAIN_METIS)
	chmod +x $(BIN_DIR)/$(MAIN_METIS)

$(MAIN): $(SRC_DIR)/workflow.base
	mkdir -p $(BIN_DIR)
	cat $< | sed s,ENV_PATOH_BIN,$(PATOHBIN),g \
		| sed s,ENV_BASE_DIR,$(BASEDIR),g \
		| sed s,ENV_TEMP_DIR,$(TEMPDIR),g > $(BIN_DIR)/$(MAIN)
	chmod +x $(BIN_DIR)/$(MAIN)
