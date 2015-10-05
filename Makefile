
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
		export ${*}=${${*}}_DEFAULT; \
	fi

NETCDF = ${NETCDF_DIR}
COMBBLAS = ${COMBBLAS_DIR}
PATOHBIN = ${PATOH_BIN}
BASEDIR = ${BASE_DIR}
TEMPDIR = ${TEMP_DIR}

CXX = mpicxx

SRC_DIR = src
BIN_DIR = bin
OPT_FLAGS = -O3
DEBUG_FLAGS = -g
CXX_FLAGS = -std=c++11

INCL1_FLAGS = -I$(NETCDF)/include
INCL2_FLAGS = -I$(COMBBLAS)
CXX_INCL = $(INCL1_FLAGS) $(INCL2_FLAGS)

LIBS1_FLAGS = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf
LIBS2_FLAGS = -L$(COMBBLAS) -lCommGridlib -lMPITypelib -lMemoryPoollib

OBJ1 = $(SRC_DIR)/netcdf_utils.o $(SRC_DIR)/extract_max_level_cell.o
OBJ2 = $(SRC_DIR)/APowers.o
DEPS_HPP = $(SRC_DIR)/netcdf_utils.h

BIN1 = depths
BIN2 = apowers
MAIN = hampton

REQUIRED = check-NETCDF_DIR check-COMBBLAS_DIR check-PATOH_BIN
OPTIONAL = test-BASE_DIR test-TEMP_DIR
TESTALL = $(REQUIRED) $(OPTIONAL)

all: CXX_FLAGS+=$(OPT_FLAGS)
all: $(TESTALL) $(BIN1) $(BIN2) $(MAIN)

debug: CXX_FLAGS+=$(DEBUG_FLAGS)
debug: $(TESTALL) $(BIN1) $(BIN2) $(MAIN)

clean:
	rm -f $(OBJ1) $(OBJ2) $(BIN_DIR)/$(BIN1) $(BIN_DIR)/$(BIN2) $(BIN_DIR)/$(MAIN)
	rmdir $(BIN_DIR)

$(BIN1): $(OBJ1)
	mkdir -p $(BIN_DIR)
	$(CXX) -o $(BIN_DIR)/$@ $^ $(CXX_FLAGS) $(LIBS1_FLAGS)

$(BIN2): $(OBJ2)
	mkdir -p $(BIN_DIR)
	$(CXX) -o $(BIN_DIR)/$@ $^ $(CXX_FLAGS) $(LIBS2_FLAGS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS_HPP)
	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)


$(MAIN): $(SRC_DIR)/workflow.base
	mkdir -p $(BIN_DIR)
	cat $< | sed s,ENV_PATOH_BIN,$(PATOHBIN),g \
		| sed s,ENV_BASE_DIR,$(BASEDIR),g \
		| sed s,ENV_TEMP_DIR,$(TEMPDIR),g > $(BIN_DIR)/$(MAIN)
	chmod +x $(BIN_DIR)/$(MAIN)
