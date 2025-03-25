CC=gcc
CFLAGS=-Wall -std=c++11 -g -fPIC -fopenmp -O3 $(FLAGS)

LIBS=-lstdc++ -lz -ldl -pthread -lbz2 -llzma -lm
INCLUDE=./htslib

SRC=.
BUILD=build


KMC_API_DIR = ../KMC/kmc_api
KMC_DUMP_DIR = ../KMC/kmc_dump

KMC_UTIL_OBJS = \
$(KMC_DUMP_DIR)/nc_utils.o

KMC_API_OBJS = \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o

HTS_LIB=-L./htslib ./htslib/libhts.a

_OBJS=anchor.o
OBJS = $(patsubst %, $(BUILD)/%, $(_OBJS))

DEPENDS := $(patsubst %.o, %.d, $(OBJS))

all: dirs run_anchor

run_anchor: build/anchor.o
	$(CC) $(KMC_API_OBJS) $(CFLAGS) $< -o $@ -lz $(HTS_LIB) $(LIBS)

-include $(DEPENDS)

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@ -I$(INCLUDE)

.PHONY: dirs
dirs: $(BUILD)/ #$(BIN)/

$(BUILD)/:
	mkdir -p $@

#$BIN)/:
#	mkdir -p $@

clean:
	rm -rf $(BUILD)
