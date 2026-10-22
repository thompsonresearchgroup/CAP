# === Paths and Settings ===
mqcroot = $(HOME)/mqcPack
FC = gfortran
TARGET = noCAPs2withCAPs
SRC = noCAPs2withCAPs.f03

# === Compilation Flags ===
OPTFLAGS ?= -O3
FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp
STRICTFLAGS = -Wall -Wextra -Wpedantic -Wconversion -Wimplicit-interface -fcheck=all
INCLUDE = -I$(mqcroot)/GNU/mod

# === Libraries ===
MKLLIB ?= $(if $(MKLROOT),$(MKLROOT)/lib/intel64,$(HOME)/miniconda3/envs/science/lib)
LIBS = $(mqcroot)/GNU/lib/libmqc.a -L$(MKLLIB) -Wl,-rpath,$(MKLLIB) -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl

# === Default Build Target ===
all: $(TARGET)

.PHONY: all strict clean test test-python

tests/read_cap: tests/read_cap.f03
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ $< $(LIBS)

test: all tests/read_cap
	python3 -m unittest discover -s tests -p 'test_*.py' -v

test-python:
	PYTHONPATH=tests python3 -m unittest test_lebedev test_smooth_voronoi test_resonance test_wavepacket -v

strict:
	$(FC) $(FCFLAGS) $(STRICTFLAGS) $(INCLUDE) -o $(TARGET) $(SRC) $(LIBS)

$(TARGET): $(SRC) makefile $(mqcroot)/GNU/lib/libmqc.a
	$(FC) $(FCFLAGS) $(INCLUDE) $(OPTFLAGS) -o $@ $(SRC) $(LIBS)

# === Clean Target ===
clean:
	rm -f $(TARGET) *.o *.mod
