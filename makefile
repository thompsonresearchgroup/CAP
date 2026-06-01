# === Paths and Settings ===
mqcroot = $(HOME)/mqcPack
FC = gfortran
TARGET = noCAPs2withCAPs
SRC = noCAPs2withCAPs.f03

# === Compilation Flags ===
FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp
INCLUDE = -I$(mqcroot)/GNU/mod

# === Libraries ===
LIBS = $(mqcroot)/GNU/lib/libmqc.a -L${MKLROOT}/lib/intel64 -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl

# === Default Build Target ===
all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FCFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# === Clean Target ===
clean:
	rm -f $(TARGET) *.o *.mod

