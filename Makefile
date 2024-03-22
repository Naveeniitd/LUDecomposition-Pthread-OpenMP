# Makefile for compiling a1_openmp.cpp and a1_pthread.cpp

# Compiler variables
CC_OPENMP = g++-13
CC_PTHREAD = g++
CC_A1 = g++
# Compiler flags
FLAGS_OPENMP = -fopenmp
FLAGS_PTHREAD = -lpthread -std=c++11
FLAGS_A1 =  -std=c++11
# Output binaries
OUT_OPENMP = a1_openmp
OUT_PTHREAD = a1_pthread
OUT_A1 = a1
# Source files
SRC_OPENMP = a1_openmp.cpp
SRC_PTHREAD = a1_pthread.cpp
SRC_A1 = a1.cpp

# Default target
all: $(OUT_OPENMP) $(OUT_PTHREAD) $(OUT_A1)
# Rule for a1
$(OUT_A1): $(SRC_A1)
	$(CC_A1) -o $@ $< $(FLAGS_A1)
# Rule for a1_openmp
$(OUT_OPENMP): $(SRC_OPENMP)
	$(CC_OPENMP) $(FLAGS_OPENMP) $< -o $@

# Rule for a1_pthread
$(OUT_PTHREAD): $(SRC_PTHREAD)
	$(CC_PTHREAD) -o $@ $< $(FLAGS_PTHREAD)

# Clean
clean:
	rm -f $(OUT_OPENMP) $(OUT_PTHREAD) $(OUT_A1)

.PHONY: all clean

