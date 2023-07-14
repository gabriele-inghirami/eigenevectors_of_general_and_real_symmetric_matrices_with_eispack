# Compiler options
COMPILER = gfortran
FLAGS = -O3 -march=native

# Source files
SOURCE_LIBRARY_FILES = eispack_real_symm_matrix.f90 eispack_general_matrix.f90
SOURCE_EXEC_FILE = test_eispack.f90
OBJECT_LIBRARY_FILES = $(SOURCE_LIBRARY_FILES:.f90=.o)
OBJECT_EXEC_FILE = $(SOURCE_EXEC_FILE:.f90=.o)

# Names of the targets
LIB_TARGET = eispack.a
EXE_TARGET = test.exe

# Symbolic targets not associated with files, but with actions
.PHONY: all clean

all: $(LIB_TARGET) $(EXE_TARGET)

# To build the static library
$(LIB_TARGET): $(OBJECT_LIBRARY_FILES)
	ar rcs $@ $^

# To build the executable for testing purposes
$(EXE_TARGET): $(OBJECT_EXEC_FILE)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LIB_TARGET)

# Rule to compile source files
%.o: %.f90
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -f $(OBJECT_LIBRARY_FILES) $(OBJECT_EXEC_FILE) $(LIB_TARGET) $(EXE_TARGET)
