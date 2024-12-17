## TFHE LIBRARY
TFHE_FOLDER_PATH = /home/ab235386/code/tfhe/tfhe32/mylibs
TFHE_LIB_PATH = $(TFHE_FOLDER_PATH)/lib
TFHE_INCLUDE_PATH = $(TFHE_FOLDER_PATH)/include

CC := g++ -std=c++14 -no-pie

# Our code paths
DIR = /home/ab235386/code/github/tfhe_id_bs
IDIR := $(DIR)/include
SDIR := $(DIR)/src
ODIR := $(DIR)/obj
BDIR := $(DIR)/bin

INCS:= -I${TFHE_INCLUDE_PATH} -I${IDIR}

LIBS:= -L${TFHE_LIB_PATH} -ltfhe-spqlios-fma

# Executable objects
_OBJ_EXEC = tfhe_bootstrapping.o main.o

OBJ_EXEC = $(patsubst %,$(ODIR)/%,$(_OBJ_EXEC))

.PHONY := all clean delete

all: $(BDIR)/exe

${ODIR}/%.o: ${SDIR}/%.cpp
	$(CC) -c -o $@ $< ${INCS}

$(BDIR)/exe: $(OBJ_EXEC) $(OBJ)
	$(CC) $^ -o $@ ${INCS} $(LIBS)

clean:
	rm -f $(ODIR)/*.o $(BDIR)/*

# delete:
# 	rm -f $(DDIR)/*
