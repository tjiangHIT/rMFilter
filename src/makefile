
BIN_PATH=bin
BUILD_PATH=build
BIN_INDEX=rMFilter-indexer
BIN_ALIGN=rMFilter-aligner
RM=rm -rf
# General compiler flags
CC=g++
CFLAGS=-O3 -c -Wall
DFLAGS=-g -Wall
LDFLAGS= -lz -lm -lpthread
DEGUG=debug

# ***************SOURCES**************
KSW_SRC := ksw.c 
ALIGN_SRC := Aln.cpp Del_local_aln.cpp mk_hash.cpp LCtrl_option.cpp 
INDEX_SRC := Index.cpp mk_hash.cpp readfl.cpp Ctrl_option.cpp

ALIGN_OBJECTS := $(ALIGN_SRC:%.cpp=%.o) $(KSW_SRC:%.c=%.o)
INDEX_OBJECTS := $(INDEX_SRC:%.cpp=%.o) $(KSW_SRC:%.c=%.o)


# ***************COMPIILE**************
# all:$(BIN_PATH)/$(BIN_INDEX) $(BIN_PATH)/$(BIN_ALIGN)
all:$(BIN_ALIGN) $(BIN_INDEX) 
	@echo "Making symlink: $(BIN_INDEX) $(BIN_ALIGN)"

$(BIN_INDEX):$(INDEX_OBJECTS)
	@echo "Linking: $@"
	$(CC) -o $@ $^ $(LDFLAGS)

$(BIN_ALIGN):$(ALIGN_OBJECTS)
	@echo "Linking: $@"
	$(CC) -o $@ $^ $(LDFLAGS)

.cpp.o:
	@echo "Compiling: $@"
	$(CC) $(CFLAGS) $< -o $@
.c.o:
	@echo "Compiling: $@"
	$(CC) $(CFLAGS) $< -o $@

# .PHONY: release
# release:dirs

# .PHONY: dirs
# dirs:
# 	@echo "Creating directories"
# 	@mkdir $(BUILD_PATH)
# 	@mkdir $(BIN_PATH)


.PHONY: clean
clean:
	@echo "Deleting rMFilter"
	@$(RM) *.o $(BIN_INDEX) $(BIN_ALIGN)
	@$(RM) $(DEGUG)

.PHONY: debug
debug:
	@$(CC) $(DFLAGS) $(ALIGN_SRC) $(KSW_SRC) $(LDFLAGS) -o $(DEGUG)

# Installs to the set path
.PHONY: install
install:
	@echo "Installing to /usr/local/bin"
	install $(BIN_INDEX) $(BIN_ALIGN) /usr/local/bin/

# Uninstalls the program
.PHONY: uninstall
uninstall:
	@echo "Removing /usr/local/bin/rMFilter"
	@$(RM) /usr/local/bin/$(BIN_INDEX)
	@$(RM) /usr/local/bin/$(BIN_ALIGN)