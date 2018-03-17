BENCHMARK_DIR=benchmarks
BIN_DIR=bin
INC_DIR=include
OBJ_DIR=obj
OBJ_DBG_DIR=obj/debug
SRC_DIR=cpp
TOOL_DIR=tool

MATRICES=resources/control_matrices.txt
DEMO_BIN=$(BIN_DIR)/demo
TB_BIN=$(BIN_DIR)/testbench
BENCHMAKR_TOOL=$(TOOL_DIR)/benchmark
PLOT_TOOL=$(TOOL_DIR)/plot

_DEMO_OBJS=blockdecoder.o ctrlmat.o demo.o
DEMO_OBJS=$(patsubst %.o, $(OBJ_DBG_DIR)/%.o, $(_DEMO_OBJS))

_TB_OBJS=blockdecoder.o ctrlmat.o testbench.o
TB_OBJS=$(patsubst %.o, $(OBJ_DIR)/%.o, $(_TB_OBJS))

DEPS=$(wildcard $(INC_DIR)/*)

CPPFLAGS=-O3 -g -Wall

all: $(DEMO_BIN) $(TB_BIN)

$(DEMO_BIN): $(DEMO_OBJS)
	g++ -o $@ $^

$(TB_BIN): $(TB_OBJS)
	g++ -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(DEPS)
	g++ -c -o $@ $< -I$(INC_DIR) $(CPPFLAGS) -DNDEBUG

$(OBJ_DBG_DIR)/%.o: $(SRC_DIR)/%.cc $(DEPS)
	g++ -c -o $@ $< -I$(INC_DIR) $(CPPFLAGS)

.PHONY: clean benchmark demo

benchmark: $(BENCHMARK_TOOL) $(TB_BIN)
	$(BENCHMAKR_TOOL) $(MATRICES) $(BENCHMARK_DIR)

clean:
	@rm -f $(OBJ_DIR)/*.o $(OBJ_DBG_DIR)/*.o $(BIN_DIR)/*
	@rm -rf $(BENCHMARK_DIR)/*

demo: $(DEMO_BIN)
	$(DEMO_BIN)
