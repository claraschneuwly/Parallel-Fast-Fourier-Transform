CXX = g++
CFLAGS = -pthread -std=c++11 -Wall

SOURCES = fft_1.cpp parallel_fft_paper3.cpp parallel_fft_paper4.cpp fft.cpp
OBJECTS = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SOURCES))

BUILD_DIR = build

.PHONY: all parallel simple

all: parallel simple

main: $(BUILD_DIR)/main

parallel: parallel3 parallel4

parallel3: $(BUILD_DIR)/parallel3

parallel4: $(BUILD_DIR)/parallel4

simple: $(BUILD_DIR)/simple

$(BUILD_DIR)/main: $(BUILD_DIR)/fft.o 
	$(CXX) $(CFLAGS) -o $@ $^

$(BUILD_DIR)/parallel3: $(BUILD_DIR)/parallel_fft_paper3.o
	$(CXX) $(CFLAGS) -o $@ $^

$(BUILD_DIR)/parallel4: $(BUILD_DIR)/parallel_fft_paper4.o
	$(CXX) $(CFLAGS) -o $@ $^

$(BUILD_DIR)/simple: $(BUILD_DIR)/fft_1.o
	$(CXX) $(CFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: %.cpp fft.hpp | $(BUILD_DIR)
	$(CXX) $(CFLAGS) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
