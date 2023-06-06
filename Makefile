CXX = g++
CFLAGS = -pthread -std=c++11 -Wall

SOURCES = fft_1.cpp parallel_fft.cpp
OBJECTS = fft_1.o parallel_fft.o 

BUILD_DIR = build

.PHONY: all parallel simple

all: parallel simple

parallel: $(BUILD_DIR)/parallel

simple: $(BUILD_DIR)/simple

$(BUILD_DIR)/parallel: $(addprefix $(BUILD_DIR)/,$(OBJECTS))
	$(CXX) $(CFLAGS) -o $@ $(filter-out $(BUILD_DIR)/fft_1.o,$^)

$(BUILD_DIR)/simple: $(BUILD_DIR)/fft_1.o
	$(CXX) $(CFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CFLAGS) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
