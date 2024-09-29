# Variables
CXX = g++
CXXFLAGS = -Wall -O2
SRC_DIR = src
TARGET = clinvar-ancestry
SRC = $(SRC_DIR)/clinvar-ancestry.cpp

# Rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)

.PHONY: all clean
