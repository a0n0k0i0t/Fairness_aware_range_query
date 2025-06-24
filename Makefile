ifeq ($(PLATFORM),windows)
    CXX = g++
    EXE_EXT = .exe
    BOOST_INCLUDE = C:/local/boost_1_82_0  
else
    CXX = clang++
    EXE_EXT =
    BOOST_INCLUDE = /opt/homebrew/include  
endif

SRC = skyline_MD.cpp
BIN = skyline_MD$(EXE_EXT)

CXXFLAGS = -std=c++17 -O2 -Wall -I$(BOOST_INCLUDE)

all: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

run: all
ifeq ($(PLATFORM),windows)
	$(BIN)
else
	./$(BIN)
endif

.PHONY: all run clean

clean:
ifeq ($(PLATFORM),windows)
	del /Q skyline_1D.exe skyline_MD.exe 2>nul || exit 0
else
	rm -f skyline_1D skyline_MD
endif
