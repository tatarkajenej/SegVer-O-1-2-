ALL=build/nice build/ugly
COMMON=segver.cpp select.cpp

all: $(ALL)

build/%: %.cpp $(COMMON)
	g++ $< -o"$@" -fopenmp -lblas -lgivaro

.PHONY: run
run: $(ALL)
	python3 ./nice.py
	python3 ./ugly.py

