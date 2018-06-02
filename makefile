CXX = g++
CXXFLAGS = -lharppi -lfftw3 -lfftw3_omp -fopenmp 
CXXOPTS = -march=native -O3

build: cic gadgetReader particle pma main.cpp
	$(CXX) $(CXXFLAGS) $(CXXOPTS) -o cipress main.cpp obj/gadgetReader.o obj/particle.o obj/pma.o obj/cic.o
	
cic: source/cic.cpp
	mkdir -p obj
	$(CXX) $(CXXOPTS) -c -o obj/cic.o source/cic.cpp
	
gadgetReader: source/gadgetReader.cpp
	mkdir -p obj
	$(CXX) $(CXXOPTS) -c -o obj/gadgetReader.o source/gadgetReader.cpp
	
particle: source/particle.cpp
	mkdir -p obj
	$(CXX) $(CXXOPTS) -c -o obj/particle.o source/particle.cpp
	
pma: source/pma.cpp
	mkdir -p obj
	$(CXX) $(CXXFLAGS) $(CXXOPTS) -c -o obj/pma.o source/pma.cpp
	
clean:
	rm obj/*.o
