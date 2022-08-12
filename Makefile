TARGET=GPURADMC

CXX=g++
CXXFLAGS=-std=c++14 -g -ggdb -Wall -O3
INCLUDE=-I./include
OBJS=obj/dust_specie.o obj/dust.o obj/star.o obj/stars.o obj/frequencies.o obj/grid.o obj/radmc.o

all:
	 make $(TARGET)

$(TARGET):$(OBJS)
		    $(CXX) $^ -o $@ $(CXXFLAGS) $(INCLUDE)

obj/dust_specie.o:src/dust_species.cc include/dust_species.hh
				$(CXX) $< -c -o $@ $(CXXFLAGS) $(INCLUDE)

obj/dust.o:src/dust.cc include/dust.hh
				$(CXX) $< -c -o $@ $(CXXFLAGS) $(INCLUDE)

obj/star.o:src/star.cc include/star.hh
				$(CXX) $< -c -o $@ $(CXXFLAGS) $(INCLUDE)

obj/stars.o:src/stars.cc include/stars.hh
				$(CXX) $< -c -o $@ $(CXXFLAGS) $(INCLUDE)

obj/frequencies.o:src/frequencies.cc include/frequencies.hh
				$(CXX) $< -c -o $@ $(CXXFLAGS) $(INCLUDE)

obj/grid.o:src/grid.cc include/grid.hh
				$(CXX) $< -c -o $@ $(CXXFLAGS) $(INCLUDE)

obj/radmc.o:src/radmc.cc
				$(CXX) $< -c -o $@ $(CXXFLAGS) $(INCLUDE)

clean:
	${RM} $(OBJS) $(TARGET)