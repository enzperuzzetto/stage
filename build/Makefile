CXX=g++
CXXFLAGS=-Wall -g -O3
LDFLAGS=-Wall
EXEC=AnalogImage

all: $(EXEC)

AnalogImage: loadSaveTiffClass.o main.o classeFeatures.o Point.o annUse.o BestCoherence.o -ltiff -lANN
	$(CXX) -o $@ $^ $(LDFLAGS)
	
%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(EXEC) loadSaveTiffClass.o main.o classeFeatures.o Point.o annUse.o BestCoherence.o