CXX	:= g++
CXXFLAGS := -O2 -fpic
    
MODULES = ini.o INIReader.o

libinih.so: $(MODULES)
	$(CXX) -shared ini.o INIReader.o -o libinih.so

%.o: %.cpp %.c
	$(CXX) $(CXXFLAGS) -c $(input) -o $(output)
