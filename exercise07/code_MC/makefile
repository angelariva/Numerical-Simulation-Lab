CXXFLAGS := $(CXXFLAGS) -std=c++11  -O3 -Wall
#it allows to put object files in a subfolder
obj=$(patsubst %.o,$(OBJPATH)%.o, random.o main.o montecarlo.o)

		#source path
VPATH=./source/
		#object files path
OBJPATH=./obj/
		#executable name
EXENAME=main

all: $(obj)
ifeq ($(wildcard frames*),) #search for frames path, create it if there isn't
	@mkdir -p frames
endif
ifeq ($(wildcard results*),) #search for results path, create it if there isn't
	@mkdir -p results
endif
	$(CXX) $(OBJPATH)*.o -o $(EXENAME) $(CXXFLAGS)

$(OBJPATH)%.o: %.cpp
ifeq ($(wildcard $(OBJPATH)*),) #search for obj path, create it if there isn't
	@mkdir -p $(OBJPATH)
endif
	$(CXX) -c $< -o $@ $(CXXFLAGS)

.PHONY: clean
clean:
	rm $(OBJPATH)*.o $(EXENAME)
