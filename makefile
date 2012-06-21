SOURCES = CEquiEnergy.cpp CMixtureModel.cpp CModel.cpp CSimpleGaussianModel.cpp CTransitionModel_SimpleGaussian.cpp CUniformModel.cpp CBoundedModel.cpp test_gaussian_mixture.cpp
OBJS = $(SOURCES: .cpp = .o)
CPP = g++
DEBUG = -g
CPPFLAGS = -c -Wall $(DEBUG)
LINKFLAGS = -Wall $(DEBUG)
GSL_LINKFLAGS = -lgsl -lgslcblas -lm 
EXECUTABLE = test_gaussian_mixture
INCLUDE_DIR := $(INCLUDE_DIR) -I/home/f1hxw01/equal_energy_hw/include
 
all : $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJS) 
	$(CPP) $(LINKFLAGS) $(OBJS) $(GSL_LINKFLAGS) -o $@

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@ 

clean:
	rm -rf *.o $(EXECUTABLE)
