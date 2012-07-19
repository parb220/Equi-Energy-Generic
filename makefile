SOURCES = CEquiEnergy.cpp CMixtureModel.cpp CModel.cpp CSimpleGaussianModel.cpp CTransitionModel_SimpleGaussian.cpp CUniformModel.cpp CBoundedModel.cpp test_gaussian_mixture.cpp
OBJS = $(SOURCES: .cpp = .o)
EXECUTABLE = test_gaussian_mixture

CPP = gcc
CPPFLAGS := $(CPPFLAGS) -g -Wall 
LINKFLAGS := $(LINKFLAGS) -lstdc++ -lgsl -lgslcblas -lm 
INCLUDE_DIR := $(INCLUDE_DIR) -I/usr/include/gsl -I/home/f1hxw01/equal_energy_hw/include
LINK_DIR := $(LINK_DIR) -L/usr/lib64
 
all : $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJS) 
	$(CPP) $(CPPFLAGS) $(LINK_DIR) $(LINKFLAGS) $(OBJS) -o $@

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@ 

clean:
	rm -rf *.o $(EXECUTABLE)
