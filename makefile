SOURCES = CMixtureModel.cpp CModel.cpp CSimpleGaussianModel.cpp CTransitionModel_SimpleGaussian.cpp CUniformModel.cpp CBoundedModel.cpp AddScaledLogs.cpp
OBJS = CMixtureModel.o CModel.o CSimpleGaussianModel.o CTransitionModel_SimpleGaussian.o CUniformModel.o CBoundedModel.o AddScaledLogs.o
#EXECUTABLE = test_gaussian_mixture

CPP = gcc
CPPFLAGS := $(CPPFLAGS) -g -Wall  
INCLUDE_DIR := $(INCLUDE_DIR) -I/usr/include/gsl -I/home/f1hxw01/equal_energy_hw/include
 
all : $(OBJS) 

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIR) -c $< -o $@ 

clean:
	rm -rf *.o $(EXECUTABLE)
