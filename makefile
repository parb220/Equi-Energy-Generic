SOURCES = CEquiEnergy.cpp CMixtureModel.cpp CModel.cpp CSimpleGaussianModel.cpp CTransitionModel_SimpleGaussian.cpp CUniformModel.cpp CBoundedModel.cpp test_gaussian_mixture.cpp
OBJS = $(SOURCES: .cpp = .o)
CPP = g++
DEBUG = -g
CPPFLAGS = -c -Wall $(DEBUG)
LINKFLAGS = -Wall $(DEBUG)
GSL_LINKFLAGS = -lgsl -lgslcblas -lm 
EXECUTABLE = test_gaussian_mixture
 
all : $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJS) 
	$(CPP) $(LINKFLAGS) $(OBJS) $(GSL_LINKFLAGS) -o $@

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $< -o $@ 

clean:
	rm -rf *.o $(EXECUTABLE)
