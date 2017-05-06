#include make.inc
#################################################
### Defining Compiling variables              ###
#################################################

CXX           = g++ -Wno-deprecated -Wall
#-m32  -Wno-deprecated
LD            = g++ -Wno-deprecated -Wall
#-m32  -Wno-deprecated 
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

#SHAREDLIBFLAGS =-L./CommonUtils/lib -lTauNtuple -L./DataFormats/lib -lDataFormats
#SHAREDCXXFLAGS =-I./CommonUtils -I./DataFormats
##SUBDIRS = CommonUtils/ 
#TAUSPINERFLAG = 

override CXXFLAGS += 
ALLCXXFLAGS = $(CXXFLAGS) -I$(ROOTSYS)/include -c -I/include $(SHAREDCXXFLAGS) -c -I./ -c  -L -lshared
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lCint -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint -lReflexDict -lReflex -lPostscript -lMatrix -lPhysics -lGui -LObj -lMinuit2 -lThread -lTreePlayer -lRIO -lMathCore -rdynamic -Wl,--rpath $(ROOTSYS)/lib

# -ldl -lpthread

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS)# $(SHAREDLIBFLAGS)

TARGETS       =	MultiplyNumbers \
		                AddNumbers


CINTTARGETS = 


HDRS          = $(TARGETS:=.h) 		
SRCS          = $(TARGETS:=.cc)  main.cc
OBJS          = $(SRCS:.cc=.o)

PROGRAM       = run.exe

$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) *.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) *.o $(LIBS) -o $(PROGRAM)
	@echo "done"
 
VPATH = utilities:i386_linux


#	@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o ../$(PROGRAM)

##########################################################
###                                                    ###
##########################################################
#	$(CC) -Wall -L -lshared -o run main.cc  
$(OBJS): %.o : %.cc
		$(CXX) $(ALLCXXFLAGS) $(DEFS) $< -o $@ 

.PHONY: clean cleanall cleandf all dataformats install sharedlib 

install:  Analysis.exe

clean:
	@rm *.o
	@rm run.exe

cleandf:
	@cd DataFormats; gmake clean; cd ../

cleanall:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone clean && cd ../../; \
	done
	@cd DataFormats; gmake clean; cd ../
	@rm *.o
	@rm run.exe

all: sharedlib dataformats install























#http://www.iram.fr/~roche/code/c++/AddNumbers.html
#Create a static library
#First the source file src/AddNumbers.cpp is turned into an object file.
#[~/workspace/C++/AddNumbers] > g++ -I ./inc -c src/AddNumbers.cpp -o obj/AddNumbers.o#
#A static library is basically a set of object files that were copied into a single file. It is created invoking the archiver ar. The library name must start with the three letters lib and have the suffix .a.
#[~/workspace/C++/AddNumbers] > ar rcs lib/libAddNumbers.a obj/AddNumbers.o
#You can also write similar rules in a makefile. See the file Makefile.static given in the AddNumbers.tar.bz2 archive.#
#Refer to useful options of g++ for details.
#reate a shared library
#The -fpic option tells g++ to create position independant code which is needed for shared libraries.
#[~/workspace/C++/AddNumbers] > g++ -I ./inc -fpic -c src/AddNumbers.cpp -o obj/AddNumbers.o
#Finally the shared library is created. Note the library name must start with the three letters lib and have the suffix .so.
#[~/workspace/C++/AddNumbers] > g++ -shared -o lib/libAddNumbers.so obj/AddNumbers.o
#As a makefile example see the file Makefile.shared given in the AddNumbers.tar.bz2 archive.
#Refer to useful options of g++ for details.