
name := he3
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = /home/busenitz/lz/AmLi/simulation/detector/he-3
endif

CPPFLAGS += -DKL_USE_ROOT=1 $(shell root-config --cflags)
EXTRALIBS += $(shell root-config --libs)


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk


visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

