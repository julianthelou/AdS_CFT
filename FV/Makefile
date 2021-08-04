# *********************************************************************************************
# README
# *********************************************************************************************
#
# Available configuration Parameters for ExaHyPE
#
# export variable  |  default-value  |  further values         |  description
#-------------------------------------------------------------------------------------------------------
# COMPILER            Intel             GNU                          Used compiler (and linker)
# MODE                Release           Debug,                       Verbosity and Debug level
#                                       Asserts,PeanoProfile,Profile          
# SHAREDMEM           None              OMP, TBB                     Shared-memory parallelisation
# DISTRIBUTEDMEM      None              MPI                          Distributed-memory parallelisation
# BOUNDARYCONDITIONS  None              Periodic                     Type of boundary conditions
# USE_IPO             Off               On                           IPO during compilation (intel only)
# ******************************************************************************************************

# Helper
SHELL = bash
tolower = $(shell echo $(1) | tr '[:upper:]' '[:lower:]')

# Default Paths and name
PEANO_KERNEL_PEANO_PATH=/home/jujuu/exahype/ExaHyPE-Engine/./Peano/peano
PEANO_KERNEL_TARCH_PATH=/home/jujuu/exahype/ExaHyPE-Engine/./Peano/tarch
PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=/home/jujuu/exahype/ExaHyPE-Engine/./Peano/multiscalelinkedcell
PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=/home/jujuu/exahype/ExaHyPE-Engine/./Peano/sharedmemoryoracles
PEANO_TOOLBOX_MPI_BLANCING_PATH=/home/jujuu/exahype/ExaHyPE-Engine/./Peano/mpibalancing
EXAHYPE_PATH=/home/jujuu/exahype/ExaHyPE-Engine/./ExaHyPE
PROJECT_PATH=/home/jujuu/exahype/ExaHyPE-Engine/./AstroApplications/CCZ4/FV
EXECUTABLE=ExaHyPE-CCZ4

# Default flags from specification file constants
PROJECT_CFLAGS+= -DDim3
ARCHITECTURE=snb
PROJECT_CFLAGS+= -DALIGNMENT=32

# Default sharedmem/distributedmem required by the specification file
ifeq ($(SHAREDMEM),)
    SHAREDMEM=None
endif
ifeq ($(DISTRIBUTEDMEM),)
    DISTRIBUTEDMEM=None
endif

# If using TBB only, link pthread (comment this if you want to link a specific pthread)
ifeq ($(call tolower,$(SHAREDMEM)),tbb)
    ifeq ($(call tolower,$(DISTRIBUTEDMEM)),none)
        PROJECT_LFLAGS+= -pthread
    endif
endif

# Default options required by the specification file



# Main Makefile include
-include /home/jujuu/exahype/ExaHyPE-Engine/./ExaHyPE/Makefile


all: 
	@echo CCZ4
	@echo =================
	@echo An ExaHyPE solver