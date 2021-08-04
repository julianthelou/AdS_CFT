/* This file is autogenerated by generate-buildinfo.sh */

#ifndef _EXAHYPE_BUILD_INFO_H_
#define _EXAHYPE_BUILD_INFO_H_

#define EXAHYPE_BUILDINFO_AVAILABLE
#define EXAHYPE_BUILD_DATE          "Tue 27 Jul 2021 04:27:00 PM CEST"
#define EXAHYPE_BUILD_HOST          "jujuu-ThinkPad-X220"

/* Strings passed by the Makefile */
#define EXAHYPE_BUILD_INFO \
    "COMPILER = GNU\n" \
    "MODE = Asserts\n" \
    "SHAREDMEM = None\n" \
    "DISTRIBUTEDMEM = None\n" \
    "------\n" \
    "ARCHITECTURE = snb\n" \
    "CC = g++\n" \
    "BOUNDARYCONDITIONS = None\n" \
    "------\n" \
    "COMPILER_CFLAGS = -DAsserts -DTrackGridStatistics -std=c++11 -pedantic -Wall -Drestrict=__restrict__ -pipe -D__assume_aligned=__builtin_assume_aligned -Wstrict-aliasing -fopenmp-simd -g3 -march=sandybridge\n" \
    "COMPILER_LFLAGS = -lm -lstdc++ -g3 -lrt -march=sandybridge\n" \
    "FCOMPILER_CFLAGS = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp -Wall -Wno-tabs -Wno-unused-variable -fopenmp-simd -g3\n" \
    "FCOMPILER_LFLAGS = \n" \
    "------\n" \
    "EXAHYPE_PATH = /home/jujuu/exahype/ExaHyPE-Engine/./ExaHyPE\n" \
    "PROJECT_PATH = /home/jujuu/exahype/ExaHyPE-Engine/./AstroApplications/CCZ4/FV\n" \
    "PEANO_KERNEL_PEANO_PATH = /home/jujuu/exahype/ExaHyPE-Engine/./Peano/peano\n" \
    "PEANO_KERNEL_TARCH_PATH = /home/jujuu/exahype/ExaHyPE-Engine/./Peano/tarch\n" \
    "------\n" \
    "PROJECT_CFLAGS = -DDim3 -DALIGNMENT=32\n" \
    "PROJECT_LFLAGS = \n" \
    ""


/*
 * ExaHyPE Git Repository information extraction
 * Extracted from git repository in /home/jujuu/exahype/ExaHyPE-Engine/./ExaHyPE
 */

/* Information collected with git version 2.25.1 */
#define EXAHYPE_GIT_INFO "master  f12ff8b56 Sun Jul 4 04:52:32 2021"

/*
 * Peano Git Repository information extraction
 * Extracted from git repository in /home/jujuu/exahype/ExaHyPE-Engine/./Peano/peano
 */

/* Information collected with git version 2.25.1 */
#define PEANO_GIT_INFO "HEAD  40b5c080 Tue Oct 13 22:36:28 2020"
#define PEANO_SVN_INFO  PEANO_GIT_INFO /* transition time */

/*
    Peano version check

    FIXME TODO This is the worst place ever to hook in version
               requirements. Please move the following to somewhere
               suitable, such as the Peano startup phase.

*/

#include "peano/version.h"

#endif /* EXAHYPE_BUILD_INFO_H */
#if PEANO_VERSION<2509 
#error Old Peano version. Version 2509 required. Please update your Peano installation.
#endif
