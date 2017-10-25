cmake_minimum_required(VERSION 2.4)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

PROJECT(mtrap)

# =========================================================================
# 1. How to install MTRAP without being root
#   Edit the OPTIONS as follows.
#   From:
#     DATA_DIR = /usr/local/share/mtrap
#   To:
#     DATA_DIR = /home/[your_home]/[somewhere]
# =========================================================================
# --- OPTIONS ---
set(mtrap_DATA_DIR "/usr/local/share/mtrap")
set(mtrap_VERSION_MAJOR 2)
set(mtrap_VERSION_MINOR 0)
set(mtrap_VERSION_PATCH 0)
set(mtrap_VERSION ${mtrap_VERSION_MAJOR}.${mtrap_VERSION_MINOR}.${mtrap_VERSION_PATCH})
# --- BUILD TYPE ---
set(CMAKE_BUILD_TYPE RELEASE)
#set(CMAKE_BUILD_TYPE DEBUG)
# 32bit compile with 64bit environment: Note that -m32 option needs g++-multilib package in ubuntu 64bit
#set (CMAKE_CXX_FLAGS "-mfpmath=sse -m32")

# =========================================================================
# build type depend flags
set(CMAKE_CXX_FLAGS_DEBUG "-msse2 -g -W -Wall -Wno-char-subscripts -Wno-sign-compare -Wno-unused-parameter -Wno-unused-function -Wno-self-assign")
set(CMAKE_CXX_FLAGS_RELEASE "-msse2 -Wno-deprecated -O3 -funroll-loops -ftree-vectorize")

# OPENMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# if defined then the probalign style alignment is used for the postprob routine, if not defined then the transiton quantity style is used.
#ADD_DEFINITIONS(-DUSE_PROBALIGN)

# TR1 LIBRARY SUPPORT
include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX (tr1/memory HAVE_SYS_TR1_MEMORY_H)
if (HAVE_SYS_TR1_MEMORY_H)
    ADD_DEFINITIONS(-DHAVE_SYS_TR1_MEMORY_H)
endif()

# isinf SUPPORT
include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS (isinf HAVE_DECL_ISINF)
if (HAVE_DECL_ISINF)
	ADD_DEFINITIONS(-DHAVE_DECL_ISINF)
endif()

# use snappy library
ADD_DEFINITIONS(-DCOMPRESS)

# use eigen library
ADD_DEFINITIONS(-DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET)

# add debug flag when debug mode
IF(CMAKE_BUILD_TYPE STREQUAL DEBUG)
        MESSAGE("Adding Debug flag...")
        ADD_DEFINITIONS(-DDEBUG)
ENDIF(CMAKE_BUILD_TYPE STREQUAL DEBUG)

# OS check
MESSAGE("OS is ${CMAKE_SYSTEM_NAME}.")
IF(${UNIX})
        MESSAGE("Adding UNIX flag...")
        ADD_DEFINITIONS(-DUNIX)
ENDIF(${UNIX})

# add include path
include_directories("${PROJECT_SOURCE_DIR}/lib/snappy")
include_directories("${PROJECT_SOURCE_DIR}/lib/eigen")
include_directories("${PROJECT_SOURCE_DIR}/lib")

# add library path
#link_directories(dir1 dir2 ...)

# option file will be generated in binary dir.
configure_file(
    "${PROJECT_SOURCE_DIR}/src/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
)
include_directories("${PROJECT_BINARY_DIR}")

# call subdirectories
SUBDIRS(src)

# === INSTALL INFORMATION ===
# cygwin support
IF(${CYGWIN})
	install(PROGRAMS "${PROJECT_BINARY_DIR}/src/mtrap.exe" DESTINATION bin)
ELSE(${CYGWIN})
	install(PROGRAMS "${PROJECT_BINARY_DIR}/src/mtrap" DESTINATION bin)
ENDIF(${CYGWIN})

# add resources
#FILE(GLOB DATA_FILES "${PROJECT_SOURCE_DIR}/resource/*[^~]")
#install(FILES ${DATA_FILES} DESTINATION ${mtrap_DATA_DIR})
install(FILES ${PROJECT_SOURCE_DIR}/resource/SABmark1.63_sup_weighted.btq DESTINATION ${mtrap_DATA_DIR})
install(FILES ${PROJECT_SOURCE_DIR}/resource/VTML200I DESTINATION ${mtrap_DATA_DIR})
install(FILES ${PROJECT_SOURCE_DIR}/resource/scop1.75-id-with-sab1.65-group.csv DESTINATION ${mtrap_DATA_DIR})
FILE(GLOB FAMILYDATABASE_FILES "${PROJECT_SOURCE_DIR}/resource/fsmatrix-sabmark1.65/*[^~]")
install(FILES ${FAMILYDATABASE_FILES} DESTINATION ${mtrap_DATA_DIR}/fsmatrix-sabmark1.65)

# add executables
#install(TARGETS mtrap executable DESTINATION "/usr/local/bin")

# build a CPack driven installer package
include(InstallRequiredSystemLibraries)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MTRAP")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION_MAJOR "${mtrap_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${mtrap_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${mtrap_VERSION_PATCH}")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Toshihide Hara")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "mtrap")
set(CPACK_SOURCE_IGNORE_FILES "/\\\\..*$;~$;${CPACK_SOURCE_IGNORE_FILES}")
#SET(CPACK_SYSTEM_NAME "i386")
SET(CPACK_GENERATOR "DEB")
#set(CPACK_SET_DESTDIR "ON")

# TODO: support for MAC OS-X

include (CPack)

