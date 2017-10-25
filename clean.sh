#!/bin/sh
find . -name "*~" -exec rm {} \;
rm Makefile 2> /dev/null
rm install_manifest.txt 2> /dev/null
rm CMakeCache* 2> /dev/null
rm CMakeLists.txt.user* 2> /dev/null
rm CPack* 2> /dev/null
rm cmake_install.cmake 2> /dev/null
rm config.h 2> /dev/null
rm CMakeFiles 2> /dev/null
rm -R _CPack_Packages 2> /dev/null
rm -R ../mtrap-build 2> /dev/null
rm -R ./src/CMakeFiles 2> /dev/null
rm ./src/config.h 2> /dev/null
rm ./src/cmake_install.cmake 2> /dev/null
rm ./src/CMakeLists.txt.user 2> /dev/null
rm ./src/test_tree 2> /dev/null
rm ./src/mtrap 2> /dev/null

