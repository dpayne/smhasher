#!/bin/sh
rm -rf build
mkdir build
cd build
# Note that gcc-mp-5 on darwin is still ~2000 MiB/sec faster than the system clang
# clang-mp-dev (4.0) also
case `uname -s` in
#Darwin)  CXX=clang++-mp-devel  CC=clang-mp-devel cmake .. ;;
#Darwin)  CXX=g++-mp-6  CC=gcc-mp-6 cmake .. ;;
Darwin)  cmake .. ;;
FreeBSD) CXX=clang++ CC=clang cmake .. ;;
*)       CXX=clang++-4.0  CC=clang-4.0 cmake ..   ;;
esac
make VERBOSE=1 V=1 -j4 $@
cd ..
