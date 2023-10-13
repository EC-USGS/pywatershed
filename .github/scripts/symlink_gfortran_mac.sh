#!/usr/bin/env bash

# get full gfortran version string
# assumes installed via brew as by
# https://github.com/fortran-lang/setup-fortran
#
# sed not head for first line, avoid ruby broken pipe issues
# (https://stackoverflow.com/a/2845541/6514033)
full_version=$(brew info gfortran | sed -n 1p | cut -d' ' -f 4)

# get major version
version=$(echo "$full_version" | cut -d'.' -f 1)

# symlink gfortran libraries
old_libdir="/usr/local/opt/gcc/lib/gcc/${version}"
new_libdir="/usr/local/lib/"
mkdir -p "$new_libdir"
if [ -d "$old_libdir" ]
then
  sudo ln -fs "$old_libdir/libgfortran.5.dylib" "$new_libdir/libgfortran.5.dylib"
  sudo ln -fs "$old_libdir/libquadmath.0.dylib" "$new_libdir/libquadmath.0.dylib"
fi