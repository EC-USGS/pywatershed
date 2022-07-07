#!/bin/bash

compiler=$1
debug=$2
clean=$3

# want to specify prms version?
cd prms_src/prms5.2.1/prms/

if [[ "$clean" == 'rebuild' ]]; then
    echo "rebuilding PRMS"
    make clean
    make
fi

if [ "$clean" == 'build' ]; then
    echo "Building PRMS"
    make
fi


exit
