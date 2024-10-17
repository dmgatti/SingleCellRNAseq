#!/usr/bin/env bash

# Only try running R to translate files if there are some files present.
# The Makefile passes in the names of files.

if [ $# -eq 2 ] ; then
    if [ "$1" -nt "$2" ]; then
	echo "Knitting $2 --> $1"
	Rscript -e "source('bin/generate_md_episodes.R')" "$@"
    else
	echo "$2 is up to date, skipping knitting"
    fi
fi
