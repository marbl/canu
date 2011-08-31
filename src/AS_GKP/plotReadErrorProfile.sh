#! /bin/sh

export READ_ERROR_PROFILE="$1"

gnuplot -persist "${0%.sh}.gp"