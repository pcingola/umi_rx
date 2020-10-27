#!/bin/bash -eu
set -o pipefail

cd src

mkdir -p ../bin
gcc -o ../bin/umi_rx -L../htslib/lib -l hts umi_rx.c
