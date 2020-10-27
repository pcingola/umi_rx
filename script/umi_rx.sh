#!/bin/bash

DIR=$( cd $(dirname "$0") ; pwd -P )
HTSLIB_PATH=$( cd "$DIR/../htslib/lib"; pwd -P )
UMI_RX_PATH=$( cd "$DIR/../bin"; pwd -P )

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HTSLIB_PATH"

"$UMI_RX_PATH/umi_rx" $*
