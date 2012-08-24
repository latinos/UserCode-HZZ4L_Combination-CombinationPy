#!/bin/bash
if [[ "$1" == "" ]]; then echo "Usage: $0 mass"; exit 1; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass"; exit 1; fi; 
cd $MASS

perl -i -npe 's/^#?(CMS_hww[os]f_stat_[01]j_(VH|ZH|WW|Wjets|ggH|ggWW)_bin1|QCDscale_(VH|ggH)_ACEPT)/#$1/' *hww[os]f_[01]j_shape.txt
perl -i -npe 's/^#?(CMS_hwwll_stat_2j_(qqH|VH|ZH|WW|ggH|ggWW))/#$1/' *hww_2j_cut.txt
