#!/bin/bash

if [[ -z $1 || -z $2 ]]; then
    echo "Usage: ./makeCardsWorkspaces.sh <massfile> <sqrts>"
    exit;
fi

# if sqrts == 0, do combined 7+8 TeV
# if sqrts == -2, do combined 7+8 TeV with tagged categories

# make combined cards
for M in $(cat $1); do bash make_combined_cards_hzz4l.sh $M $2 ; done

# make combined workspaces
for M in $(cat $1); do bash make_binary_workspaces_hzz4l.sh $M comb_hzz4l.txt; done

