Quick instructions:

All scripts produce configuration file fragments under the directory
CardFragments.
The merger of these fragments 

All common configurable parameters are specified in the file Config.h.

1. SIGNAL EFFICIENCIES
Run:

root -q -b signalEfficiency_w.C+

Output parameters are written in:
CardFragments/signalEfficiency_[sqrts]TeV_[fs].txt
Efficiency plots are saved in the subdirectories sigFigs7TeV, sigFigs8TeV.

2. ZZ BACKGROUND RATES
Run:

root -q -b ZZbackgroundRate.C+

Output parameters are written in:
CardFragments/ZZRates_[sqrts]TeV_[fs].txt

3. ZJET RATES AND SHAPES
(come from external code)
Parameters to be used in the config files are in:
CardFragments/zjetRate_[sqrts]TeV_[fs].txt
CardFragments/zjetShape_[sqrts]TeV_[fs].txt


4. 1D SIGNAL SHAPES 
Run:

root -q signalFits.C

Parameters to be used in the config files are in:
CardFragments/signalFunctions_[sqrts]TeV_[fs].txt


5. 1D BACKGROUND SHAPES
Run:

root -q -b backgroundFits_qqzz_1Dw.C
root -q -b backgroundFits_ggzz_1Dw.C

Output parameters are written in:
CardFragments//[xx]zzBackgroundFit_[sqrts]TeV_[fs].txt
Shape plots are saved in the subdirectories bkgFigs7TeV, bkgFigs8TeV.


6. MERGE CARD FRAGMENTS
Run:

root -q -b mergeFragments.C+

The full config files are written under 
../CreateDatacards/SM_inputs_?TeV_tagged/ 



7. DATA FILES
Run:

root -q -b prepareData.C+

This creates all files for 3 final states, 7 and 8 TeV, and for the tagged/untagged and stores them
in the final destination directory, i.e. 
../CreateDatacards/CMSdata/


8. SIGNAL AND BACKGROUND 2D TEMPLATES
Run:

root -q -b loadMELA.C generateTemplates.C+

The resulting templates are written into the final destination directory, i.e. 
../CreateDatacards/templates2D
