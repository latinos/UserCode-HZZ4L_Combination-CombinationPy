Quick instructions:

All scripts produce configuration file fragments under the directory
CardFragments.

All common configurable parameters are specified in the file Config.h.

DEPENDENCIES:

cd ZZ4L_Combination/CombinationPy/CreateDatacards
cmsrel CMSSW_6_1_1
cd CMSSW_6_1_1/src
cvs co -r V03-01-02 HiggsAnalysis/CombinedLimit 
cvs co -r V00-03-01 -d Higgs/Higgs_CS_and_Width UserCode/Snowball/Higgs/Higgs_CS_and_Width 

To run californiaSignalShapes.C you need to download the shapes. In CreateDataCards/CMSSW_6_1_1/src
cvs co -r 1.38 -d WWAnalysis/TreeModifiers/macro UserCode/Mangano/WWAnalysis/TreeModifiers/macro/SignalInterpolationStrings.h


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
Yields number come from the wiki.
Use the script:

root -q -b zjetRates.cc

that creates the fragments:
CardFragments/zjetRate_[sqrts]TeV_[fs].txt

Shapes parameters taken from the wiki should be updated by hand here:
CardFragments/zjetShape_[sqrts]TeV_[fs].txt


4. 1D SIGNAL SHAPES 
Run:

root -q -b californiaSignalShapes.C 

Parameters to be used in the config files are in:
CardFragments/signalFunctions_[sqrts]TeV_[fs].txt

For mass/width measurement, the shapes are slightly different; use:

root -q -b californiaSignalShapesMassWidth.C

which creates the files
CardFragments/signalFunctionsMW_[sqrts]TeV_[fs].txt


4.1  signal EBE  pdf 
Run:

root -l -b -q californiaSignalEBE.C

5. 1D BACKGROUND SHAPES
Run:

root -q -b backgroundFits_qqzz_1Dw.C
root -q -b backgroundFits_ggzz_1Dw.C

Output parameters are written in:
CardFragments//[xx]zzBackgroundFit_[sqrts]TeV_[fs].txt
Shape plots are saved in the subdirectories bkgFigs7TeV, bkgFigs8TeV.


5.1. BACKGROUND EBE SHAPES
(come from external code)
Parameters to be used in the config files are in:
CardFragments/bkgEBE_[sqrts]TeV_[fs].txt


6. MERGE CARD FRAGMENTS
Run:

root -q -b mergeFragments.C+

The full config files are written under 
../CreateDatacards/SM_inputs_?TeV_tagged/ 

For mass/width measurement, edit mergeFragments.C and set
bool forMass=true;
before running.


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

9. TAGGED 2D TEMPLATES
Run:

root -q -b generateFisherTemplates.C+

The resulting templates are written into the final destination directory, i.e.
../CreateDatacards/templates2D
(for ggH, VBF, qqZZ four systematic templates are created, two for JEC, two for shape systematics. Use shape systematics (_alt, _alt2) for analysis)
