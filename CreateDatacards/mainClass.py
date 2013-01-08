#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from datacardClass import *
from kdClass import *

class mainClass():
    
    ID_4mu = 1
    ID_4e  = 2
    ID_2e2mu = 3
    

    def __init__(self):

        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
        ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")
        

    def makeCardsWorkspaces(self, theMH, theOutputDir, theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD):


        ## ------------------ CHECK CHANNEL ------------------- ##
        print theInputs['decayChannel']

        if (theInputs['decayChannel'] == self.ID_4mu): appendName = '4mu'
        elif (theInputs['decayChannel'] == self.ID_4e): appendName = '4e'
        elif (theInputs['decayChannel'] == self.ID_2e2mu): appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"
        
        #if (theInputs['decayChannel'] != self.ID_4mu and theInputs['decayChannel'] != self.ID_4mu and theInputs['decayChannel'] != self.ID_4mu):
        #    raise RuntimeError, "Error: channel is not an option (4mu,4e,2e2mu)" 

        
        ## ----------------- WIDTH AND RANGES ----------------- ##
        myCSW = HiggsCSandWidth()
        widthHVal =  myCSW.HiggsWidth(0,theMH)

        self.windowVal = max( widthHVal, 1.0)
        self.windowVal = max( widthHVal, 1.0)
        lowside = 100.0
        highside = 1000.0
        if (theMH >= 275):
            lowside = 180.0
            highside = 650.0
        if (theMH >= 350):
            lowside = 200.0
            highside = 900.0
        if (theMH >= 500):
            lowside = 250.0
            highside = 1000.0
        if (theMH >= 700):
            lowside = 350.0
            highside = 1400.0
                        
        self.low_M = max( (theMH - 20.*self.windowVal), lowside)
        self.high_M = min( (theMH + 15.*self.windowVal), highside)
               
        print "Higgs Width: ",widthHVal, " Window --- Low: ", self.low_M, " High: ", self.high_M

        ## add to the inputs
        theInputs['low_M'] = self.low_M
        theInputs['high_M'] = self.high_M
        
        

        ## ------------- MELA for Signal Separation ----------- ##
        self.isAltSig = False
        if (theInputs['doHypTest']):
            self.isAltSig = True
            
        if self.isAltSig and not self.all_chan :
            raise RuntimeError, "You asked to prepare DC and WS for Hyp Test but you did not want to sum over all signal channels. This is forbidden. Check inputs ! (it should have already send you this error message, strange that  you are here...)"
        
        if (self.isAltSig and not (self.is2D==1)):
            raise RuntimeError, "Cannot perform hypothesis testing without a 2D analysis, feature not supported yet. Exiting."
        
        
        self.appendHypType = theInputs['altHypLabel']
        if self.isAltSig and self.appendHypType=="" :
            self.appendHypType = "_ALT"


        myDatacardClass = datacardClass()

        if( theis2D == 0 ):

            myDatacardClass.makeMassShapesYields(theMH,theOutputDir,theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD)
            myDatacardClass.fetchDataset()
            myDatacardClass.writeWorkspace()
            myDatacardClass.prepareDatacard()

        if( theis2D == 1 ):

            myDatacardClass2D = kdClass()
            myDatacardClass2D.makeMassShapesYields(theMH,theOutputDir,theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD)

            myDatacardClass2D.setKD()
            myDatacardClass2D.fetchDatasetKD()            
            myDatacardClass2D.makeKDAnalysis()
            myDatacardClass2D.writeWorkspaceKD()
            myDatacardClass2D.prepareDatacardKD()
            
            

            
        

