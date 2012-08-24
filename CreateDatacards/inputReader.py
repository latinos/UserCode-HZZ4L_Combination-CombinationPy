#!/usr/bin/python
import os
import re
import math
import collections
from ROOT import *
from array import array

## ---------------------------------------------------------------
## card reader class
## ---------------------------------------------------------------

class inputReader:

    def __init__(self, inputTextFile):

        if not os.path.exists(inputTextFile):
            raise RuntimeError, "File {0} does not exist!!!".format(inputTextFile)
        
        # input file
        self.theInput = inputTextFile
        # model
        self.model = ""
        # decay channel
        self.decayChan = -999.9
        # lumi
        self.lumi = -999.9
        # sqrts
        self.sqrts = -999.9
        # channels
        self.ggH_chan = False
        self.qqH_chan = False
        self.WH_chan = False
        self.ZH_chan = False
        self.ttH_chan = False
        self.qqZZ_chan = False
        self.ggZZ_chan = False
        self.zjets_chan = False
        self.ttbar_chan = False
        self.zbb_chan = False
        # rates
        self.qqZZ_rate = -999.9
        self.ggZZ_rate = -999.9
        self.zjets_rate = -999.9
        self.ttbar_rate = -999.9
        self.zbb_rate = -999.9
        self.qqZZ_lumi = -999.9
        self.ggZZ_lumi = -999.9
        self.zjets_lumi = -999.9
        self.ttbar_lumi = -999.9
        self.zbb_lumi = -999.9
        # signal shapes
        self.n_CB_shape = -999.9
        self.alpha_CB_shape = -999.9
        self.mean_CB_shape = -999.9
        self.sigma_CB_shape = -999.9
        # signal efficiency params
        self.sigeff_a1 = -999.9
        self.sigeff_a2 = -999.9
        self.sigeff_a3 = -999.9
        self.sigeff_a4 = -999.9
        self.sigeff_b1 = -999.9
        self.sigeff_b2 = -999.9
        self.sigeff_b3 = -999.9
        # qqZZ shape
        self.qqZZshape_a0 = -999.9
        self.qqZZshape_a1 = -999.9
        self.qqZZshape_a2 = -999.9
        self.qqZZshape_a3 = -999.9
        self.qqZZshape_a4 = -999.9
        self.qqZZshape_a5 = -999.9
        self.qqZZshape_a6 = -999.9
        self.qqZZshape_a7 = -999.9
        self.qqZZshape_a8 = -999.9
        self.qqZZshape_a9 = -999.9
        self.qqZZshape_a10 = -999.9
        self.qqZZshape_a11 = -999.9
        self.qqZZshape_a12 = -999.9
        self.qqZZshape_a13 = -999.9
        # ggZZ shape
        self.ggZZshape_a0 = -999.9
        self.ggZZshape_a1 = -999.9
        self.ggZZshape_a2 = -999.9
        self.ggZZshape_a3 = -999.9
        self.ggZZshape_a4 = -999.9
        self.ggZZshape_a5 = -999.9
        self.ggZZshape_a6 = -999.9
        self.ggZZshape_a7 = -999.9
        self.ggZZshape_a8 = -999.9
        self.ggZZshape_a9 = -999.9
        # zjets shape
        self.zjetsShape_mean = -999.9
        self.zjetsShape_sigma = -999.9
        # systematics 
        self.zjetsKappaLow = -999.9
        self.zjetsKappaHigh = -999.9


    def goodEntry(self,variable):
        if variable == -999.9:
            return False
        else:
            return True

    def readInputs(self):
        
        for line in open(self.theInput,'r'):
            f = line.split()
            if len(f) < 1: continue

            if f[0].startswith("#"): continue
            
            if f[0].lower().startswith("model"):
                
                if f[1].upper() == "SM": self.model = "SM"
                elif f[1].upper() == "SM4": self.model = "SM4"
                elif f[1].upper() == "FF" or f[1].upper() == "FP": self.model = "FF"
                else : raise RuntimeError, "Unknow model {0}, choices are SM, SM4, or FF".format(f[1].upper()) 

            if f[0].lower().startswith("decay"):

                if f[1] == "4mu": self.decayChan = 1
                elif f[1] == "4e": self.decayChan = 2
                elif f[1] == "2e2mu": self.decayChan = 3
                elif f[1] == "2mu2e": self.decayChan = 3
                else : raise RuntimeError, "Unknown decay channel {0}, choices are 4mu, 4e, or 2e2mu".format(f[1])
                
            if f[0].lower().startswith("channels"):
                for chan in f:
                    if chan == f[0]: continue
                    if chan.lower().startswith("ggh"):     self.ggH_chan = True
                    elif chan.lower().startswith("qqh"):   self.qqH_chan = True
                    elif chan.lower().startswith("wh"):    self.WH_chan = True
                    elif chan.lower().startswith("zh"):    self.ZH_chan = True
                    elif chan.lower().startswith("tth"):   self.ttH_chan = True
                    elif chan.lower().startswith("qqzz"):  self.qqZZ_chan = True
                    elif chan.lower().startswith("ggzz"):  self.ggZZ_chan = True
                    elif chan.lower().startswith("zjets"): self.zjets_chan = True
                    elif chan.lower().startswith("ttbar"): self.ttbar_chan = True
                    elif chan.lower().startswith("zbb"):   self.zbb_chan = True
                    else : raise RuntimeError, "Unknown channel {0}, choices are ggH, qqH, WH, ZH, ttH, qqZZ, ggZZ, zjets".format(chan)
                    
            if f[0].lower().startswith("rate"):
                
                if f[1].lower().startswith("qqzz"):
                    self.qqZZ_rate = float(f[2])
                    if len(f) == 4: self.qqZZ_lumi = float(f[3])
                if f[1].lower().startswith("ggzz"):
                    self.ggZZ_rate = float(f[2])
                    if len(f) == 4: self.ggZZ_lumi = float(f[3])
                if f[1].lower().startswith("zjets"):
                    self.zjets_rate = float(f[2])
                    if len(f) == 4: self.zjets_lumi = float(f[3])
                if f[1].lower().startswith("ttbar"):
                    self.ttbar_rate = float(f[2])
                    if len(f) == 4: self.ttbar_lumi = float(f[3])
                if f[1].lower().startswith("zbb"):
                    self.zbb_rate = float(f[2])
                    if len(f) == 4: self.zbb_lumi = float(f[3])
                    
            if f[0].lower().startswith("signalshape"):

                if f[1].lower().startswith("n_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.n_CB_shape = f[2]
                if f[1].lower().startswith("alpha_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.alpha_CB_shape = f[2]
                if f[1].lower().startswith("mean_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mean_CB_shape = f[2]
                if f[1].lower().startswith("sigma_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.sigma_CB_shape = f[2]
                    
            if f[0].lower().startswith("signaleff"):

                if f[1].lower().startswith("a1"): self.sigeff_a1 = float(f[2])
                if f[1].lower().startswith("a2"): self.sigeff_a2 = float(f[2])
                if f[1].lower().startswith("a3"): self.sigeff_a3 = float(f[2])
                if f[1].lower().startswith("a4"): self.sigeff_a4 = float(f[2])
                if f[1].lower().startswith("b1"): self.sigeff_b1 = float(f[2])
                if f[1].lower().startswith("b2"): self.sigeff_b2 = float(f[2])
                if f[1].lower().startswith("b3"): self.sigeff_b3 = float(f[2])

            if f[0].lower().startswith("qqzzshape"):

                if f[1].lower().startswith("a0"): self.qqZZshape_a0 = float(f[2])
                if f[1].lower().startswith("a1_") or f[1].lower().startswith("a1 "): self.qqZZshape_a1 = float(f[2])
                if f[1].lower().startswith("a2"): self.qqZZshape_a2 = float(f[2])
                if f[1].lower().startswith("a3"): self.qqZZshape_a3 = float(f[2])
                if f[1].lower().startswith("a4"): self.qqZZshape_a4 = float(f[2])
                if f[1].lower().startswith("a5"): self.qqZZshape_a5 = float(f[2])
                if f[1].lower().startswith("a6"): self.qqZZshape_a6 = float(f[2])
                if f[1].lower().startswith("a7"): self.qqZZshape_a7 = float(f[2])
                if f[1].lower().startswith("a8"): self.qqZZshape_a8 = float(f[2])
                if f[1].lower().startswith("a9"): self.qqZZshape_a9 = float(f[2])
                if f[1].lower().startswith("a10"): self.qqZZshape_a10 = float(f[2])
                if f[1].lower().startswith("a11"): self.qqZZshape_a11 = float(f[2])
                if f[1].lower().startswith("a12"): self.qqZZshape_a12 = float(f[2])
                if f[1].lower().startswith("a13"): self.qqZZshape_a13 = float(f[2])
                
            if f[0].lower().startswith("ggzzshape"):

                if f[1].lower().startswith("a0"): self.ggZZshape_a0 = float(f[2])
                if f[1].lower().startswith("a1"): self.ggZZshape_a1 = float(f[2])
                if f[1].lower().startswith("a2"): self.ggZZshape_a2 = float(f[2])
                if f[1].lower().startswith("a3"): self.ggZZshape_a3 = float(f[2])
                if f[1].lower().startswith("a4"): self.ggZZshape_a4 = float(f[2])
                if f[1].lower().startswith("a5"): self.ggZZshape_a5 = float(f[2])
                if f[1].lower().startswith("a6"): self.ggZZshape_a6 = float(f[2])
                if f[1].lower().startswith("a7"): self.ggZZshape_a7 = float(f[2])
                if f[1].lower().startswith("a8"): self.ggZZshape_a8 = float(f[2])
                if f[1].lower().startswith("a9"): self.ggZZshape_a9 = float(f[2])
               
            if f[0].lower().startswith("zjetsshape"):

                if f[1].lower().startswith("mean"):  self.zjetsShape_mean = f[2]
                if f[1].lower().startswith("sigma"): self.zjetsShape_sigma = f[2]
                

            if f[0].lower().startswith("systematic"):
                
                if f[1].lower().startswith("zjet") and f[1].lower().find("kappalow") >= 0 :
                    self.zjetsKappaLow = f[2]
                if f[1].lower().startswith("zjet") and f[1].lower().find("kappahigh") >= 0 :
                    self.zjetsKappaHigh = f[2]

            if f[0].lower().startswith("lumi"):

                self.lumi = float(f[1])

            if f[0].lower().startswith("sqrts"):

                self.sqrts = float(f[1])


    



    def getInputs(self):

        dict = {}

        ## check settings ##

        if not self.goodEntry(self.sqrts): raise RuntimeError, "{0} is not set.  Check inputs!".format("sqrts")

        if self.qqZZ_chan and not self.goodEntry(self.qqZZ_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZ_rate")
        if self.ggZZ_chan and not self.goodEntry(self.ggZZ_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZ_rate")
        if self.zjets_chan and not self.goodEntry(self.zjets_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("zjets_rate")
        if self.zbb_chan and not self.goodEntry(self.zbb_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("zbb_rate")
        if self.ttbar_chan and not self.goodEntry(self.ttbar_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("ttbar_rate")

        if not self.goodEntry(self.n_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("n_CB_shape")
        if not self.goodEntry(self.alpha_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("alpha_CB_shape")
        if not self.goodEntry(self.mean_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("mean_CB_shape")
        if not self.goodEntry(self.sigma_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigma_CB_shape")

        if not self.goodEntry(self.sigeff_a1): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_a1")
        if not self.goodEntry(self.sigeff_a2): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_a2")
        if not self.goodEntry(self.sigeff_a3): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_a3")
        if not self.goodEntry(self.sigeff_a4): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_a4")
        if not self.goodEntry(self.sigeff_b1): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_b1")
        if not self.goodEntry(self.sigeff_b2): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_b2")
        if not self.goodEntry(self.sigeff_b3): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_b3")

        if not self.goodEntry(self.qqZZshape_a0): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a0")
        if not self.goodEntry(self.qqZZshape_a1): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a1")
        if not self.goodEntry(self.qqZZshape_a2): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a2")
        if not self.goodEntry(self.qqZZshape_a3): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a3")
        if not self.goodEntry(self.qqZZshape_a4): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a4")
        if not self.goodEntry(self.qqZZshape_a5): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a5")
        if not self.goodEntry(self.qqZZshape_a6): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a6")
        if not self.goodEntry(self.qqZZshape_a7): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a7")
        if not self.goodEntry(self.qqZZshape_a8): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a8")
        if not self.goodEntry(self.qqZZshape_a9): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a9")
        if not self.goodEntry(self.qqZZshape_a10): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a10")
        if not self.goodEntry(self.qqZZshape_a11): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a11")
        if not self.goodEntry(self.qqZZshape_a12): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a12")
        if not self.goodEntry(self.qqZZshape_a13): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a13")

        if not self.goodEntry(self.ggZZshape_a0): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a0")
        if not self.goodEntry(self.ggZZshape_a1): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a1")
        if not self.goodEntry(self.ggZZshape_a2): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a2")
        if not self.goodEntry(self.ggZZshape_a3): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a3")
        if not self.goodEntry(self.ggZZshape_a4): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a4")
        if not self.goodEntry(self.ggZZshape_a5): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a5")
        if not self.goodEntry(self.ggZZshape_a6): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a6")
        if not self.goodEntry(self.ggZZshape_a7): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a7")
        if not self.goodEntry(self.ggZZshape_a8): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a8")
        if not self.goodEntry(self.ggZZshape_a9): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a9")

        if not self.goodEntry(self.zjetsShape_mean): raise RuntimeError, "{0} is not set.  Check inputs!".format("zjetsShape_mean")
        if not self.goodEntry(self.zjetsShape_sigma): raise RuntimeError, "{0} is not set.  Check inputs!".format("zjetsShape_sigma")
        
        if not self.goodEntry(self.zjetsKappaLow): raise RuntimeError, "{0} is not set.  Check inputs!".format("self.zjetsKappaLow")
        if not self.goodEntry(self.zjetsKappaHigh): raise RuntimeError, "{0} is not set.  Check inputs!".format("self.zjetsKappaHigh")

        if not self.goodEntry(self.qqZZ_lumi):  self.qqZZ_lumi = self.lumi
        if not self.goodEntry(self.ggZZ_lumi):  self.ggZZ_lumi = self.lumi
        if not self.goodEntry(self.zjets_lumi): self.zjets_lumi = self.lumi
        if not self.goodEntry(self.zbb_lumi):   self.zbb_lumi = self.lumi
        if not self.goodEntry(self.ttbar_lumi): self.ttbar_lumi = self.lumi
        
        ## Set dictionary entries to be passed to datacard class ##
        
        dict['decayChannel'] = self.decayChan
        dict['model'] = self.model
        dict['lumi'] = self.lumi
        dict['sqrts'] = self.sqrts

        dict['ggH'] = self.ggH_chan
        dict['qqH'] = self.qqH_chan
        dict['WH'] = self.WH_chan
        dict['ZH'] = self.ZH_chan
        dict['ttH'] = self.ttH_chan
        dict['qqZZ'] = self.qqZZ_chan
        dict['ggZZ'] = self.ggZZ_chan
        dict['zjets'] = self.zjets_chan
        dict['ttbar'] = self.ttbar_chan
        dict['zbb'] = self.zbb_chan

        dict['qqZZ_rate'] = self.qqZZ_rate
        dict['ggZZ_rate'] = self.ggZZ_rate 
        dict['zjets_rate'] = self.zjets_rate
        dict['ttbar_rate'] = self.ttbar_rate
        dict['zbb_rate'] = self.zbb_rate 

        dict['qqZZ_lumi'] = self.qqZZ_lumi
        dict['ggZZ_lumi'] = self.ggZZ_lumi 
        dict['zjets_lumi'] = self.zjets_lumi
        dict['ttbar_lumi'] = self.ttbar_lumi
        dict['zbb_lumi'] = self.zbb_lumi

        dict['n_CB_shape'] = self.n_CB_shape
        dict['alpha_CB_shape'] = self.alpha_CB_shape
        dict['mean_CB_shape'] = self.mean_CB_shape
        dict['sigma_CB_shape'] = self.sigma_CB_shape

        dict['sigEff_a1'] = self.sigeff_a1
        dict['sigEff_a2'] = self.sigeff_a2
        dict['sigEff_a3'] = self.sigeff_a3
        dict['sigEff_a4'] = self.sigeff_a4
        dict['sigEff_b1'] = self.sigeff_b1
        dict['sigEff_b2'] = self.sigeff_b2
        dict['sigEff_b3'] = self.sigeff_b3

        dict['qqZZshape_a0'] = self.qqZZshape_a0
        dict['qqZZshape_a1'] = self.qqZZshape_a1
        dict['qqZZshape_a2'] = self.qqZZshape_a2
        dict['qqZZshape_a3'] = self.qqZZshape_a3
        dict['qqZZshape_a4'] = self.qqZZshape_a4
        dict['qqZZshape_a5'] = self.qqZZshape_a5
        dict['qqZZshape_a6'] = self.qqZZshape_a6
        dict['qqZZshape_a7'] = self.qqZZshape_a7
        dict['qqZZshape_a8'] = self.qqZZshape_a8
        dict['qqZZshape_a9'] = self.qqZZshape_a9
        dict['qqZZshape_a10'] = self.qqZZshape_a10
        dict['qqZZshape_a11'] = self.qqZZshape_a11
        dict['qqZZshape_a12'] = self.qqZZshape_a12
        dict['qqZZshape_a13'] = self.qqZZshape_a13

        dict['ggZZshape_a0'] = self.ggZZshape_a0
        dict['ggZZshape_a1'] = self.ggZZshape_a1
        dict['ggZZshape_a2'] = self.ggZZshape_a2
        dict['ggZZshape_a3'] = self.ggZZshape_a3
        dict['ggZZshape_a4'] = self.ggZZshape_a4
        dict['ggZZshape_a5'] = self.ggZZshape_a5
        dict['ggZZshape_a6'] = self.ggZZshape_a6
        dict['ggZZshape_a7'] = self.ggZZshape_a7
        dict['ggZZshape_a8'] = self.ggZZshape_a8
        dict['ggZZshape_a9'] = self.ggZZshape_a9

        dict['zjetsShape_mean'] = self.zjetsShape_mean
        dict['zjetsShape_sigma'] = self.zjetsShape_sigma

        dict['zjetsKappaLow'] = self.zjetsKappaLow
        dict['zjetsKappaHigh'] = self.zjetsKappaHigh
        
        return dict

