#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from systematicsClass import *
from inputReader import *

## ------------------------------------
##  card and workspace class
## ------------------------------------

class datacardClass:

    def __init__(self):
    
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3    
        self.isFSR = True
        self.sigMorph = False
        self.bkgMorph = True

    def loadIncludes(self):
        
        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        #ROOT.gSystem.AddIncludePath("-Iinclude/FSR/")
        ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
        ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        

    # cross section filter for 7 TeV efficiency
    def csFilter(self,hmass):

        a = 80.85
        b = 50.42
        
        f = 0.5 + 0.5*erf( (hmass - a)/b )
        
        return f

    # cs x br function 
    def makeXsBrFunction(self,signalProc,rrvMH):
            
        procName = "ggH"
        if(signalProc == 1): procName = "ggH"
        if(signalProc == 2): procName = "qqH"
        if(signalProc == 3): procName = "WH"
        if(signalProc == 4): procName = "ZH"
        if(signalProc == 5): procName = "ttH"
        
        
        channelName = ""
        if (self.channel == self.ID_4mu): channelName = "4mu"
        elif (self.channel == self.ID_4e): channelName = "4e"
        elif (self.channel == self.ID_2e2mu): channelName = "2e2mu"
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)" 
        
        
        myCSWrhf = HiggsCSandWidth()
        
        histXsBr = ROOT.TH1F("hsmxsbr_{0}_{1}".format(procName,channelName),"", 981, 109.75, 600.25)
        
        for i in range(1,981):
            
            mHVal = histXsBr.GetBinCenter(i)
            BR = 0.0 
            if (self.channel == self.ID_2e2mu):
                BR = myCSWrhf.HiggsBR(13,mHVal)
            else:
                BR = myCSWrhf.HiggsBR(12,mHVal)
            histXsBr.SetBinContent(i, myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts) * BR)

            
            
        rdhname = "rdhXsBr_{0}_{1}_{2}".format(procName,self.channel,self.sqrts)
        rdhXsBr = RooDataHist(rdhname,rdhname, ROOT.RooArgList(rrvMH), histXsBr)  
        
        return rdhXsBr
    
    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement): 
            return trueVar
        else:
            return falseVar
    
    # main datacard and workspace function
    def makeCardsWorkspaces(self, theMH, theis2D, theOutputDir, theInputs):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = theMH
        self.lumi = theInputs['lumi']
        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.is2D = theis2D
        self.outputDir = theOutputDir
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3    
        self.isFSR = True
        self.sigMorph = False
        self.bkgMorph = True
        FactorizedShapes = False

        ## ---------------- SET PLOTTING STYLE ----------------  ## 
        ROOT.setTDRStyle(True)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)        

        self.bUseCBnoConvolution = False
        ForXSxBR = False

        myCSW = HiggsCSandWidth()
                
        self.widthHVal =  myCSW.HiggsWidth(0,self.mH)
        if(self.widthHVal < 0.12):
            self.bUseCBnoConvolution = True
    
        if(DEBUG): print "width: ",self.widthHVal
        
        self.windowVal = max( self.widthHVal, 1.0)
        lowside = 100.0
        if (self.mH >= 275):
            lowside = 180.0
        else:
            lowside = 100.0
        
        self.low_M = max( (self.mH - 20.*self.windowVal), lowside)
        self.high_M = min( (self.mH + 15.*self.windowVal), 800.0)
       
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"
	
        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##
    
        systematics = systematicsClass( self.sqrts, self.channel, self.mH, False, self.isFSR, "SM")
        systematics_forXSxBR = systematicsClass( self.sqrts, self.channel, self.mH, True, self.isFSR, "SM")

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##
    
        bins = 1000
        if(self.bUseCBnoConvolution): bins = 200

        CMS_zz4l_mass = ROOT.RooRealVar("CMS_zz4l_mass","CMS_zz4l_mass",self.low_M,self.high_M)
        CMS_zz4l_mass.setBins(bins,"fft") 

        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        self.LUMI.setConstant(True)
    
        self.MH = ROOT.RooRealVar("MH","MH",self.mH)
        self.MH.setConstant(True)
    
        n_CB_d = 0.0
        alpha_CB_d = 0.0
        mean_CB_d = 0.0
        sigma_CB_d = 0.0
        mean_BW_d = self.mH
        
        rdhXsBrFuncV_1 = self.makeXsBrFunction(1,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ggH",self.channel,self.sqrts)
        rhfXsBrFuncV_1 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_1, 1)
        
        rdhXsBrFuncV_2 = self.makeXsBrFunction(2,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("VBF",self.channel,self.sqrts)
        rhfXsBrFuncV_2 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_2, 1)
        
        rdhXsBrFuncV_3 = self.makeXsBrFunction(3,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("WH",self.channel,self.sqrts)
        rhfXsBrFuncV_3 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_3, 1)
        
        rdhXsBrFuncV_4 = self.makeXsBrFunction(4,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ZH",self.channel,self.sqrts)
        rhfXsBrFuncV_4 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_4, 1)
        
        rdhXsBrFuncV_5 = self.makeXsBrFunction(5,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ttH",self.channel,self.sqrts)
        rhfXsBrFuncV_5 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_5, 1)

    
        ## -------- Variable Definitions -------- ##
    
        name = "CMS_zz4l_mean_e_sig"
        CMS_zz4l_mean_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_sigma_e_sig"
        CMS_zz4l_sigma_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
        
        name = "CMS_zz4l_mean_m_sig"
        CMS_zz4l_mean_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_sigma_m_sig"
        CMS_zz4l_sigma_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
        
        name = "CMS_zz4l_alpha_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha = ROOT.RooRealVar(name,"CMS_zz4l_alpha",1.,-10.,10.)
        name = "CMS_zz4l_n_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n = ROOT.RooRealVar(name,"CMS_zz4l_n",2.,-10.,10.)
        name = "CMS_zz4l_mean_BW_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_BW = ROOT.RooRealVar(name,"CMS_zz4l_mean_BW",self.mH,self.low_M,self.high_M)
        name = "CMS_zz4l_gamma_BW_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_gamma_BW = ROOT.RooRealVar(name,"CMS_zz4l_gamma_BW",self.widthHVal)
        name = "CMS_zz4l_widthScale_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_widthScale = ROOT.RooRealVar(name,"CMS_zz4l_widthScale",1.0)
        one = ROOT.RooRealVar("one","one",1.0)
        one.setConstant(True)
    
        CMS_zz4l_mean_BW.setVal( mean_BW_d )
        CMS_zz4l_gamma_BW.setVal( self.widthHVal )
        CMS_zz4l_mean_e_sig.setVal(0)
        CMS_zz4l_sigma_e_sig.setVal(0)
        CMS_zz4l_mean_m_sig.setVal(0)
        CMS_zz4l_sigma_m_sig.setVal(0)
        CMS_zz4l_alpha.setVal(0)
        CMS_zz4l_n.setVal(0)
    
        CMS_zz4l_widthScale.setConstant(True)
        CMS_zz4l_alpha.setConstant(True)
        CMS_zz4l_mean_BW.setConstant(True)
        CMS_zz4l_gamma_BW.setConstant(True)

        ## -------------------- RooFormulaVar's -------------------- ##
        
        name = "CMS_zz4l_n_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
        name = "CMS_zz4l_alpha_{0:.0f}_centralValue".format(self.channel)
        rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape'], ROOT.RooArgList(one))
        name = "CMS_zz4l_mean_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig))
        name = "CMS_zz4l_sigma_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))

        CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"@0+@1", ROOT.RooArgList(rfv_mean_CB, self.MH))

        ## --------------------- SHAPE FUNCTIONS ---------------------- ##
    
        signalCB_ggH = ROOT.RooCBShape("signalCB_ggH","signalCB_ggH",CMS_zz4l_mass, self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB, self.bUseCBnoConvolution) ,rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB)
        signalBW_ggH = ROOT.RooRelBWUFParam("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ggH =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH,signalCB_ggH, 2)
        
        signalCB_VBF = ROOT.RooCBShape("signalCB_VBF","signalCB_VBF",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB)
        signalBW_VBF = ROOT.RooRelBWUFParam("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_VBF = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF,signalCB_VBF, 2)
        
        signalCB_WH = ROOT.RooCBShape("signalCB_WH","signalCB_WH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB)
        signalBW_WH = ROOT.RooRelBWUFParam("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale);
        sig_WH = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH,signalCB_WH, 2);
        
        signalCB_ZH = ROOT.RooCBShape("signalCB_ZH","signalCB_ZH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB)
        signalBW_ZH = ROOT.RooRelBWUFParam("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ZH = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH,signalCB_ZH, 2)
        
        signalCB_ttH = ROOT.RooCBShape("signalCB_ttH","signalCB_ttH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),rfv_sigma_CB,rfv_alpha_CB,rfv_n_CB)
        signalBW_ttH = ROOT.RooRelBWUFParam("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ttH = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH,signalCB_ttH, 2) 
        
        ## Buffer fraction for cyclical behavior
        sig_ggH.setBufferFraction(0.2)
        sig_VBF.setBufferFraction(0.2)
        sig_WH.setBufferFraction(0.2)
        sig_ZH.setBufferFraction(0.2)
        sig_ttH.setBufferFraction(0.2)
        
        ## --------------------------- MELA 2D PDFS ------------------------- ##

        D = ROOT.RooRealVar("melaLD","melaLD",0,1);
    
        templateSigName = "templates2D/Dsignal_{0}.root".format(self.appendName)
        
        sigTempFile = ROOT.TFile(templateSigName)
        sigTemplate = sigTempFile.Get("h_mzzD")
        sigTemplate_Up = sigTempFile.Get("h_mzzD_up")
        sigTemplate_Down = sigTempFile.Get("h_mzzD_dn")
        
        TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate)
        TemplateName = "sigTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Up)
        TemplateName = "sigTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Down)
        
        TemplateName = "sigTemplatePdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ggH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ggH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_VBF = ROOT.RooHistPdf(TemplateName,TemplateName,RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_VBF_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_VBF_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_VBF_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_VBF_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_WH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_WH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_WH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_WH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_WH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ZH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ZH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ZH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ttH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ttH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ttH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)

        funcList_ggH = ROOT.RooArgList()  
        funcList_VBF = ROOT.RooArgList()
        funcList_WH  = ROOT.RooArgList()
        funcList_ZH  = ROOT.RooArgList()
        funcList_ttH = ROOT.RooArgList()

        if(self.sigMorph):
            
            funcList_ggH.add(sigTemplatePdf_ggH)
            funcList_ggH.add(sigTemplatePdf_ggH_Up)
            funcList_ggH.add(sigTemplatePdf_ggH_Down)  
            
            funcList_VBF.add(sigTemplatePdf_VBF)
            funcList_VBF.add(sigTemplatePdf_VBF_Up)
            funcList_VBF.add(sigTemplatePdf_VBF_Down)  
            
            funcList_WH.add(sigTemplatePdf_WH)
            funcList_WH.add(sigTemplatePdf_WH_Up)
            funcList_WH.add(sigTemplatePdf_WH_Down)  
            
            funcList_ZH.add(sigTemplatePdf_ZH)
            funcList_ZH.add(sigTemplatePdf_ZH_Up)
            funcList_ZH.add(sigTemplatePdf_ZH_Down)  
            
            funcList_ttH.add(sigTemplatePdf_ttH)
            funcList_ttH.add(sigTemplatePdf_ttH_Up)
            funcList_ttH.add(sigTemplatePdf_ttH_Down)  
            
        else:
            
            funcList_ggH.add(sigTemplatePdf_ggH)
            funcList_VBF.add(sigTemplatePdf_VBF)
            funcList_WH.add(sigTemplatePdf_WH)
            funcList_ZH.add(sigTemplatePdf_ZH)
            funcList_ttH.add(sigTemplatePdf_ttH)
      
    
        morphSigVarName = "CMS_zz4l_sigMELA_{0:.0f}".format(self.channel)
        alphaMorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        if(self.sigMorph): alphaMorphSig.setConstant(False)
        else: alphaMorphSig.setConstant(True)
        
        morphVarListSig = ROOT.RooArgList()
    
        if(self.sigMorph): morphVarListSig.add(alphaMorphSig);  ## just one morphing for all signal processes (fully correlated systematics)
        
        TemplateName = "sigTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ggH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_VBF = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_VBF,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_WH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_WH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_ZH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ZH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ttH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_ttH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ttH,morphVarListSig,1.0,1)
    

    
        sig2d_ggH = ROOT.RooProdPdf("sig2d_ggH","sig2d_ggH",ROOT.RooArgSet(sig_ggH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(D)))
        sig2d_VBF = ROOT.RooProdPdf("sig2d_VBF","sig2d_VBF",ROOT.RooArgSet(sig_VBF),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(D)))
        sig2d_WH = ROOT.RooProdPdf("sig2d_WH","sig2d_WH",ROOT.RooArgSet(sig_WH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(D)))
        sig2d_ZH = ROOT.RooProdPdf("sig2d_ZH","sig2d_ZH",ROOT.RooArgSet(sig_ZH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(D)))
        sig2d_ttH = ROOT.RooProdPdf("sig2d_ttH","sig2d_ttH",ROOT.RooArgSet(sig_ttH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(D)))
        
        sigCB2d_ggH = ROOT.RooProdPdf("sigCB2d_ggH","sigCB2d_ggH",ROOT.RooArgSet(signalCB_ggH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(D)))
        sigCB2d_VBF = ROOT.RooProdPdf("sigCB2d_VBF","sigCB2d_VBF",ROOT.RooArgSet(signalCB_VBF),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(D)))
        sigCB2d_WH = ROOT.RooProdPdf("sigCB2d_WH","sigCB2d_WH",ROOT.RooArgSet(signalCB_WH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(D)))
        sigCB2d_ZH = ROOT.RooProdPdf("sigCB2d_ZH","sigCB2d_ZH",ROOT.RooArgSet(signalCB_ZH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(D)))
        sigCB2d_ttH = ROOT.RooProdPdf("sigCB2d_ttH","sigCB2d_ttH",ROOT.RooArgSet(signalCB_ttH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(D)))
        


        ## -------------------------- BACKGROUND SHAPES ---------------------------------- ##
    
        ## qqZZ contribution
        name = "CMS_qqzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a0 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a0",115.3,0.,200.)
        name = "CMS_qqzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a1 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a1",21.96,0.,200.)
        name = "CMS_qqzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a2 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a2",122.8,0.,200.)
        name = "CMS_qqzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a3 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a3",0.03479,0.,1.)
        name = "CMS_qqzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a4 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a4",185.5,0.,200.)
        name = "CMS_qqzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a5 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a5",12.67,0.,200.)
        name = "CMS_qqzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a6 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a6",34.81,0.,100.)
        name = "CMS_qqzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a7 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a7",0.1393,0.,1.)
        name = "CMS_qqzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a8 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a8",66.,0.,200.)
        name = "CMS_qqzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a9 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a9",0.07191,0.,1.)
        name = "CMS_qqzzbkg_a10_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a10 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a10",94.11,0.,200.)
        name = "CMS_qqzzbkg_a11_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a11 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a11",-5.111,-100.,100.)
        name = "CMS_qqzzbkg_a12_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a12 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a12",4834,0.,10000.)
        name = "CMS_qqzzbkg_a13_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a13 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a13",0.2543,0.,1.)
        

        if (DEBUG) :
            print "qqZZshape_a0 = ",theInputs['qqZZshape_a0']
            print "qqZZshape_a1 = ",theInputs['qqZZshape_a1']
            print "qqZZshape_a2 = ",theInputs['qqZZshape_a2']
            print "qqZZshape_a3 = ",theInputs['qqZZshape_a3']
            print "qqZZshape_a4 = ",theInputs['qqZZshape_a4']
            print "qqZZshape_a5 = ",theInputs['qqZZshape_a5']
            print "qqZZshape_a6 = ",theInputs['qqZZshape_a6']
            print "qqZZshape_a7 = ",theInputs['qqZZshape_a7']
            print "qqZZshape_a8 = ",theInputs['qqZZshape_a8']
            print "qqZZshape_a9 = ",theInputs['qqZZshape_a9']
            print "qqZZshape_a10 = ",theInputs['qqZZshape_a10']
            print "qqZZshape_a11 = ",theInputs['qqZZshape_a11']
            print "qqZZshape_a12 = ",theInputs['qqZZshape_a12']
            print "qqZZshape_a13 = ",theInputs['qqZZshape_a13']

        
        CMS_qqzzbkg_a0.setVal(theInputs['qqZZshape_a0'])
        CMS_qqzzbkg_a1.setVal(theInputs['qqZZshape_a1'])
        CMS_qqzzbkg_a2.setVal(theInputs['qqZZshape_a2'])
        CMS_qqzzbkg_a3.setVal(theInputs['qqZZshape_a3'])
        CMS_qqzzbkg_a4.setVal(theInputs['qqZZshape_a4'])
        CMS_qqzzbkg_a5.setVal(theInputs['qqZZshape_a5'])
        CMS_qqzzbkg_a6.setVal(theInputs['qqZZshape_a6'])
        CMS_qqzzbkg_a7.setVal(theInputs['qqZZshape_a7'])
        CMS_qqzzbkg_a8.setVal(theInputs['qqZZshape_a8'])
        CMS_qqzzbkg_a9.setVal(theInputs['qqZZshape_a9'])
        CMS_qqzzbkg_a10.setVal(theInputs['qqZZshape_a10'])
        CMS_qqzzbkg_a11.setVal(theInputs['qqZZshape_a11'])
        CMS_qqzzbkg_a12.setVal(theInputs['qqZZshape_a12'])
        CMS_qqzzbkg_a13.setVal(theInputs['qqZZshape_a13'])
        
        CMS_qqzzbkg_a0.setConstant(True)
        CMS_qqzzbkg_a1.setConstant(True)
        CMS_qqzzbkg_a2.setConstant(True)
        CMS_qqzzbkg_a3.setConstant(True)
        CMS_qqzzbkg_a4.setConstant(True)
        CMS_qqzzbkg_a5.setConstant(True)
        CMS_qqzzbkg_a6.setConstant(True)
        CMS_qqzzbkg_a7.setConstant(True)
        CMS_qqzzbkg_a8.setConstant(True)
        CMS_qqzzbkg_a9.setConstant(True)
        CMS_qqzzbkg_a10.setConstant(True)
        CMS_qqzzbkg_a11.setConstant(True)
        CMS_qqzzbkg_a12.setConstant(True)
        CMS_qqzzbkg_a13.setConstant(True)
        
        bkg_qqzz = ROOT.RooqqZZPdf_v2("bkg_qqzz","bkg_qqzz",CMS_zz4l_mass,CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,CMS_qqzzbkg_a9,CMS_qqzzbkg_a10,CMS_qqzzbkg_a11,CMS_qqzzbkg_a12,CMS_qqzzbkg_a13);
        
        ## ggZZ contribution
        name = "CMS_ggzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a0 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a0",115.3,0.,200.)
        name = "CMS_ggzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a1 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a1",21.96,0.,200.)
        name = "CMS_ggzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a2 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a2",122.8,0.,200.)
        name = "CMS_ggzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a3 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a3",0.03479,0.,1.)
        name = "CMS_ggzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        CMS_ggzzbkg_a4 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a4",185.5,0.,200.)
        name = "CMS_ggzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a5 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a5",12.67,0.,200.)
        name = "CMS_ggzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a6 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a6",34.81,0.,100.)
        name = "CMS_ggzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a7 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a7",0.1393,0.,1.)
        name = "CMS_ggzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a8 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a8",66.,0.,200.)
        name = "CMS_ggzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        CMS_ggzzbkg_a9 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a9",0.07191,0.,1.)
        
        CMS_ggzzbkg_a0.setVal(theInputs['ggZZshape_a0'])
        CMS_ggzzbkg_a1.setVal(theInputs['ggZZshape_a1'])
        CMS_ggzzbkg_a2.setVal(theInputs['ggZZshape_a2'])
        CMS_ggzzbkg_a3.setVal(theInputs['ggZZshape_a3'])
        CMS_ggzzbkg_a4.setVal(theInputs['ggZZshape_a4'])
        CMS_ggzzbkg_a5.setVal(theInputs['ggZZshape_a5'])
        CMS_ggzzbkg_a6.setVal(theInputs['ggZZshape_a6'])
        CMS_ggzzbkg_a7.setVal(theInputs['ggZZshape_a7'])
        CMS_ggzzbkg_a8.setVal(theInputs['ggZZshape_a8'])
        CMS_ggzzbkg_a9.setVal(theInputs['ggZZshape_a9'])
        
        CMS_ggzzbkg_a0.setConstant(True)
        CMS_ggzzbkg_a1.setConstant(True)
        CMS_ggzzbkg_a2.setConstant(True)
        CMS_ggzzbkg_a3.setConstant(True)
        CMS_ggzzbkg_a4.setConstant(True)
        CMS_ggzzbkg_a5.setConstant(True)
        CMS_ggzzbkg_a6.setConstant(True)
        CMS_ggzzbkg_a7.setConstant(True)
        CMS_ggzzbkg_a8.setConstant(True)
        CMS_ggzzbkg_a9.setConstant(True)

        if (DEBUG) :
            print "ggZZshape_a0 = ",theInputs['ggZZshape_a0']
            print "ggZZshape_a1 = ",theInputs['ggZZshape_a1']
            print "ggZZshape_a2 = ",theInputs['ggZZshape_a2']
            print "ggZZshape_a3 = ",theInputs['ggZZshape_a3']
            print "ggZZshape_a4 = ",theInputs['ggZZshape_a4']
            print "ggZZshape_a5 = ",theInputs['ggZZshape_a5']
            print "ggZZshape_a6 = ",theInputs['ggZZshape_a6']
            print "ggZZshape_a7 = ",theInputs['ggZZshape_a7']
            print "ggZZshape_a8 = ",theInputs['ggZZshape_a8']
            print "ggZZshape_a9 = ",theInputs['ggZZshape_a9']
                   
        
        bkg_ggzz = ROOT.RooggZZPdf_v2("bkg_ggzz","bkg_ggzz",CMS_zz4l_mass,CMS_ggzzbkg_a0,CMS_ggzzbkg_a1,CMS_ggzzbkg_a2,CMS_ggzzbkg_a3,CMS_ggzzbkg_a4,CMS_ggzzbkg_a5,CMS_ggzzbkg_a6,CMS_ggzzbkg_a7,CMS_ggzzbkg_a8,CMS_ggzzbkg_a9)
    
        ## Reducible backgrounds
        val_meanL = float(theInputs['zjetsShape_mean'])
        val_sigmaL = float(theInputs['zjetsShape_sigma'])
 
        name = "mlZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",val_meanL)
        name = "slZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",val_sigmaL)
        bkg_zjets = ROOT.RooLandau("bkg_zjets","bkg_zjets",CMS_zz4l_mass,mlZjet,slZjet) 

        ## ----------------- 2D BACKGROUND SHAPES --------------- ##
        
        templateBkgName = "templates2D/Dbackground_qqZZ_{0}.root".format(self.appendName)
        bkgTempFile = ROOT.TFile(templateBkgName)
        bkgTemplate = bkgTempFile.Get("h_mzzD")
        bkgTemplate_Up = bkgTempFile.Get("h_mzzD_up")
        bkgTemplate_Down = bkgTempFile.Get("h_mzzD_dn")
        
        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),bkgTemplate)
        TemplateName = "bkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),bkgTemplate_Up)
        TemplateName = "bkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),bkgTemplate_Down)
        
        templateggBkgName = "templates2D/Dbackground_ggZZ_{0}.root".format(self.appendName)
        ggbkgTempFile = ROOT.TFile(templateggBkgName)
        ggbkgTemplate = ggbkgTempFile.Get("h_mzzD")
        TemplateName = "ggbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        ggbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),ggbkgTemplate)
        
        TemplateName = "bkgTemplatePdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_qqzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist)
        TemplateName = "bkgTemplatePdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_ggzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),ggbkgTempDataHist)
        TemplateName = "bkgTemplatePdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist)
        TemplateName,"bkgTemplatePdf_zjets_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist_Up)
        TemplateName = "bkgTemplatePdf_zjets_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist_Down)
        
        funcList_zjets = ROOT.RooArgList()  
        morphBkgVarName = "CMS_zz4l_bkgMELA"    
        alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()
        
        if(self.bkgMorph):
            funcList_zjets.add(bkgTemplatePdf_zjets)
            funcList_zjets.add(bkgTemplatePdf_zjets_Up)
            funcList_zjets.add(bkgTemplatePdf_zjets_Down)  
            alphaMorphBkg.setConstant(False)
            morphVarListBkg.add(alphaMorphBkg)  
        else:
            funcList_zjets.add(bkgTemplatePdf_zjets)
            alphaMorphBkg.setConstant(True)


        TemplateName = "bkgTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,ROOT.RooArgList(bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,ROOT.RooArgList(bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_zjets,morphVarListBkg,1.0,1)

        bkg2d_qqzz = ROOT.RooProdPdf("bkg2d_qqzz","bkg2d_qqzz",ROOT.RooArgSet(bkg_qqzz),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_qqzz),ROOT.RooArgSet(D)))
        bkg2d_ggzz = ROOT.RooProdPdf("bkg2d_ggzz","bkg2d_ggzz",ROOT.RooArgSet(bkg_ggzz),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_ggzz),ROOT.RooArgSet(D)))
        bkg2d_zjets = ROOT.RooProdPdf("bkg2d_zjets","bkg2d_zjets",ROOT.RooArgSet(bkg_zjets),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_zjets),ROOT.RooArgSet(D)))
      
        ## ----------------------- PLOTS FOR SANITY CHECKS -------------------------- ##
        
        czz = ROOT.TCanvas( "czz", "czz", 750, 700 )
        czz.cd()
        zzframe_s = CMS_zz4l_mass.frame(45);
        if self.bUseCBnoConvolution: super(ROOT.RooCBShape,signalCB_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        else: super(ROOT.RooFFTConvPdf,sig_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        super(ROOT.RooqqZZPdf_v2,bkg_qqzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(4) )
        super(ROOT.RooggZZPdf_v2,bkg_ggzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(6) )
        super(ROOT.RooLandau,bkg_zjets).plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(6) )
        zzframe_s.Draw()
        figName = "{0}/figs/mzz_{1}_{2}.png".format(self.outputDir, self.mH, self.appendName)
        czz.SaveAs(figName)
        del czz
        
        ## ------------------- LUMI -------------------- ##
        
        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  
        
        ## ----------------------- SIGNAL RATES ----------------------- ##
        
        CMS_zz4l_mass.setRange("shape",self.low_M,self.high_M)
        
        fr_low_M = self.low_M
        fr_high_M = self.high_M        
        if (self.mH >= 450): 
            fr_low_M = 100.
            fr_high_M = 1000.

        CMS_zz4l_mass.setRange("fullrangesignal",fr_low_M,fr_high_M)
        CMS_zz4l_mass.setRange("fullrange",100.,1000.)
        
        rfvCsFilter = RooFormulaVar()
        filterName = "cmshzz4l_csFilter_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.sqrts == 7): 
            rfvCsFilter = ROOT.RooFormulaVar(filterName,"0.5+0.5*TMath::Erf((@0 - 80.85)/50.42)", ROOT.RooArgList(self.MH) )
        else:
            rfvCsFilter = ROOT.RooFormulaVar(filterName,"@0",ROOT.RooArgList(one))

        if(DEBUG):
            print "@@@@@@@ rfvCsFilter = ",rfvCsFilter.getVal()

        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        rrva1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        rrva2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        rrva3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a3'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        rrva4 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a4'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        rrvb1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        rrvb2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        rrvb3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b3'])

        if(DEBUG):
            print "sigEff_a1 = ",theInputs['sigEff_a1']
            print "sigEff_a2 = ",theInputs['sigEff_a2']
            print "sigEff_a3 = ",theInputs['sigEff_a3']
            print "sigEff_a4 = ",theInputs['sigEff_a4']
            print "sigEff_b1 = ",theInputs['sigEff_b1']
            print "sigEff_b2 = ",theInputs['sigEff_b2']
            print "sigEff_b3 = ",theInputs['sigEff_b3']
           
    
        sigEffName = "hzz4lsignaleff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        rfvSigEff = ROOT.RooFormulaVar(sigEffName,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)", ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH))
        ## following printout is needed ,  dont remove it
        print " @@@@@@@@ sigeff ",rfvSigEff.getVal()
    
        CS_ggH = myCSW.HiggsCS(1,self.mH,self.sqrts)
        CS_VBF = myCSW.HiggsCS(2,self.mH,self.sqrts)
        CS_WH = myCSW.HiggsCS(3,self.mH,self.sqrts)
        CS_ZH = myCSW.HiggsCS(4,self.mH,self.sqrts)
        CS_ttH = myCSW.HiggsCS(5,self.mH,self.sqrts)
    
        BRH2e2mu = myCSW.HiggsBR(13,self.mH)
        BRH4mu = myCSW.HiggsBR(12,self.mH)
        BRH4e = myCSW.HiggsBR(12,self.mH)
        BR = 0.0
        if( self.channel == self.ID_4mu ): BR = BRH4mu
        if( self.channel == self.ID_4e ): BR = BRH4e
        if( self.channel == self.ID_2e2mu ): BR = BRH2e2mu
    
        sigEfficiency = rfvSigEff.getVal()

        if(DEBUG):
            print "CS_ggH: ",CS_ggH,", CS_VBF: ",CS_VBF,", CS_WH: ",CS_WH,", CS_ZH: ",CS_ZH
            print ", CS_ttH: ",CS_ttH,", BRH2e2mu: ",BRH2e2mu,", BRH4mu: ",BRH4mu,", BRH4e: ",BRH4e

        csCorrection = 1.0
        if(self.sqrts == 7): csCorrection = self.csFilter(self.mH)

        ## SIG YIELDS
        sigRate_ggH = csCorrection*CS_ggH*BR*sigEfficiency*1000.*self.lumi
        sigRate_VBF = csCorrection*CS_VBF*BR*sigEfficiency*1000.*self.lumi
        sigRate_WH = csCorrection*CS_WH*BR*sigEfficiency*1000.*self.lumi
        sigRate_ZH = csCorrection*CS_ZH*BR*sigEfficiency*1000.*self.lumi
        sigRate_ttH = csCorrection*CS_ttH*BR*sigEfficiency*1000.*self.lumi
       
        tmpNormSigNoConv = signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        tmpNormSigConv = sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()

        normalizationSignal = self.getVariable(tmpNormSigNoConv,tmpNormSigConv,self.bUseCBnoConvolution)

        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()

        
        sclFactorSig_ggH = sigRate_ggH/normalizationSignal
        sclFactorSig_VBF = sigRate_VBF/normalizationSignal
        sclFactorSig_WH = sigRate_WH/normalizationSignal
        sclFactorSig_ZH = sigRate_ZH/normalizationSignal
        sclFactorSig_ttH = sigRate_ttH/normalizationSignal
    
        integral_ggH = self.getVariable(signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        integral_VBF = self.getVariable(signalCB_VBF.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_VBF.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        integral_WH = self.getVariable(signalCB_WH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_WH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        integral_ZH = self.getVariable(signalCB_ZH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ZH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        integral_ttH = self.getVariable(signalCB_ttH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ttH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        
        sigRate_ggH_Shape = sclFactorSig_ggH*integral_ggH
        sigRate_VBF_Shape = sclFactorSig_VBF*integral_VBF
        sigRate_WH_Shape = sclFactorSig_WH*integral_WH
        sigRate_ZH_Shape = sclFactorSig_ZH*integral_ZH
        sigRate_ttH_Shape = sclFactorSig_ttH*integral_ttH
        
        
        normSigName = "cmshzz4l_normalizationSignal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, self.getVariable(signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),self.bUseCBnoConvolution))
        rrvNormSig.setConstant(True)

        rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_1))
        
        
        rfvSigRate_VBF = ROOT.RooFormulaVar("VBF_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_VBF.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_VBF.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_2))
                         

        rfvSigRate_WH = ROOT.RooFormulaVar("WH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_WH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_WH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_3))
                         

        rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_ZH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ZH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_4))
                         

        rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_ttH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ttH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_5))
                         


        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal()
        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal()
        print " @@@@@@@ norm sig = ",rrvNormSig.getVal()
        print " @@@@@@@ rfvSigRate_ggH = ",rfvSigRate_ggH.getVal()
        print " sigRate_ggH_Shape=",sigRate_ggH_Shape
        print " @@@@@@@ rfvSigRate_VBF = ",rfvSigRate_VBF.getVal()
        print " sigRate_VBF_Shape=",sigRate_VBF_Shape
        print " @@@@@@@ rfvSigRate_WH = ",rfvSigRate_WH.getVal()
        print " sigRate_WH_Shape=",sigRate_WH_Shape
        print " @@@@@@@ rfvSigRate_ZH = ",rfvSigRate_ZH.getVal()
        print " sigRate_ZH_Shape=",sigRate_ZH_Shape
        print " @@@@@@@ rfvSigRate_ttH = ",rfvSigRate_ttH.getVal()
        print " sigRate_ttH_Shape=",sigRate_ttH_Shape
        
        ## SET RATES TO 1 
        ## DC RATES WILL BE MULTIPLIED
        ## BY RATES IMPORTED TO WS
        sigRate_ggH_Shape = 1
        sigRate_VBF_Shape = 1
        sigRate_WH_Shape = 1
        sigRate_ZH_Shape = 1
        sigRate_ttH_Shape = 1

             
        ## ----------------------- BACKGROUND RATES ----------------------- ##

        ## rates per lumi for scaling
        bkgRate_qqzz = theInputs['qqZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_ggzz = theInputs['ggZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_zjets = theInputs['zjets_rate']/theInputs['zjets_lumi']
        
        ## Get Normalizations
        normalizationBackground_qqzz = bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_ggzz = bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_zjets = bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        
        sclFactorBkg_qqzz = self.lumi*bkgRate_qqzz/normalizationBackground_qqzz
        sclFactorBkg_ggzz = self.lumi*bkgRate_ggzz/normalizationBackground_ggzz
        sclFactorBkg_zjets = self.lumi*bkgRate_zjets/normalizationBackground_zjets
               
        bkgRate_qqzz_Shape = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_ggzz_Shape = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_zjets_Shape = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        
        if(DEBUG):
            print "Shape signal rate: ",sigRate_ggH_Shape,", background rate: ",bkgRate_qqzz_Shape,", ",bkgRate_zjets_Shape," in ",low_M," - ",high_M
            CMS_zz4l_mass.setRange("lowmassregion",100.,160.)
            bkgRate_qqzz_lowmassregion = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_ggzz_lowmassregion = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_zjets_lowmassregion = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            lowmassyield = bkgRate_qqzz_lowmassregion + bkgRate_ggzz_lowmassregion + bkgRate_zjets_lowmassregion
            print "low mass yield: ",lowmassyield
        
        ## --------------------------- DATASET --------------------------- ##

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs" 
        dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.lumi)
        if (DEBUG): print dataFileName," ",dataTreeName 
        data_obs_file = ROOT.TFile(dataFileName)

        print data_obs_file.Get(dataTreeName)
        
        if not (data_obs_file.Get(dataTreeName)):
            print "File, \"",dataFileName,"\", or tree, \"",dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
        
        data_obs_tree = data_obs_file.Get(dataTreeName)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{0}".format(self.appendName)
        
        if not self.is2D:
            data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass))
        else:
            data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,D))
            
            
        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.mH
        if ( math.fabs(math.floor(tmpMH)-self.mH) == 0.5): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5

        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWS2 = ""
        name_ShapeWSXSBR = ""
        
        if (endsInP5): name_Shape = "{0}/HCG/{1}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        if (endsInP5): name_ShapeWS = "{0}/HCG/{1}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        if (endsInP5): name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)

        if(DEBUG): print name_Shape,"  ",name_ShapeWS2
        
        w = ROOT.RooWorkspace("w","w")
        
        w.importClassCode(RooqqZZPdf_v2.Class(),True)
        w.importClassCode(RooggZZPdf_v2.Class(),True)
        w.importClassCode(RooRelBWUFParam.Class(),True)
        w.importClassCode(RooCBShape.Class(),True)
        w.importClassCode(RooFormulaVar.Class(),True)

        if( FactorizedShapes ):
            if( self.channel == self.ID_4mu ):
                w.importClassCode(RooFourMuMassShapePdf2.Class(),True)
                w.importClassCode(RooFourMuMassRes.Class(),True)
            elif( self.channel == self.ID_4e ):
                w.importClassCode(RooFourEMassShapePdf2.Class(),True)
                w.importClassCode(RooFourEMassRes.Class(),True)
            elif( self.channel == self.ID_2e2mu ):
                w.importClassCode(RooTwoETwoMuMassShapePdf2.Class(),True)
                w.importClassCode(RooTwoETwoMuMassRes.Class(),True)
            
                
                
        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?
    
        if(self.bUseCBnoConvolution) :
            if not self.is2D :
                signalCB_ggH.SetNameTitle("ggH","ggH")
                signalCB_VBF.SetNameTitle("qqH","qqH")
                signalCB_WH.SetNameTitle("WH","WH")
                signalCB_ZH.SetNameTitle("ZH","ZH")
                signalCB_ttH.SetNameTitle("ttH","ttH")
                
                getattr(w,'import')(signalCB_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(signalCB_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(signalCB_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(signalCB_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(signalCB_ttH, ROOT.RooFit.RecycleConflictNodes())
                
            else:
                sigCB2d_ggH.SetNameTitle("ggH","ggH")
                sigCB2d_VBF.SetNameTitle("qqH","qqH")
                sigCB2d_WH.SetNameTitle("WH","WH")
                sigCB2d_ZH.SetNameTitle("ZH","ZH")
                sigCB2d_ttH.SetNameTitle("ttH","ttH")
                
                getattr(w,'import')(sigCB2d_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigCB2d_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigCB2d_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigCB2d_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigCB2d_ttH, ROOT.RooFit.RecycleConflictNodes())
        else:
                
            if not self.is2D :
                sig_ggH.SetNameTitle("ggH","ggH")
                sig_VBF.SetNameTitle("qqH","qqH")
                sig_WH.SetNameTitle("WH","WH")
                sig_ZH.SetNameTitle("ZH","ZH")
                sig_ttH.SetNameTitle("ttH","ttH")
                
                getattr(w,'import')(sig_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig_ttH, ROOT.RooFit.RecycleConflictNodes())
                    
            else:
                sig2d_ggH.SetNameTitle("ggH","ggH")
                sig2d_VBF.SetNameTitle("qqH","qqH")
                sig2d_WH.SetNameTitle("WH","WH")
                sig2d_ZH.SetNameTitle("ZH","ZH")
                sig2d_ttH.SetNameTitle("ttH","ttH")     
                
                getattr(w,'import')(sig2d_ggH,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig2d_VBF,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig2d_WH,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig2d_ZH,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig2d_ttH,ROOT.RooFit.RecycleConflictNodes())            
        
        
        if not self.is2D :
            getattr(w,'import')(bkg_qqzz, ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkg_ggzz, ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkg_zjets, ROOT.RooFit.RecycleConflictNodes())
            
        else:
            getattr(w,'import')(bkg2d_qqzz,ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkg2d_ggzz,ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkg2d_zjets,ROOT.RooFit.RecycleConflictNodes())
        
        getattr(w,'import')(rfvSigRate_ggH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_VBF, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_WH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_ZH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_ttH, ROOT.RooFit.RecycleConflictNodes())

        w.writeToFile(name_ShapeWS)
        w.writeToFile(name_ShapeWSXSBR)
        
        ## --------------------------- DATACARDS -------------------------- ##

        rates_Shape = array('d',[sigRate_ggH_Shape,sigRate_VBF_Shape,sigRate_WH_Shape,sigRate_ZH_Shape,sigRate_ttH_Shape,bkgRate_qqzz_Shape,bkgRate_ggzz_Shape,bkgRate_zjets_Shape])

        systematics.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape, bkgRate_zjets_Shape)
        systematics_forXSxBR.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape,bkgRate_zjets_Shape)

        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo, name_ShapeWS2, rates_Shape, data_obs.numEntries(), self.is2D )
        systematics.WriteSystematics(fo, self.bkgMorph, self.sigMorph)
        systematics.WriteShapeSystematics(fo)
        fo.close()
        
        ## forXSxBR
        if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)	
        else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)
            
        fo = open( name_Shape, "wb" )
        self.WriteDatacard(fo, name_ShapeWS2, rates_Shape, data_obs.numEntries(), self.is2D )
        systematics_forXSxBR.WriteSystematics(fo, self.bkgMorph, self.sigMorph)
        systematics_forXSxBR.WriteShapeSystematics(fo)
        fo.close();
        


    def WriteDatacard(self,file,nameWS,rates,obsEvents,is2D):

        file.write("imax 1\n")
        file.write("jmax 7\n")
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS \n".format(nameWS))
        file.write("------------\n")
        
        file.write("bin a{0} \n".format(self.channel))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        
        file.write("bin ")
        for i in range(0,8): file.write("a{0} ".format(self.channel))
        file.write("\n")
        
        if not self.is2D:
            file.write("process ggH qqH WH ZH ttH bkg_qqzz bkg_ggzz bkg_zjets\n")
        else:
            file.write("process ggH qqH WH ZH ttH bkg2d_qqzz bkg2d_ggzz bkg2d_zjets\n")
            
        file.write("process -4 -3 -2 -1 0 1 2 3\n")
            
        file.write("rate ")
        for i in range(0,8): file.write("{0:.4f} ".format(rates[i]))
        file.write("\n")
        file.write("------------\n")






