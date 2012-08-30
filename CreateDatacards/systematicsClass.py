#! /usr/bin/env python
import os
import re
import math
from ROOT import *
import ROOT
from array import array


## ------------------------------------
##  systematics class
## ------------------------------------

class systematicsClass:

    def __init__(self,theMass,theForXSxBR,theisFSR,theInputs):

        self.ID_4mu = 1
        self.ID_4e = 2
        self.ID_2e2mu = 3

        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.mH = theMass
        self.isForXSxBR = theForXSxBR
        self.isFSR = theisFSR
        self.model = theInputs['model']


        self.muSelError = 0.0
        self.eSelError = 0.0
        self.muSelError2e2mu = 0.0
        self.eSelError2e2mu = 0.0
        self.muSelErrorZZ2e2mu = 0.0
        self.eSelErrorZZ2e2mu = 0.0

        self.qqVV_scaleSys = 0.0
        self.qqVV_pdfSys = 0.0
        self.ggVV_scaleSys = 0.0
        self.ggVV_pdfSys = 0.0

        self.lumiUncertainty = theInputs['lumiUnc']
        self.sel_muontrig = theInputs['muonTrigUnc']
        self.sel_muonfull = theInputs['muonFullUnc']

        self.sel_eletrig = theInputs['elecTrigUnc']
        self.sel_elefull = theInputs['elecFullUnc']

        self.zjetKappaLow = theInputs['zjetsKappaLow']
        self.zjetKappaHigh = theInputs['zjetsKappaHigh']

        self.theoryHighMass = 1
        
        if theInputs['muonTrigCutoff'] > 100 and theInputs['muonTrigUnc_HM'] > 0:
            if self.mH > theInputs['muonTrigCutoff']:
                self.sel_muontrig = theInputs['muonTrigUnc_HM']
                
        if theInputs['muonFullCutoff'] > 100 and theInputs['muonFullUnc_HM'] > 0:
            if self.mH > theInputs['muonFullCutoff']:
                self.sel_muonfull = theInputs['muonFullUnc_HM']
 
        if theInputs['elecTrigCutoff'] > 100 and theInputs['elecTrigUnc_HM'] > 0:
            if self.mH > theInputs['elecTrigCutoff']:
                self.sel_eletrig = theInputs['elecTrigUnc_HM']
                
        if theInputs['elecFullCutoff'] > 100 and theInputs['elecFullUnc_HM'] > 0:
            if self.mH > theInputs['elecFullCutoff']:
                self.sel_elefull = theInputs['elecFullUnc_HM']

        self.qqVV_scaleSys = 1. + 0.01*math.sqrt((self.mH - 20.)/13.)
        self.qqVV_pdfSys = 1. + 0.0035*math.sqrt(self.mH - 30.)
        self.ggVV_scaleSys = 1.04 + 0.10*math.sqrt((self.mH + 40.)/40.)
        self.ggVV_pdfSys = 1. + 0.0066*math.sqrt(self.mH - 10.)
            
          


    def setSystematics(self,theRateBkg_qqZZ,theRateBkg_ggZZ,theRateBkg_zjets ):

        self.rateBkg_qqZZ = theRateBkg_qqZZ
        self.rateBkg_ggZZ = theRateBkg_ggZZ
        self.rateBkg_zjets = theRateBkg_zjets

        if not self.model == "SM4" and not self.model == "SM" and not self.model == "FF" :
            
            print "In setSystematics in Systematics -------> Unknown model ",self.model
            print "Choices are SM, SM4, or FF"
            
        
        #ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        #ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")
        
        self.myCSW = HiggsCSandWidth()
        self.myCSWSM4 = HiggsCSandWidthSM4()

        if not self.isForXSxBR:
        
            self.CSpdfErrPlus_gg = self.myCSW.HiggsCSpdfErrPlus(1,self.mH,self.sqrts)
            self.CSpdfErrMinus_gg = self.myCSW.HiggsCSpdfErrMinus(1,self.mH,self.sqrts)
            self.CSpdfErrPlus_vbf = self.myCSW.HiggsCSpdfErrPlus(2,self.mH,self.sqrts)
            self.CSpdfErrMinus_vbf = self.myCSW.HiggsCSpdfErrMinus(2,self.mH,self.sqrts)
            self.CSpdfErrPlus_wh = self.myCSW.HiggsCSpdfErrPlus(3,self.mH,self.sqrts)
            self.CSpdfErrMinus_wh = self.myCSW.HiggsCSpdfErrMinus(3,self.mH,self.sqrts)
            self.CSpdfErrPlus_zh = self.myCSW.HiggsCSpdfErrPlus(4,self.mH,self.sqrts)
            self.CSpdfErrMinus_zh = self.myCSW.HiggsCSpdfErrMinus(4,self.mH,self.sqrts)
            self.CSpdfErrPlus_tth = self.myCSW.HiggsCSpdfErrPlus(5,self.mH,self.sqrts)
            self.CSpdfErrMinus_tth = self.myCSW.HiggsCSpdfErrMinus(5,self.mH,self.sqrts)
            
            self.CSscaleErrPlus_gg = self.myCSW.HiggsCSscaleErrPlus(1,self.mH,self.sqrts)
            self.CSscaleErrMinus_gg = self.myCSW.HiggsCSscaleErrMinus(1,self.mH,self.sqrts)
            self.CSscaleErrPlus_vbf = self.myCSW.HiggsCSscaleErrPlus(2,self.mH,self.sqrts)
            self.CSscaleErrMinus_vbf = self.myCSW.HiggsCSscaleErrMinus(2,self.mH,self.sqrts)
            self.CSscaleErrPlus_wh = self.myCSW.HiggsCSscaleErrPlus(3,self.mH,self.sqrts)
            self.CSscaleErrMinus_wh = self.myCSW.HiggsCSscaleErrMinus(3,self.mH,self.sqrts)
            self.CSscaleErrPlus_zh = self.myCSW.HiggsCSscaleErrPlus(4,self.mH,self.sqrts)
            self.CSscaleErrMinus_zh = self.myCSW.HiggsCSscaleErrMinus(4,self.mH,self.sqrts)
            self.CSscaleErrPlus_tth = self.myCSW.HiggsCSscaleErrPlus(5,self.mH,self.sqrts)
            self.CSscaleErrMinus_tth = self.myCSW.HiggsCSscaleErrMinus(5,self.mH,self.sqrts)
      

            if( self.mH >= 200): self.theoryHighMass = 1 + 1.5*(self.mH/1000)*(self.mH/1000)*(self.mH/1000)
      
            self.BRErr_Hff = self.myCSWSM4.HiggsBRErr_Hff(11,self.mH,self.sqrts)
            self.BRErr_HVV = self.myCSWSM4.HiggsBRErr_HVV(11,self.mH,self.sqrts)
            self.BRErr_Hgg = self.myCSWSM4.HiggsBRErr_Hgluglu(11,self.mH,self.sqrts)
      

        
        self.muSelError = 1 + math.sqrt( self.sel_muonfull*self.sel_muonfull + self.sel_muontrig*self.sel_muontrig )
        self.eSelError = 1 + math.sqrt( self.sel_elefull*self.sel_elefull + self.sel_eletrig*self.sel_eletrig )
        self.muSelError2e2mu = 1 + math.sqrt( self.sel_muonfull*self.sel_muonfull/4 + self.sel_muontrig*self.sel_muontrig/4 )
        self.eSelError2e2mu = 1 + math.sqrt( self.sel_elefull*self.sel_elefull/4 + self.sel_eletrig*self.sel_eletrig/4 )
        
        

    def WriteSystematics(self,theFile,theInputs):

        if theInputs['useLumiUnc']:
            if(self.sqrts == 7):
                theFile.write("lumi_7TeV lnN ")
            elif (self.sqrts == 8):
                theFile.write("lumi_8TeV lnN ")
            else:
                raise RuntimeError, "Unknown sqrts in systematics!"
            for i in range(0,7):
                theFile.write("{0} ".format(self.lumiUncertainty))
            theFile.write(" -\n")
        
		
        if not self.isForXSxBR:
            
            if not self.model == "FF" and theInputs['usePdf_gg']:
                theFile.write("pdf_gg lnN {0:.4f} - - - {1:.4f} - {2:.4f} - \n".format( 1 + (self.CSpdfErrPlus_gg-self.CSpdfErrMinus_gg)/2., 1. + (self.CSpdfErrPlus_tth-self.CSpdfErrMinus_tth)/2, self.ggVV_pdfSys))

            if theInputs['usePdf_qqbar']:
                theFile.write("pdf_qqbar lnN - {0:.4f} {1:.4f} {2:.4f} - {3:.4f} - - \n".format( 1. + (self.CSpdfErrPlus_vbf-self.CSpdfErrMinus_vbf)/2.,1. + (self.CSpdfErrPlus_wh-self.CSpdfErrMinus_wh)/2.,1. + (self.CSpdfErrPlus_zh-self.CSpdfErrMinus_zh)/2.,self.qqVV_pdfSys))

	    if theInputs['usePdf_hzz4l_accept']:
                theFile.write("pdf_hzz4l_accept lnN 1.02 1.02 1.02 1.02 1.02 - - - \n")
	    
            if not self.model == "FF" and theInputs['useQCDscale_ggH']:
                theFile.write("QCDscale_ggH lnN {0:.4f} - - - - - - - \n".format(1. + (self.CSscaleErrPlus_gg-self.CSscaleErrMinus_gg)/2.))

	    if theInputs['useQCDscale_qqH']:
                theFile.write("QCDscale_qqH lnN - {0:.4f} - - - - - - \n".format(1. + (self.CSscaleErrPlus_vbf-self.CSscaleErrMinus_vbf)/2.))

	    if theInputs['useQCDscale_VH']:
                theFile.write("QCDscale_VH lnN - - {0:.4f} {1:.4f} - - - - \n".format(1. + (self.CSscaleErrPlus_wh-self.CSscaleErrMinus_wh)/2.,1. + (self.CSscaleErrPlus_zh-self.CSscaleErrMinus_zh)/2.))
	    
            if not self.model == "FF" and theInputs['useQCDscale_ttH']:
                theFile.write("QCDscale_ttH lnN - - - - {0:.4f} - - - \n".format(1. + (self.CSscaleErrPlus_tth-self.CSscaleErrMinus_tth)/2.))

	    if theInputs['useTheoryUncXS_HighMH']:
                theFile.write("theoryUncXS_HighMH lnN {0:.3f} {0:.3f} {0:.3f} {0:.3f} {0:.3f} - - - \n".format(self.theoryHighMass))

	elif self.isForXSxBR:
            if theInputs['usePdf_gg']:
                theFile.write("pdf_gg lnN - - - - - - {0:.4f} - \n".format(self.ggVV_pdfSys))
            if theInputs['usePdf_qqbar']:
                theFile.write("pdf_qqbar lnN - - - - - {0:.4f} - - \n".format(self.qqVV_pdfSys))
		
        if theInputs['useQCDscale_ggVV']:
            theFile.write("QCDscale_ggVV lnN - - - - - - {0:.4f} - \n".format(self.ggVV_scaleSys))
        if theInputs['useQCDscale_VV']:
            theFile.write("QCDscale_VV lnN - - - - - {0:.4f} - - \n".format(self.qqVV_scaleSys))
	
	## Higgs BR
        if(self.model == "SM" or self.model == "FF") and theInputs['useBRhiggs_hzz4l']:
            theFile.write("BRhiggs_hzz4l lnN 1.02 1.02 1.02 1.02 1.02 - - -\n")
	  
	elif(self.model == "SM4"):
	  
	    theFile.write("gamma_Hff lnN ")
            for i in range(0,5): theFile.write("{0:.4f} ".format(self.BRErr_Hff))
	    theFile.write("- - - \n")

	    theFile.write("gamma_HVV lnN ")
            for i in range(0,5): theFile.write("{0:.4f} ".format(self.BRErr_HVV))
	    theFile.write("- - - \n")
            
	    theFile.write("gamma_Hgluglu lnN ")
            for i in range(0,5): theFile.write("{0:.4f} ".format(self.BRErr_Hgg))
	    theFile.write("- - - \n")
	  

	
	##  ----------- SELECTION EFFICIENCIES ----------
            
        if theInputs['useCMS_eff']:
            if ( self.channel == self.ID_4mu ):
	  
                theFile.write("CMS_eff_m lnN ")
                for i in range(0,7): theFile.write("{0:.3f} ".format(self.muSelError))
                theFile.write("-\n")

            if ( self.channel == self.ID_4e ):
	  
                theFile.write("CMS_eff_e lnN ")
                for i in range(0,7): theFile.write("{0:.3f} ".format(self.eSelError))
                theFile.write("-\n")

            if( self.channel == self.ID_2e2mu ):
	  
                theFile.write("CMS_eff_m lnN ")
                for i in range(0,7): theFile.write("{0:.3f} ".format(self.muSelError2e2mu))
                theFile.write("-\n")
        
                theFile.write("CMS_eff_e lnN ")
                for i in range(0,7): theFile.write("{0:.3f} ".format(self.eSelError2e2mu))
                theFile.write("-\n")

	if (self.channel == self.ID_4mu):

            if theInputs['useCMS_hzz4l_Zjets']:
                theFile.write("CMS_hzz4mu_Zjets lnN - - - - - - - {0}/{1} \n".format(self.zjetKappaLow,self.zjetKappaHigh))
            if theInputs['useCMS_zz4l_bkgMELA']:
                theFile.write("CMS_zz4l_bkgMELA param 0  1  [-3,3]\n")
	    if theInputs['useCMS_zz4l_sigMELA']:
                theFile.write("CMS_zz4l_sigMELA param 0  1  [-3,3]\n")
	  
	if (self.channel == self.ID_4e):

            if theInputs['useCMS_hzz4l_Zjets']:
                theFile.write("CMS_hzz4e_Zjets lnN - - - - - - - {0}/{1} \n".format(self.zjetKappaLow,self.zjetKappaHigh))
            if theInputs['useCMS_zz4l_bkgMELA']:
                theFile.write("CMS_zz4l_bkgMELA param 0  1  [-3,3]\n")
            if theInputs['useCMS_zz4l_sigMELA']:
                theFile.write("CMS_zz4l_sigMELA param 0  1  [-3,3]\n")
	  
	if (self.channel == self.ID_2e2mu):

            if theInputs['useCMS_hzz4l_Zjets']:
                theFile.write("CMS_hzz2e2mu_Zjets lnN - - - - - - - {0}/{1} \n".format(self.zjetKappaLow,self.zjetKappaHigh))
            if theInputs['useCMS_zz4l_bkgMELA']:
                theFile.write("CMS_zz4l_bkgMELA param 0  1  [-3,3]\n")
            if theInputs['useCMS_zz4l_sigMELA']:
                theFile.write("CMS_zz4l_sigMELA param 0  1  [-3,3]\n")

    

    def WriteShapeSystematics(self,theFile,theInputs):
  
        meanCB_e_errPerCent = theInputs['CMS_zz4l_mean_e_sig']
        sigmaCB_e_errPerCent = theInputs['CMS_zz4l_sigma_e_sig']
        N_CB_errPerCent = theInputs['CMS_zz4l_n_sig']
        meanCB_m_errPerCent = theInputs['CMS_zz4l_mean_m_sig']
        sigmaCB_m_errPerCent = theInputs['CMS_zz4l_sigma_m_sig']
        
        if( self.channel == self.ID_4mu):

            if theInputs['useCMS_zz4l_mean']:
                theFile.write("CMS_zz4l_mean_m_sig param 0.0 {0} \n".format(meanCB_m_errPerCent))
            if theInputs['useCMS_zz4l_sigma']:
                theFile.write("CMS_zz4l_sigma_m_sig param 0.0 {0} \n".format(sigmaCB_m_errPerCent))
            if theInputs['useCMS_zz4l_n']:
                theFile.write("CMS_zz4l_n_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,N_CB_errPerCent))

        if( self.channel == self.ID_4e):

            if theInputs['useCMS_zz4l_mean']:
                theFile.write("CMS_zz4l_mean_e_sig param 0.0 {0} \n".format(meanCB_e_errPerCent))
            if theInputs['useCMS_zz4l_sigma']:
                theFile.write("CMS_zz4l_sigma_e_sig param 0.0 {0} \n".format(sigmaCB_e_errPerCent))
            if theInputs['useCMS_zz4l_n']:
                theFile.write("CMS_zz4l_n_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,N_CB_errPerCent))
            
        if( self.channel == self.ID_2e2mu):

            if theInputs['useCMS_zz4l_mean']:
                theFile.write("CMS_zz4l_mean_m_sig param 0.0 {0} \n".format(meanCB_m_errPerCent))
                theFile.write("CMS_zz4l_mean_e_sig param 0.0 {0} \n".format(meanCB_e_errPerCent))
            if theInputs['useCMS_zz4l_sigma']:
                theFile.write("CMS_zz4l_sigma_m_sig param 0.0 {0} \n".format(sigmaCB_m_errPerCent))
                theFile.write("CMS_zz4l_sigma_e_sig param 0.0 {0} \n".format(sigmaCB_e_errPerCent))
            if theInputs['useCMS_zz4l_n']:
                theFile.write("CMS_zz4l_n_sig_{0}_{1:.0f} param 0.0 {2} \n".format(self.channel,self.sqrts,N_CB_errPerCent))
                                                





