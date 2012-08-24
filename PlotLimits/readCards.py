#! /usr/bin/env python
import os
import re
import math
from decimal import *
import ROOT
from ROOT import RooFit
from array import array

from optparse import OptionParser

############################################
#            Job steering                  #
############################################
parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-c', '--channel',
                  action='store', type='string', dest='channel',default='4mu',
                  help='string: channel to inspect; 4mu, 4e, 2e2mu')
parser.add_option('-i', '--inputDir',
                  action='store', type='string', dest='inputDir',default='/home/ntran/XZZAnalysis/HZZ4L_2012/v6_combo/ULtaskForce/inputs_UL_5fb',
                  help='string: dir location of datacards/workspaces')
parser.add_option('-o', '--outputDir',
                  action='store', type='string', dest='outputDir',default='figs',
                  help='string: where you want your plots to go')
parser.add_option('-w', '--workingArea',
                  action='store', type='string', dest='workingArea',default='/home/ntran/XZZAnalysis/HZZ4L_2012/v6_combo/HZZ4L_Combination',
                  help='obsolete!!! ... HZZ4L_Combination working area, needed for some libraries')

(options, args) = parser.parse_args()
############################################
############################################

ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/");
ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
from ROOT import setTDRStyle
ROOT.setTDRStyle(True)
ROOT.gStyle.SetPalette(1)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so");

getcontext().prec = 2


## ---------------------------------------------------------------
## card reader class
## ---------------------------------------------------------------
class cardreader:

    #-----------------        
    def __init__(self, dcname, wsfilename, oDir):
        self.theODir = oDir
            
        self.signals = ['ggH','qqH','WH','ZH','ttH']
        self.backgrounds = ['qqZZ','ggZZ','Zjets']

        fin = ROOT.TFile(wsfilename)
#        fin.ls()
        self.theCard = open(dcname)
        self.theWorkspace = fin.Get("w")
        self.theWorkspace.Print()
        print "mH: ",self.theWorkspace.var("MH").getVal()
        self.theMH = self.theWorkspace.var("MH").getVal()
        self.parseCard()
        self.allSysUp = [ ]*8
        self.allSysDn = [ ]*8        
        self.computeSystematicsSum()
        print self.allSysUp
        print self.allSysDn
        
        # define the channel    
        self.theChannel = ""
        if dcname.find( "4mu" ) >= 0: self.theChannel = "4mu"
        elif dcname.find( "4e" ) >= 0: self.theChannel = "4e"
        elif dcname.find( "2e2mu" ) >= 0: self.theChannel = "2e2mu"
        else: print "sorry this is just weird"

    #-----------------    
    def parseCard(self):
        startSystematics = False
        self.theSystematics = []
        ratelist = []
        for line in self.theCard:
            if line.find("rate") >= 0: ratelist = line.split()
            if line.find("lumi") >= 0: startSystematics = True
            if startSystematics: self.theSystematics.append( line )
        self.bkgYields = [ float(ratelist[6]),float(ratelist[7]),float(ratelist[8])]
#        print self.theSystematics
    #-----------------
    def computeSystematicsSum(self):
        curSysUp = []
        curSysDn = []        
        nSys = 0
        for line in self.theSystematics:
            rowSysUp = []
            rowSysDn = []            
            if line.find("lnN") >= 0:
#                print line
                nSys = nSys + 1
                linelist = line.split()
                print linelist
                for i in range(2,10):
                    if linelist[i] == "-": 
                        rowSysUp.append( 0 )
                        rowSysDn.append( 0 )                    
                    elif linelist[i].find('/') >= 0:
#                        print "hey special error"
                        asymSplit = linelist[i].split('/')
                        rowSysDn.append( float(asymSplit[0]) - 1. )
                        rowSysUp.append( float(asymSplit[1]) - 1. )                                          
                    else: 
                        rowSysUp.append( float(linelist[i]) - 1. )
                        rowSysDn.append( float(linelist[i]) - 1. )                                          
                curSysUp.append(rowSysUp)
                curSysDn.append(rowSysDn)                                          
#        print curSys
        sumSigSysSquaredUp = [0.]*8
        sumSigSysSquaredDn = [0.]*8                                          
#        print sumSigSysSquared
        for ii in range(nSys):
            for jj in range(8):
                sumSigSysSquaredUp[jj] = sumSigSysSquaredUp[jj] + curSysUp[ii][jj]**2
                sumSigSysSquaredDn[jj] = sumSigSysSquaredDn[jj] + curSysDn[ii][jj]**2                                          
        self.allSysUp = [ math.sqrt( sumSigSysSquaredUp[x] ) for x in xrange(len( sumSigSysSquaredUp )) ]    
        self.allSysDn = [ math.sqrt( sumSigSysSquaredDn[x] ) for x in xrange(len( sumSigSysSquaredDn )) ]                                              
    #-----------------    
    def getSignalYields(self):
        ggH_norm = self.theWorkspace.function("ggH_norm").getVal()
        qqH_norm = self.theWorkspace.function("qqH_norm").getVal()
        WH_norm = self.theWorkspace.function("WH_norm").getVal()
        ZH_norm = self.theWorkspace.function("ZH_norm").getVal()
        ttH_norm = self.theWorkspace.function("ttH_norm").getVal()        
        self.syields = [ggH_norm,qqH_norm,WH_norm,ZH_norm,ttH_norm]
        return self.syields
    #-----------------
    def getBackgroundYields(self):
        self.byields = [ self.bkgYields[0],
                   self.bkgYields[1],
                   self.bkgYields[2]
                   ]
        return self.byields
    #-----------------
    def getSignalErrorsUp(self):
        self.serrorsUp = [ self.allSysUp[0],
                    self.allSysUp[1],
                    self.allSysUp[2],
                    self.allSysUp[3],
                    self.allSysUp[4]
                   ]
        return self.serrorsUp
    #-----------------
    def getSignalErrorsDn(self):
        self.serrorsDn = [ self.allSysDn[0],
                self.allSysDn[1],
                self.allSysDn[2],
                self.allSysDn[3],
                self.allSysDn[4]
                ]
        return self.serrorsDn
    #-----------------
    def getBackgroundErrorsUp(self):
        self.berrorsUp = [ self.allSysUp[5],
                   self.allSysUp[6],
                   self.allSysUp[7]
                   ]
        return self.berrorsUp
    #-----------------
    def getBackgroundErrorsDn(self):
        self.berrorsDn = [ self.allSysDn[5],
                    self.allSysDn[6],
                    self.allSysDn[7]
                    ]
        return self.berrorsDn
    #-----------------
    def getMHWindow(self):
        lo = self.theWorkspace.var("CMS_zz4l_mass").getMin()
        hi = self.theWorkspace.var("CMS_zz4l_mass").getMax()
        self.bounds = [lo,hi]
        return self.bounds
    #-----------------
    def PrintLatex(self):
        _syields = self.getSignalYields()
        _serrors = self.getSignalErrors()
        _byields = self.getBackgroundYields()
        _berrors = self.getBackgroundErrors()            
        
        print "mH: ",self.theWorkspace.var("MH").getVal()
        for x in xrange(len(signals)):
            print self.signals[x]," = ",_syields[x]," +/- ",_syields[x]*_serrors[x]
        for x in xrange(len(backgrounds)):
            print self.backgrounds[x]," = ",_byields[x]," +/- ",_byields[x]*_berrors[x]

    def getYieldsAll(self):
        
        _syields = self.getSignalYields()
        _byields = self.getBackgroundYields()
        _serrorsUp = self.getSignalErrorsUp()                                          
        _berrorsUp = self.getBackgroundErrorsUp()            
        _serrorsDn = self.getSignalErrorsDn()                                          
        _berrorsDn = self.getBackgroundErrorsDn()            
        
        ## -------------------
        sumsig = 0.
        err2Up = 0.
        err2Dn = 0.        
        for i in xrange(len(self.signals)):
            sumsig = sumsig + _syields[i]
            err2Up = err2Up + _serrorsUp[i]**2
            err2Dn = err2Dn + _serrorsDn[i]**2

        ## -------------------
        CMS_zz4l_mass = self.theWorkspace.var("CMS_zz4l_mass")
        CMS_zz4l_mass.setRange("unblinded",121.5,130.5);
        CMS_zz4l_mass.setRange("fullrange",100.,800.);
        
        bkg_qqzz = self.theWorkspace.pdf("bkg_qqzz")
        bkg_ggzz = self.theWorkspace.pdf("bkg_ggzz")
        bkg_zjets = self.theWorkspace.pdf("bkg_zjets")

        curInt_qqzz = bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass) ).getVal();
        unblindedInt_qqzz = bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), "unblinded" ).getVal();
        fullInt_qqzz = bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), "fullrange" ).getVal();
          
        curInt_ggzz = bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass) ).getVal();
        unblindedInt_ggzz = bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), "unblinded" ).getVal();
        fullInt_ggzz = bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), "fullrange" ).getVal();

        curInt_zjets = bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass) ).getVal();
        unblindedInt_zjets = bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), "unblinded" ).getVal();
        fullInt_zjets = bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), "fullrange" ).getVal();

        yieldfull_qqzz = _byields[0]*fullInt_qqzz/curInt_qqzz
        yieldfull_ggzz = _byields[1]*fullInt_ggzz/curInt_ggzz
        yieldfull_zjets = _byields[2]*fullInt_zjets/curInt_zjets
        yieldUB_qqzz = _byields[0]*unblindedInt_qqzz/curInt_qqzz
        yieldUB_ggzz = _byields[1]*unblindedInt_ggzz/curInt_ggzz
        yieldUB_zjets = _byields[2]*unblindedInt_zjets/curInt_zjets

        print "yield signal: ",round(sumsig,2)," + ",round(math.sqrt(err2Up)*sumsig,2)," - ",round(math.sqrt(err2Dn)*sumsig,2)

        bkgyieldfull_zz = yieldfull_qqzz+yieldfull_ggzz
        bkgyieldfull_zj = yieldfull_zjets
        bkgerrorfullUp_zz = math.sqrt((_berrorsUp[0]*yieldfull_qqzz)**2 + (_berrorsUp[1]*yieldfull_ggzz)**2)
        bkgerrorfullUp_zj = (_berrorsUp[2]*yieldfull_zjets)
        bkgerrorfullDn_zz = math.sqrt((_berrorsDn[0]*yieldfull_qqzz)**2 + (_berrorsDn[1]*yieldfull_ggzz)**2)
        bkgerrorfullDn_zj = (_berrorsDn[2]*yieldfull_zjets)
        
        bkgyieldUB_zz = yieldUB_qqzz+yieldUB_ggzz
        bkgyieldUB_zj = yieldUB_zjets
        bkgerrorUBUp_zz = math.sqrt((_berrorsUp[0]*yieldUB_qqzz)**2 + (_berrorsUp[1]*yieldUB_ggzz)**2)                                          
        bkgerrorUBUp_zj = (_berrorsUp[2]*yieldUB_zjets)
        bkgerrorUBDn_zz = math.sqrt((_berrorsDn[0]*yieldUB_qqzz)**2 + (_berrorsDn[1]*yieldUB_ggzz)**2)                                          
        bkgerrorUBDn_zj = (_berrorsDn[2]*yieldUB_zjets)

        print "yield ggzz, 100-1000: ", round((yieldfull_ggzz),2), " + ", round((_berrorsUp[1]*yieldfull_ggzz),2), " - ", round((_berrorsDn[1]*yieldfull_ggzz),2)
        print "yield full range (zz), 100-1000: ", round((bkgyieldfull_zz),2), " + ", round((bkgerrorfullUp_zz),2), " - ", round((bkgerrorfullDn_zz),2)
        print "yield full range (zjets), 100-1000: ", round((bkgyieldfull_zj),2), " + ", round((bkgerrorfullUp_zj),2), " - ", round((bkgerrorfullDn_zj),2)

        print "yield unblinded (zz), 120-130: ", round((bkgyieldUB_zz),2), " + ", round((bkgerrorUBUp_zz),2), " - ", round((bkgerrorUBDn_zz),2)
        print "yield unblinded (zjets), 120-130: ", round((bkgyieldUB_zj),2), " + ", round((bkgerrorUBUp_zj),2), " - ", round((bkgerrorUBDn_zj),2)

        # 10 values in output
        # sigyield, sigerror
        # bkgyield, bkgerror (zz, full range)
        # bkgyield, bkgerror (zjets, full range)
        # bkgyield, bkgerror (zz, unblinded range)
        # bkgyield, bkgerror (zjets, unblinded range)
        output = [
                  sumsig,math.sqrt(err2Up)*sumsig,math.sqrt(err2Dn)*sumsig,
                  bkgyieldfull_zz,bkgerrorfullUp_zz,bkgerrorfullDn_zz,
                  bkgyieldfull_zj,bkgerrorfullUp_zj,bkgerrorfullDn_zj,
                  bkgyieldUB_zz,bkgerrorUBUp_zz,bkgerrorUBDn_zz,
                  bkgyieldUB_zj,bkgerrorUBUp_zj,bkgerrorUBDn_zj
                  ]
        return output
        
    def plotShapes( self ):
        _syields = self.getSignalYields()
        _byields = self.getBackgroundYields()
        
        ggH = self.theWorkspace.pdf("ggH")
        bkg_qqzz = self.theWorkspace.pdf("bkg_qqzz")
        bkg_ggzz = self.theWorkspace.pdf("bkg_ggzz")
        bkg_zjets = self.theWorkspace.pdf("bkg_zjets")
        CMS_zz4l_mass = self.theWorkspace.var("CMS_zz4l_mass")

        frameM4l = CMS_zz4l_mass.frame() ;
        super(ROOT.RooAbsPdf,ggH).plotOn( frameM4l, RooFit.LineColor( 1 ), RooFit.Normalization( _syields[0] ) )
#        bkg_qqzz.plotOn( frameM4l, RooFit.LineColor( 2 ) )
        super(ROOT.RooAbsPdf,bkg_qqzz).plotOn( frameM4l, RooFit.LineColor( 2 ), RooFit.Normalization( _byields[0] ) )
        super(ROOT.RooAbsPdf,bkg_ggzz).plotOn( frameM4l, RooFit.LineColor( 3 ), RooFit.Normalization( _byields[1] ) )
        super(ROOT.RooAbsPdf,bkg_zjets).plotOn( frameM4l, RooFit.LineColor( 6 ), RooFit.Normalization( _byields[2] ) )

        leg = ROOT.TLegend(0.25,0.7,0.45,0.9)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        h_ggH = ROOT.TH1F("h_ggH","h_ggH",100,0,1); h_ggH.SetLineColor( 1 );
        leg.AddEntry( h_ggH, "ggH sig", "l" ); 
        h_bkg_qqzz = ROOT.TH1F("h_bkg_qqzz","h_bkg_qqzz",100,0,1); h_bkg_qqzz.SetLineColor( 2 );
        leg.AddEntry( h_bkg_qqzz, "qqzz bkg", "l" ); 
        h_bkg_ggzz = ROOT.TH1F("h_bkg_ggzz","h_bkg_ggzz",100,0,1); h_bkg_ggzz.SetLineColor( 3 );
        leg.AddEntry( h_bkg_ggzz, "ggzz bkg", "l" ); 
        h_bkg_zjets = ROOT.TH1F("h_bkg_zjets","h_bkg_zjets",100,0,1); h_bkg_zjets.SetLineColor( 6 );
        leg.AddEntry( h_bkg_zjets, "zjets bkg", "l" ); 


        can = ROOT.TCanvas( "can", "can", 800, 800 )
        frameM4l.Draw()
        leg.Draw()
        can.Print( self.theODir + "/shapes_m" + str(self.theMH) + "_" + self.theChannel + ".eps", "eps" )

## ---------------------------------------------------------------


###########################################################
## M A I N   M A I N   M A I N   M A I N   M A I N
###########################################################

if __name__ == '__main__':
    
    print "workspace validation..."
#    ROOT.gSystem.AddIncludePath("-I"+options.workingArea+"/CreateDatacards/PDFs");
#    ROOT.gSystem.Load(options.workingArea+"/CreateDatacards/PDFs/HZZ4LRooPdfs_cc")

    channel = options.channel
    theDirectory = options.inputDir
    oDir=options.outputDir
    cardsAll = os.listdir(theDirectory+"/datacards/")
    workspacesAll = os.listdir(theDirectory+"/workspaces/")
    
    mH = []
    cards = {} #dictionary
    workspaces = {} #dictionary
    ## -------------------
    for i in xrange(len(cardsAll)): 
        if cardsAll[i].find(channel) >= 0 and cardsAll[i].startswith("hzz4l_"): 
            val = re.findall("[1-9][0-9][0-9].[0-9]",cardsAll[i])
            mHf = float( val[0] )
            mH.append( mHf )
            cards[ mHf ] = theDirectory+"/datacards/"+cardsAll[i]
            # find corresponding workspaces
            for j in xrange(len(workspacesAll)): 
                if workspacesAll[j].find(channel) >= 0 and workspacesAll[j].find(val[0]) >= 0: 
                    workspaces[ mHf ] = theDirectory+"/workspaces/"+workspacesAll[j]                        
    mH.sort()
    ## -------------------

#    print mH
#    print cards
#    print workspaces

###########################################################
    
    # masses on which we will plot the shapes (not all of them!)
    plotmass = [123,147,200,250,350,500]

    # graphs for arrays
    aMH = array( 'f', [ ] )
    # signal yields x5
    signals = ['ggH','qqH','WH','ZH','ttH']
    aSignalYields = []
    aSignalErrors = []
    for x in xrange(len(signals)):
        aSignalYields.append( array( 'f', [ ] ) )
        aSignalErrors.append( array( 'f', [ ] ) )
    # signal efficiency
    # background yields
    backgrounds = ['qqZZ','ggZZ','Zjets']
    aBackgroundYields = []
    aBackgroundErrors = []
    for x in xrange(len(signals)):
        aBackgroundYields.append( array( 'f', [ ] ) )
        aBackgroundErrors.append( array( 'f', [ ] ) )
    # mh windows
    loMH = array( 'f', [ ] )
    hiMH = array( 'f', [ ] )
    # xerrors
    ex = array( 'f', [ ] )
    # signal shapes
    ## ---------------------
    for ii in xrange(len(mH)):
        if mH[ii] > 100 and mH[ii] < 1000:
            print " --------- "
            print "mh: ",mH[ii],", cards: ",cards[mH[ii]],", ws: ",workspaces[mH[ii]]
            aMH.append( mH[ii] )
            reader1 = cardreader(cards[mH[ii]],workspaces[mH[ii]], oDir)    
            syields = reader1.getSignalYields()
            byields = reader1.getBackgroundYields()                                          
            # errors are only + errors for now
            serrors = reader1.getSignalErrorsUp()
            berrors = reader1.getBackgroundErrorsUp()            
            mhbounds = reader1.getMHWindow()
            allYields = reader1.getYieldsAll()
#            reader1.PrintLatex()
            
            for aa in range(len(plotmass)):
                if mH[ii] == plotmass[aa]: reader1.plotShapes()
#            print syields
#            print serrors
#            print byields
#            print berrors
#            print mhbounds

            for jj in xrange(len(syields)): 
                aSignalYields[jj].append( syields[jj] )
                aSignalErrors[jj].append( serrors[jj]*syields[jj] )
            for jj in xrange(len(byields)): 
                aBackgroundYields[jj].append( byields[jj] )
                aBackgroundErrors[jj].append( berrors[jj]*byields[jj] )
            print "mhbounds[0]: ",mhbounds[0],", mhbounds[1]:",mhbounds[1]
            loMH.append( mhbounds[0] )
            hiMH.append( mhbounds[1] )    
            ex.append(0.)

    grSignalYields = []
    grBackgroundYields = []            
    for x in xrange(len(signals)):
        grSignalYields.append( ROOT.TGraphErrors( len(aMH), aMH, aSignalYields[x], ex, aSignalErrors[x] ) )
        grSignalYields[x].SetName("grSignal_"+signals[x])
    for x in xrange(len(backgrounds)):
        grBackgroundYields.append( ROOT.TGraphErrors( len(aMH), aMH, aBackgroundYields[x], ex, aBackgroundErrors[x] ) )
        grBackgroundYields[x].SetName("grBackground_"+backgrounds[x])
    grLoMH = ROOT.TGraph( len(mH), aMH, loMH ); grLoMH.SetName("grLoMH")
    grHiMH = ROOT.TGraph( len(mH), aMH, hiMH ); grHiMH.SetName("grHiMH")


###########################################################
## --------- D R A W ---------

    # ----------
    leg = ROOT.TLegend(0.6,0.6,0.85,0.85)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    cSY = ROOT.TCanvas("cSY","cSY",1600,800)
    cSY.Divide(2,1)
    cSY.cd(1)   
    ROOT.gPad.SetGrid()
    hr = ROOT.gPad.DrawFrame(100, 0, 600, 10)
    hr.SetTitle("; mass; signal yield")
    for x in xrange(len(signals)):
        grSignalYields[x].SetLineColor( x+1 )
        grSignalYields[x].SetMarkerColor( x+1 )
        grSignalYields[x].Draw("pl")
        leg.AddEntry( grSignalYields[x], signals[x], "pl" )
    leg.Draw()
    cSY.cd(2)
    ROOT.gPad.SetGrid()
    hr2 = ROOT.gPad.DrawFrame(100, 1e-4, 600, 30)
    hr2.SetTitle("; mass; signal yield")
    for x in xrange(len(signals)):
        grSignalYields[x].SetLineColor( x+1 )
        grSignalYields[x].Draw("pl")
    ROOT.gPad.SetLogy()
    cSY.SaveAs(oDir+"/signalYields_"+channel+".eps")

    # ----------
    leg2 = ROOT.TLegend(0.6,0.6,0.85,0.85)
    leg2.SetFillColor(0)
    leg2.SetBorderSize(0)
    cBY = ROOT.TCanvas("cBY","cBY",1600,800)
    cBY.Divide(2,1)
    cBY.cd(1)   
    ROOT.gPad.SetGrid()
    hr = ROOT.gPad.DrawFrame(100, 0, 600, 50)
    hr.SetTitle("; mass; background yield")
    for x in xrange(len(backgrounds)):
        grBackgroundYields[x].SetLineColor( x+1 )
        grBackgroundYields[x].SetMarkerColor( x+1 )
        grBackgroundYields[x].Draw("pl")
        leg2.AddEntry( grBackgroundYields[x], backgrounds[x], "pl" )
    leg2.Draw()
    cBY.cd(2)
    ROOT.gPad.SetGrid()
    hr2 = ROOT.gPad.DrawFrame(100, 1e-4, 600, 100)
    hr2.SetTitle("; mass; background yield")
    for x in xrange(len(backgrounds)):
        grBackgroundYields[x].SetLineColor( x+1 )
        grBackgroundYields[x].Draw("pl")
    ROOT.gPad.SetLogy()
    cBY.SaveAs(oDir+"/backgroundYields_"+channel+".eps")

    # ----------
    cBounds = ROOT.TCanvas("cBounds","cBounds",800,800)
    ROOT.gPad.SetGrid()
    hr2 = ROOT.gPad.DrawFrame(100, 50, 600, 1000)
    grLoMH.SetLineColor(2); grLoMH.SetLineWidth(4);
    grHiMH.SetLineColor(4); grHiMH.SetLineWidth(4);
    grLoMH.Draw("l")
    grHiMH.Draw("l")    
    cBounds.SaveAs(oDir+"/mHBounds_"+channel+".eps")
###########################################################


    fout = ROOT.TFile(oDir+"/graphs_"+channel+".root","RECREATE")
    fout.cd()
    for x in xrange(len(signals)): grSignalYields[x].Write()
    for x in xrange(len(backgrounds)): grBackgroundYields[x].Write()
    grLoMH.Write()
    grHiMH.Write()
    fout.Close()
