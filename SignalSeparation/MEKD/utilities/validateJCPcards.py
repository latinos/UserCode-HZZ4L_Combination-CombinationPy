#! /usr/bin/env python
import sys, os, pwd, commands
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array

m_cmssw_base=commands.getoutput("echo $CMSSW_BASE")
m_curdir=commands.getoutput("pwd")


ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
ROOT.gSystem.AddIncludePath("-Iutilities/")
ROOT.gROOT.ProcessLine(".L utilities/tdrstyle.cc")
ROOT.gSystem.Load("libRooFit")
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")


def parseOptions():
    
    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-i', '--input', dest='inputDir', type='string', default="",    help='directory containing all cards/ws')
    parser.add_option('-o', '--output', dest='outputDir', type='string', default="validation", help='output directory')
    parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
    
    parser.add_option('-x', '--xkd', dest='xKD', type='string', default="CMS_zz4l_smd",    help='name of x axis kd')
    parser.add_option('-y', '--ykd', dest='yKD', type='string', default="CMS_zz4l_pseudoKD", help='name of y axis kd')
    parser.add_option('-q', '--quick',action='store_true', dest='quickRun', default=False, help='Skip 1D slices to run quick')
    

    global opt, args
    (opt, args) = parser.parse_args()
        

    opt.xKDbinning = 'smd_binning'
    opt.yKDbinning = 'psD_binning'

    if opt.inputDir == '':
        raise RuntimeError, 'Please choose and input directory!'

            
def makeDirectory(subDirName):
    if (not os.path.exists(subDirName)):
        cmd = 'mkdir -p '+subDirName
        status, output = commands.getstatusoutput(cmd)
        if status !=0:
            print 'Error in creating submission dir '+subDirName+'. Exiting...'
            sys.exit()
    else:
        print 'Directory '+subDirName+' already exists. Exiting...'
        sys.exit()
                                                                        
def processCmd(cmd):
    #    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()



def printBinning(hist):

    nBinsX = hist.GetXaxis().GetNbins()
    nBinsY = hist.GetYaxis().GetNbins()


    print ">>>> "+str(hist.GetName())+" binning[x,y]: ["+str(nBinsX)+","+str(nBinsY)+"]"

    binString = '  '
    spacer = '   --'
    for x in xrange(1,nBinsX+1):
        if x > 9:
            binString += str('    {0}'.format(x))
            spacer += '------'
        else:
            binString += str('     {0}'.format(x))
            spacer += '------'
        

    binCounterY = nBinsY

    print ">>>> "+str(hist.GetName())+" bin areas*100 (for formatting ease):"


    for y in range(nBinsY,0,-1):
        widthString = '{0} | '.format(y)
        if y < 10:
            widthString = '{0}  | '.format(y)
        for x in range(1,nBinsX+1):
            Xsize = hist.GetXaxis().GetBinWidth(x)
            Ysize = hist.GetYaxis().GetBinWidth(y)
            area = Xsize*Ysize*100
            widthString += str("{0:.3f}".format(area))+" "
        print widthString

    print spacer
    print binString

    print "Integral: ",hist.Integral()
    print "Integral(\"width\"): ",hist.Integral("width")
        

def validate():
    global opt, args

    parseOptions()
        
    os.chdir(opt.inputDir)

    makeDirectory(opt.outputDir)

    masses = [126]
    channels = ['4e','4mu','2e2mu']
    chanNum = [2,1,3]
    energies = [7,8]
    

    for mass in masses:
        print ">>>> MASS: ",mass

        for sqrts in energies:
            print ">>>> SqrtS: ", sqrts,"TeV"
            i=0
            for chan in channels:
                print ">>>> Channel: ",chan

                filename = "hzz4l_{0}S_{1}TeV.input.root".format(chan,sqrts)
                print filename
                wsFile = ROOT.TFile(filename,"READ")

                ws = wsFile.Get("w")
                ws.var("CMS_zz4l_mass").setVal(mass)

                print ">>>> Normalizations: "

                filename = "hzz4l_{0}S_{1}TeV_ALT.txt".format(chan,sqrts)
                if not os.path.exists(filename):
                    raise RuntimeError, "File {0} does not exist!!!".format(filename)
                for line in open(filename,'r'):
                    f = line.split()
                    if len(f) < 1: continue
                    if f[0].startswith("#"): continue
                    if len(f) > 1 and f[1].startswith("ggH"):
                        print line
                    if f[0].startswith("rate"):
                        print line


                print "ggH: ",ws.function("ggH_norm").getVal()
                print "ggH_ALT: ",ws.function("ggH_ALT_norm").getVal()
                print "\n\n"
                
                ##### PDFs #####
                print ">>>> PDFs: "
                
                xKDvar = ws.var(opt.xKD)
                yKDvar = ws.var(opt.yKD)
                
                xkd_binning = ws.var(opt.xKD).getBinning(opt.xKDbinning)
                ykd_binning = ws.var(opt.yKD).getBinning(opt.yKDbinning)

                xkd_binning.Print()
                ykd_binning.Print()
                
                ggzz_2d = ws.pdf("bkg2d_ggzz").createHistogram("bkg2d_ggzz_pdf",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))
                qqzz_2d = ws.pdf("bkg2d_qqzz").createHistogram("bkg2d_qqzz_pdf",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))
                zjets_2d = ws.pdf("bkg2d_zjets").createHistogram("bkg2d_zjets_pdf",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))
                ggH_2d = ws.pdf("ggH").createHistogram("ggH_pdf",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))
                ggH_ALT_2d = ws.pdf("ggH_ALT").createHistogram("ggH_ALT_pdf",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))

                canv = ROOT.TCanvas("c","c",700,700)
                canv.cd()
                
                ggzz_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/ggzz_pdf_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))
                
                qqzz_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/qqzz_pdf_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))
                
                zjets_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/zjets_pdf_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))
                
                ggH_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/ggH_pdf_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))
                
                ggH_ALT_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/ggH_ALT_pdf_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))

                printBinning(ggH_2d)
                printBinning(ggH_ALT_2d)
                printBinning(qqzz_2d)
                printBinning(ggzz_2d)
                printBinning(zjets_2d)
 

                ##### Templates #####
                print ">>>> Templates:"
                Tggzz_2d = ws.pdf("bkgTemplatePdf_ggzz_{0}_{1}".format(chanNum[i],sqrts)).createHistogram("bkg2d_ggzz_template",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))
                Tqqzz_2d = ws.pdf("bkgTemplatePdf_qqzz_{0}_{1}".format(chanNum[i],sqrts)).createHistogram("bkg2d_qqzz_template",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))
                Tzjets_2d = ws.pdf("bkgTemplatePdf_zjets_{0}_{1}".format(chanNum[i],sqrts)).createHistogram("bkg2d_zjets_template",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))
                TggH_2d = ws.pdf("sigTemplatePdf_ggH_{0}_{1}".format(chanNum[i],sqrts)).createHistogram("ggH_template",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))
                TggH_ALT_2d = ws.pdf("sigTemplatePdf_ggH_ALT_{0}_{1}".format(chanNum[i],sqrts)).createHistogram("ggH_ALT_template",xKDvar,ROOT.RooFit.Binning(xkd_binning),ROOT.RooFit.YVar(yKDvar,ROOT.RooFit.Binning(ykd_binning)))


                canv = ROOT.TCanvas("c","c",700,700)
                canv.cd()
                
                Tggzz_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/ggzz_template_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))
                
                Tqqzz_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/qqzz_template_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))
                
                Tzjets_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/zjets_template_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))
                
                TggH_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/ggH_template_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))
                
                TggH_ALT_2d.Draw("PROF COLZ")
                canv.SaveAs("{0}/ggH_ALT_template_{1}_{2}TeV.png".format(opt.outputDir,chan,sqrts))

                printBinning(TggH_2d)
                printBinning(TggH_ALT_2d)
                printBinning(Tqqzz_2d)
                printBinning(Tggzz_2d)
                printBinning(Tzjets_2d)
 

                #### Dataset ####
                print ">>>> Dataset: "
                
                data = ws.data("data_obs")
                              
                for x in xrange(data.numEntries()):
                    xKDval = data.get(x).find(opt.xKD).getVal()
                    yKDval = data.get(x).find(opt.yKD).getVal()
                    print "{0}: {1}   {2}: {3}".format(opt.xKD,xKDval,opt.yKD,yKDval)
                    
                
                if not opt.quickRun:
                #### Plot 1D slices ####
                    
                    xKDvar.setBins(ggH_2d.GetXaxis().GetNbins())
                    yKDvar.setBins(ggH_2d.GetYaxis().GetNbins())
                    
                    xkd_frame = xKDvar.frame()
                    ykd_frame = yKDvar.frame()

                    dummyHist_ggzz = ROOT.TH1F()
                    dummyHist_qqzz = ROOT.TH1F()
                    dummyHist_zjets = ROOT.TH1F()
                    dummyHist_ggH = ROOT.TH1F()
                    dummyHist_ggH_ALT = ROOT.TH1F()

                    dummyHist_ggzz.SetLineColor(2)
                    dummyHist_qqzz.SetLineColor(3)
                    dummyHist_zjets.SetLineColor(4)
                    dummyHist_ggH.SetLineColor(1)
                    dummyHist_ggH_ALT.SetLineColor(6)

                    dummyHist_ggzz.SetLineStyle(2)

                    leg = ROOT.TLegend(0.75,0.75,1,1)
                    leg.SetBorderSize(0)
                    leg.SetTextFont(42)
                    leg.AddEntry(dummyHist_ggzz,"ggZZ","L")
                    leg.AddEntry(dummyHist_qqzz,"qqZZ","L")
                    leg.AddEntry(dummyHist_zjets,"Zjets","L")
                    leg.AddEntry(dummyHist_ggH,"ggH","L")
                    leg.AddEntry(dummyHist_ggH_ALT,"ggH_ALT","L")

                    
                    data.plotOn(xkd_frame,ROOT.RooFit.DataError(ROOT.RooAbsData.None))
                    ws.pdf("bkg2d_ggzz").plotOn(xkd_frame,ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
                    ws.pdf("bkg2d_qqzz").plotOn(xkd_frame,ROOT.RooFit.LineStyle(1),ROOT.RooFit.LineColor(3))
                    ws.pdf("bkg2d_zjets").plotOn(xkd_frame,ROOT.RooFit.LineStyle(1),ROOT.RooFit.LineColor(4))
                    ws.pdf("ggH").plotOn(xkd_frame,ROOT.RooFit.LineStyle(1),ROOT.RooFit.LineColor(1))
                    ws.pdf("ggH_ALT").plotOn(xkd_frame,ROOT.RooFit.LineStyle(1),ROOT.RooFit.LineColor(6))
                    xkd_frame.Draw()
                    leg.Draw()
                    canv.SaveAs("{0}/xKDslice_{1}_{2}_{3}TeV.png".format(opt.outputDir,mass,chan,sqrts))
                    
                    data.plotOn(ykd_frame,ROOT.RooFit.DataError(ROOT.RooAbsData.None))
                    ws.pdf("bkg2d_ggzz").plotOn(ykd_frame,ROOT.RooFit.LineStyle(2),ROOT.RooFit.LineColor(2))
                    ws.pdf("bkg2d_qqzz").plotOn(ykd_frame,ROOT.RooFit.LineStyle(1),ROOT.RooFit.LineColor(3))
                    ws.pdf("bkg2d_zjets").plotOn(ykd_frame,ROOT.RooFit.LineStyle(1),ROOT.RooFit.LineColor(4))
                    ws.pdf("ggH").plotOn(ykd_frame,ROOT.RooFit.LineStyle(1),ROOT.RooFit.LineColor(1))
                    ws.pdf("ggH_ALT").plotOn(ykd_frame,ROOT.RooFit.LineStyle(1),ROOT.RooFit.LineColor(6))
                    ykd_frame.Draw()
                    leg.Draw()
                    canv.SaveAs("{0}/yKDslice_{1}_{2}_{3}TeV.png".format(opt.outputDir,mass,chan,sqrts))
                    
                i+=1
                










if __name__ == "__main__":
    validate()
        
