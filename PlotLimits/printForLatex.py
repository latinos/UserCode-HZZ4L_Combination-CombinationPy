#! /usr/bin/env python
import os
import re
import math
import ROOT
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
parser.add_option('-t', '--textFile',
                  action='store', type='string', dest='textFile',default='test_1.txt',
                  help='spits out a text file with the table information as well')

(options, args) = parser.parse_args()
############################################
############################################

ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/");
ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
from ROOT import setTDRStyle
ROOT.setTDRStyle(True)
ROOT.gStyle.SetPalette(1)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so");

from readCards import cardreader

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
    cards_4e = {} #dictionary
    workspaces_4e = {} #dictionary
    channel = "4e"
    ## -------------------
    for i in xrange(len(cardsAll)): 
        if cardsAll[i].find(channel) >= 0 and cardsAll[i].startswith("hzz4l_"): 
            val = re.findall("[1-9][0-9][0-9].[0-9]",cardsAll[i])
            mHf = float( val[0] )
            mH.append( mHf )
            cards_4e[ mHf ] = theDirectory+"/datacards/"+cardsAll[i]
            # find corresponding workspaces
            for j in xrange(len(workspacesAll)): 
                if workspacesAll[j].find(channel) >= 0 and workspacesAll[j].find(val[0]) >= 0: 
                    workspaces_4e[ mHf ] = theDirectory+"/workspaces/"+workspacesAll[j]                        
    mH.sort()
    ## -------------------
    cards_4mu = {} #dictionary
    workspaces_4mu = {} #dictionary
    channel = "4mu"
    ## -------------------
    for i in xrange(len(cardsAll)): 
        if cardsAll[i].find(channel) >= 0 and cardsAll[i].startswith("hzz4l_"): 
            val = re.findall("[1-9][0-9][0-9].[0-9]",cardsAll[i])
            mHf = float( val[0] )
            cards_4mu[ mHf ] = theDirectory+"/datacards/"+cardsAll[i]
            # find corresponding workspaces
            for j in xrange(len(workspacesAll)): 
                if workspacesAll[j].find(channel) >= 0 and workspacesAll[j].find(val[0]) >= 0: 
                    workspaces_4mu[ mHf ] = theDirectory+"/workspaces/"+workspacesAll[j]                        
    ## -------------------  
    cards_2e2mu = {} #dictionary
    workspaces_2e2mu = {} #dictionary
    channel = "2e2mu"
    ## -------------------
    for i in xrange(len(cardsAll)): 
        if cardsAll[i].find(channel) >= 0 and cardsAll[i].startswith("hzz4l_"): 
            val = re.findall("[1-9][0-9][0-9].[0-9]",cardsAll[i])
            mHf = float( val[0] )
            cards_2e2mu[ mHf ] = theDirectory+"/datacards/"+cardsAll[i]
            # find corresponding workspaces
            for j in xrange(len(workspacesAll)): 
                if workspacesAll[j].find(channel) >= 0 and workspacesAll[j].find(val[0]) >= 0: 
                    workspaces_2e2mu[ mHf ] = theDirectory+"/workspaces/"+workspacesAll[j]                        

#    print mH
#    print cards
#    print workspaces

###########################################################
    
    #massToPrint = [125.,126.,130.,200.,350.,500.]
    massToPrint = [125.]
    nMasses = len(massToPrint)

    # signal shapes
    yields_4e = []
    yields_4mu = []
    yields_2e2mu = []
    ## ---------------------
    for ii in xrange(len(mH)):
        for jj in xrange(len(massToPrint)):
            if mH[ii] == massToPrint[jj]:
                print " --------- "
                print "mh: ",mH[ii],", cards: ",cards_4mu[mH[ii]],", ws: ",workspaces_4mu[mH[ii]]
                reader1 = cardreader(cards_4mu[mH[ii]],workspaces_4mu[mH[ii]], oDir)    
                allYields_4mu = reader1.getYieldsAll()
                reader2 = cardreader(cards_4e[mH[ii]],workspaces_4e[mH[ii]], oDir)    
                allYields_4e = reader2.getYieldsAll()
                reader3 = cardreader(cards_2e2mu[mH[ii]],workspaces_2e2mu[mH[ii]], oDir)    
                allYields_2e2mu = reader3.getYieldsAll()

#                print allYields_4e
                yields_4e.append(allYields_4e)
                yields_4mu.append(allYields_4mu)
                yields_2e2mu.append(allYields_2e2mu)

    #### ------ printing
    isy_cv = 0; isy_eu = 1; isy_ed = 2;
    iby_qq = 3; ibeu_qq = 4; ibed_qq = 5;
    iby_zj = 6; ibeu_zj = 7; ibed_zj = 8;
    iby_qqub = 9; ibeu_qqub = 10; ibed_qqub = 11;
    iby_zjub = 12; ibeu_zjub = 13; ibed_zjub = 14;

    ## from is [central, err Up, err Dn]
    totbkg_4e = yields_4e[nMasses-1][iby_qq] + yields_4e[nMasses-1][iby_zj]
    totbkgerrUp_4e =  math.sqrt( yields_4e[nMasses-1][ibeu_qq]**2 + yields_4e[nMasses-1][ibeu_zj]**2 )
    totbkgerrDn_4e =  math.sqrt( yields_4e[nMasses-1][ibed_qq]**2 + yields_4e[nMasses-1][ibed_zj]**2 )
    totbkg_4mu = yields_4mu[nMasses-1][iby_qq] + yields_4mu[nMasses-1][iby_zj]
    totbkgerrUp_4mu =  math.sqrt( yields_4mu[nMasses-1][ibeu_qq]**2 + yields_4mu[nMasses-1][ibeu_zj]**2 )
    totbkgerrDn_4mu =  math.sqrt( yields_4mu[nMasses-1][ibed_qq]**2 + yields_4mu[nMasses-1][ibed_zj]**2 )
    totbkg_2e2mu = yields_2e2mu[nMasses-1][iby_qq] + yields_2e2mu[nMasses-1][iby_zj]
    totbkgerrUp_2e2mu =  math.sqrt( yields_2e2mu[nMasses-1][ibeu_qq]**2 + yields_2e2mu[nMasses-1][ibeu_zj]**2 )
    totbkgerrDn_2e2mu =  math.sqrt( yields_2e2mu[nMasses-1][ibed_qq]**2 + yields_2e2mu[nMasses-1][ibed_zj]**2 )

    UBbkg_4e = yields_4e[nMasses-1][iby_qqub] + yields_4e[nMasses-1][iby_zjub]
    UBbkgerrUp_4e =  math.sqrt( yields_4e[nMasses-1][ibeu_qqub]**2 + yields_4e[nMasses-1][ibeu_zjub]**2 )
    UBbkgerrDn_4e =  math.sqrt( yields_4e[nMasses-1][ibed_qqub]**2 + yields_4e[nMasses-1][ibed_zjub]**2 )
    UBbkg_4mu = yields_4mu[nMasses-1][iby_qqub] + yields_4mu[nMasses-1][iby_zjub]
    UBbkgerrUp_4mu =  math.sqrt( yields_4mu[nMasses-1][ibeu_qqub]**2 + yields_4mu[nMasses-1][ibeu_zjub]**2 )
    UBbkgerrDn_4mu =  math.sqrt( yields_4mu[nMasses-1][ibed_qqub]**2 + yields_4mu[nMasses-1][ibed_zjub]**2 )
    UBbkg_2e2mu = yields_2e2mu[nMasses-1][iby_qqub] + yields_2e2mu[nMasses-1][iby_zjub]
    UBbkgerrUp_2e2mu =  math.sqrt( yields_2e2mu[nMasses-1][ibeu_qqub]**2 + yields_2e2mu[nMasses-1][ibeu_zjub]**2 )
    UBbkgerrDn_2e2mu =  math.sqrt( yields_2e2mu[nMasses-1][ibed_qqub]**2 + yields_2e2mu[nMasses-1][ibed_zjub]**2 )
        
    print "\n"
    print "\n"
    print "\n"
    ## print individual backgrounds
    print "ZZ background & ",round(yields_4e[nMasses-1][iby_qq],2)," \pm ",round(yields_4e[nMasses-1][ibeu_qq],2)," & ",round(yields_4mu[nMasses-1][iby_qq],2)," \pm ",round(yields_4mu[nMasses-1][ibeu_qq],2)," & ",round(yields_2e2mu[nMasses-1][iby_qq],2)," \pm ",round(yields_2e2mu[nMasses-1][ibeu_qq],2)," \\\\"
    
    if yields_4e[nMasses-1][ibeu_qq] == yields_4e[nMasses-1][ibed_qq] and yields_4mu[nMasses-1][ibeu_qq] == yields_4e[nMasses-1][ibed_qq] and yields_2e2mu[nMasses-1][ibeu_qq] == yields_4e[nMasses-1][ibed_qq]:
        print "Z + X & ",round(yields_4e[nMasses-1][iby_zj],2)," \pm ",round(yields_4e[nMasses-1][ibeu_zj],2)," & ",round(yields_4mu[nMasses-1][iby_zj],2)," \pm ",round(yields_4mu[nMasses-1][ibeu_zj],2)," & ",round(yields_2e2mu[nMasses-1][iby_zj],2)," \pm ",round(yields_2e2mu[nMasses-1][ibeu_zj],2)," \\\\"
    else:
        print "Z + X & ",round(yields_4e[nMasses-1][iby_zj],2)," + ",round(yields_4e[nMasses-1][ibeu_zj],2)," - ",round(yields_4e[nMasses-1][ibed_zj],2)," & ",round(yields_4mu[nMasses-1][iby_zj],2)," + ",round(yields_4mu[nMasses-1][ibeu_zj],2)," - ",round(yields_4mu[nMasses-1][ibed_zj],2)," & ",round(yields_2e2mu[nMasses-1][iby_zj],2)," + ",round(yields_2e2mu[nMasses-1][ibeu_zj],2)," - ",round(yields_2e2mu[nMasses-1][ibed_zj],2)," \\\\"

    print "\\hline"
    ## print all backgrounds    
    if totbkgerrUp_4e == totbkgerrDn_4e and totbkgerrUp_4mu == totbkgerrDn_4mu and totbkgerrUp_2e2mu == totbkgerrDn_2e2mu:
        print "All background expected & ",round(totbkg_4e,2)," \pm ",round(totbkgerrUp_4e,2)," & ",round(totbkg_4mu,2)," \pm ",round(totbkgerrUp_4mu,2)," & ",round(totbkg_2e2mu,2)," \pm ",round(totbkgerrUp_2e2mu,2)," \\\\"
    else:
        print "All background expected & ",round(totbkg_4e,2)," + ",round(totbkgerrUp_4e,2)," - ",round(totbkgerrDn_4e,2)," & ",round(totbkg_4mu,2)," + ",round(totbkgerrUp_4mu,2)," - ",round(totbkgerrDn_4mu,2)," & ",round(totbkg_2e2mu,2)," \pm "," + ",round(totbkgerrUp_2e2mu,2)," - ",round(totbkgerrDn_2e2mu,2)," \\\\"    
    if UBbkgerrUp_4e == UBbkgerrDn_4e and UBbkgerrUp_4mu == UBbkgerrDn_4mu and UBbkgerrUp_2e2mu == UBbkgerrDn_2e2mu:
        print "All background (140 < $m_{4\\ell}$} < 300 GeV) & ",round(UBbkg_4e,2)," \pm ",round(UBbkgerr_4e,2)," & ",round(UBbkg_4mu,2)," \pm ",round(UBbkgerr_4mu,2)," & ",round(UBbkg_2e2mu,2)," \pm ",round(UBbkgerr_2e2mu,2)," \\\\"
    else:
        print "All background (140 < $m_{4\\ell}$} < 300 GeV) & ",round(UBbkg_4e,2)," + ",round(UBbkgerrUp_4e,2)," - ",round(UBbkgerrDn_4e,2)," & ",round(UBbkg_4mu,2)," + ",round(UBbkgerrUp_4mu,2)," - ",round(UBbkgerrDn_4mu,2)," & ",round(UBbkg_2e2mu,2)," + ",round(UBbkgerrUp_2e2mu,2)," - ",round(UBbkgerrDn_2e2mu,2)," \\\\"

    ## print signals
    print "\\hline"    
    for i in range(nMasses):
        print "$m_H$ = ",int(massToPrint[i]),"GeV & ",round(yields_4e[i][isy_cv],2)," \pm ",round(yields_4e[i][isy_eu],2)," & ",round(yields_4mu[i][isy_cv],2)," \pm ",round(yields_4mu[i][isy_eu],2)," & ",round(yields_2e2mu[i][isy_cv],2)," \pm ",round(yields_2e2mu[i][isy_eu],2)," \\\\" 

    ###################

    fout = open(options.textFile,'w')

    print >> fout, yields_4e[nMasses-1][iby_qq]," ",yields_4e[nMasses-1][ibeu_qq]," ",yields_4e[nMasses-1][ibeu_qq]," ", yields_4mu[nMasses-1][iby_qq]," ",yields_4mu[nMasses-1][ibeu_qq]," ",yields_4mu[nMasses-1][ibeu_qq]," ", yields_2e2mu[nMasses-1][iby_qq]," ",yields_2e2mu[nMasses-1][ibeu_qq]," ",yields_2e2mu[nMasses-1][ibeu_qq]

    print >> fout, yields_4e[nMasses-1][iby_zj]," ",yields_4e[nMasses-1][ibeu_zj]," ",yields_4e[nMasses-1][ibed_zj]," ",yields_4mu[nMasses-1][iby_zj]," ",yields_4mu[nMasses-1][ibeu_zj]," ",yields_4mu[nMasses-1][ibed_zj]," ",yields_2e2mu[nMasses-1][iby_zj]," ",yields_2e2mu[nMasses-1][ibeu_zj]," ",yields_2e2mu[nMasses-1][ibed_zj]

    print >> fout, totbkg_4e," ",totbkgerrUp_4e," ",totbkgerrDn_4e," ",totbkg_4mu," ",totbkgerrUp_4mu," ",totbkgerrDn_4mu," ",totbkg_2e2mu," ",totbkgerrUp_2e2mu," ",totbkgerrDn_2e2mu

    print >> fout, UBbkg_4e," ",UBbkgerrUp_4e," ",UBbkgerrDn_4e," ",UBbkg_4mu," ",UBbkgerrUp_4mu," ",UBbkgerrDn_4mu," ",UBbkg_2e2mu," ",UBbkgerrUp_2e2mu," ",UBbkgerrDn_2e2mu

    for i in range(nMasses):
        print >> fout,yields_4e[i][isy_cv]," ",yields_4e[i][isy_eu]," ",yields_4e[i][isy_ed]," ",yields_4mu[i][isy_cv]," ",yields_4mu[i][isy_eu]," ",yields_4mu[i][isy_ed]," ",yields_2e2mu[i][isy_cv]," ",yields_2e2mu[i][isy_eu]," ",yields_2e2mu[i][isy_ed] 









