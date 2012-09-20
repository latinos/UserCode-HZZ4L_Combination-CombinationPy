#! /usr/bin/env python
import os
import re
import math
import ROOT
from array import array
import commands
import subprocess 
import glob
from optparse import OptionParser

ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
from ROOT import setTDRStyle
ROOT.setTDRStyle(True)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0000)

m_curdir=commands.getoutput("pwd")

parser = OptionParser()

parser.add_option('--stem',
                  action='store', type='string', dest='stem',default='_testXX_testYY',
                  help='string: stem of file')
parser.add_option('--seed',
                  action='store', type='int', dest='seed',default='12345',
                  help='int: seed')
parser.add_option('--toysPerJob',
                  action='store', type='int', dest='toys',default='5',
                  help='int: toys per job')
(options, args) = parser.parse_args()

#################### MAIN ####################

##_testPS_testSM_12344_0_maxllfit.root

if __name__ == '__main__':

    stem=options.stem
    seed=options.seed
    toys=options.toys
    
    print 'Adding files for LandS'
    fo = ROOT.TFile( stem+"_"+str(seed)+"_maxllfit.root", "recreate" )
    to = ROOT.TTree("T", "T")    
    limit_ = array( 'f', [ 0. ] )
    rmean_ = array( 'f', [ 0. ] )
    to.Branch("limit", limit_ , "limit/F")
    to.Branch("rmean", rmean_ , "rmean/F")

   # print 'Adding this listname: ',listnames

    for jj in range(toys):
     #   seedtot=seedd+offset
        #    if(jj==1): print seedtot,"  ",
        curfile = stem+"_"+str(seed)+"_"+str(jj)+"_maxllfit.root"
    
###loop over all available files    
###    searchPath=glob.glob( outputDir+"/"+listnames+"*_[0-9]*_maxllfit.root" )
###    #    print 'searchPath (1): ',searchPath  #os.path.join(outputDir,"/",listnames, '*_maxllfit.root')
###    cmd = "ls -l "+outputDir+"/"+listnames+"*_[0-9]*_maxllfit.root"
###    for curfile in searchPath:
###        print "Loading file: ",curfile

###list of files to exclude:
###            if (seedtot==12370 or seedtot==12371 or seedtot==12372 or seedtot==12373): continue
        fcur = ROOT.TFile( curfile )
        tcur = fcur.Get("T")
        nentries = int( tcur.GetEntries() )
        if(nentries==0): print "File with seed=",seedtot," has 0 entries"
        for ientry in range(nentries):
            tcur.GetEntry(ientry)
            limit_[0] = tcur.limit
            rmean_[0] = tcur.rmean
            #print "Limit and rmean: ",limit_[0],"  ",rmean_[0]
            to.Fill()
            fcur.Close()
            #end loop on i
        #end loop on jj
        
    print "haddLands: final tree has ",to.GetEntries()," entries"    
    fo.cd()
    to.Write()
    fo.Close()

    
