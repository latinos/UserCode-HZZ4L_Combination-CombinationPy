import os
import ROOT
higgsMass=126.0

for finalstate in ['4e','4mu','2e2mu']:
    for period in ['7TeV','8TeV']:
        #stri='hzz4l_PARAM_m'+str(higgsMass)+'_'+finalstate+'_'+period+'_'+category+'.root'
        stri='hzz4l_'+str(finalstate)+'S_'+str(period)+'.input.root'
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        f=ROOT.TFile(stri)
        f.cd()
        w=f.Get('w')
        w.var("CMS_zz4l_mass").setVal(126.)
        #w.var("").Print()
        #if category == 'vbf':
        #    w.var("H_Fisher").Print()
        #if category == 'novbf':
        #    w.var("H_PtOM").Print()
                
        print 'FINALSTATE----',finalstate+"_"+period    
        print 'GGH:',w.function("ggH_norm").getVal()
        print 'QQH',w.function("qqH_norm").getVal()
        print 'WH',w.function("WH_norm").getVal()
        print 'ZH',w.function("ZH_norm").getVal()
        print 'ttH',w.function("ttH_norm").getVal()
        
        f.Close()
    

