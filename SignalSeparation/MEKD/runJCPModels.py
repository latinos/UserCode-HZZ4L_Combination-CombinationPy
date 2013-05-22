#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2012.05.21
# by Predrag Milenovic,Matt Snowball
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import math
import ROOT




cmssw_base=commands.getoutput("echo $CMSSW_BASE")
curdir=commands.getoutput("pwd")
TOYSPERDIR=50000

def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-b', action='store_true', dest='noX', default=True, help='no X11 windows')
    parser.add_option('', '--fitNuis', action='store_true', dest='FITNUIS',  default=False, help='fit nuissances (default=False)')
    parser.add_option('', '--generateToys',action='store_true', dest='generateToys',default=False,help='generate toys')
    parser.add_option('', '--haddToys',action='store_true', dest='haddToys',default=False,help='hadd toys')
    parser.add_option('', '--plot',action='store_true', dest='plot',default=False,help='make separation plot')
    parser.add_option('', '--tool',  dest='TOOL', type='string', default='combine',    help='Tool: combine or lands')
    parser.add_option('-t', '--toys',  dest='NTOYS', type='int', default=1000000,    help='Number of total toys, will be paralelised to have 50K per job')
    parser.add_option('-d', '--dir',   dest='SOURCEDIR', type='string', default='', help='SOURCEDIR, skip if SOURCEDIR is empty')
    parser.add_option('-n', '--name',  dest='DIRNAME', type='string', default='submission_', help='submission dir names - submission_$i')
    parser.add_option('-M', '--model', dest='MODEL', type='string', default='', help='model name')
    parser.add_option('-m', '--mh',    dest='MASS', type='float', default=126.0, help='mass of Higgs hypothesis')
    parser.add_option('-e', '--sqrts', dest='ENERGY', type='string', default='7p8', help='sqrts: 7, 8, 7p8')
    parser.add_option('-u', '--mu',    dest='MUTYPE', type='string', default='fixed', help='type of mu: fixed[default] or float')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if opt.haddToys and opt.generateToys:
        print 'Can not choose hadd and generate toys! Please choose 1.'
        sys.exit()
        
    
        
    if opt.MODEL == '':
        print 'MODEL must be set (-M)'
        sys.exit()

    if opt.MUTYPE != 'float' and opt.MUTYPE != 'fixed':
        print 'Please choose float or fixed for --mu!'
        sys.exit()

    if opt.ENERGY != '7' and opt.ENERGY != '8' and opt.ENERGY != '7p8':
        print 'Please choose either 7, 8, or 7p8 for --sqrts/-e!'
        sys.exit()

    if opt.SOURCEDIR == '' and not (opt.haddToys or opt.plot) :
        print 'Please pass a source dir (-d) containing all cards and ws!'
        sys.exit()

    

#define function for processing of os command
def processCmd(cmd):
    #print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()
    return output


def runMixingModels():

    # parse the arguments and options
    global opt, args
    parseOptions()

    jcpDir = os.getcwd()

#    mixAnlges = [10,15,20,25,30,35,40,45,60,75,90]
#    mixAnlges = [100,105,110,115,120,125,130,135,140,145,150,155,160,165,170]
    mixAnlges = [10,15,20,25,30,35,40,45,60,75,90,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170]
#    mixAnlges = [45]

    strAngles = 'CP-mixing angles: '
    for angle in mixAnlges:
        strAngles = strAngles + str(angle) + ' '
    print '['+strAngles+']'

    # tag to indicate normal or PS calcs.
    KDappend = '_PS'

    print '[Preparing templates and running makeParallelToys]'
    for angle in mixAnlges:
        #create the datracards from the template in the SOURCEDIR unless name is an empty string
        cmd = 'cat ./sourceDir'+KDappend+'/StatCard_JCPmixAngle_template.txt | sed "s%MYANGLE%'+str(angle)+'%g" | sed "s%MYKDTYPE%'+KDappend+'%g" > ./sourceDir'+KDappend+'/StatCard_JCPmixAngle'+str(angle)+'.txt'
        processCmd(cmd)

        #get/prepare the templates from SOURCEDIR if name is not an empty string
        tmpDir = './sourceDir'+KDappend+'/JCPmixAngle'+str(angle)
        if not os.path.isdir(tmpDir):
            cmd = 'mkdir '+tmpDir
            processCmd(cmd)
        cmd = 'cp -r '+str(opt.SOURCEDIR)+'/Hist_1D_JCPmixAngle'+str(angle)+'_*.root '+tmpDir+'/'
        processCmd(cmd)

        # prepare parallel toys
        cmd = 'python makeParallelToys.py -t '+str(opt.NTOYS)+' -n submission_mix'+str(angle)+'_ -c StatCard_JCPmixAngle'+str(angle)+'.txt -d sourceDir'+KDappend
        processCmd(cmd)

    # run parallel toys for each model
    print '[Running toys]'
    for angle in mixAnlges:
        os.chdir(jcpDir+'/submission_mix'+str(angle)+'_1')
        print 'Running toys for angle '+str(angle)+'...'
        cmd = '. run.sh'
        processCmd(cmd)

    print '[Extracting stat. results]'
    for angle in mixAnlges:
        os.chdir(jcpDir+'/submission_mix'+str(angle)+'_1')
        #print 'Extracting stats for angle '+str(angle)+'...'
        cmd = '. extract_source.sh JCPmixAngle'+str(angle)
        output = processCmd(cmd)
        listOutput = list(shlex.shlex(output))
        print 'angle: '+str(angle)+', exp. sep.: ' + str(listOutput[len(listOutput)-8])+str(listOutput[len(listOutput)-7])+str(listOutput[len(listOutput)-6]) + ', CLs: ' + str(listOutput[len(listOutput)-3])+str(listOutput[len(listOutput)-2])+str(listOutput[len(listOutput)-1])




def combineCards(dir,sqrts):

    os.chdir(dir)
    
    if sqrts == '7':
        cmd = 'combineCards.py hzz4l_4mu_7TeV=hzz4l_4muS_7TeV_ALT.txt hzz4l_4e_7TeV=hzz4l_4eS_7TeV_ALT.txt hzz4l_2e2mu_7TeV=hzz4l_2e2muS_7TeV_ALT.txt > comb.txt'
        processCmd(cmd)
    elif sqrts == '8':
        cmd = 'combineCards.py hzz4l_4mu_8TeV=hzz4l_4muS_8TeV_ALT.txt hzz4l_4e_8TeV=hzz4l_4eS_8TeV_ALT.txt hzz4l_2e2mu_8TeV=hzz4l_2e2muS_8TeV_ALT.txt > comb.txt'
        processCmd(cmd)
    elif sqrts == '7p8':
        cmd = 'combineCards.py hzz4l_4mu_7TeV=hzz4l_4muS_7TeV_ALT.txt hzz4l_4e_7TeV=hzz4l_4eS_7TeV_ALT.txt hzz4l_2e2mu_7TeV=hzz4l_2e2muS_7TeV_ALT.txt hzz4l_4mu_8TeV=hzz4l_4muS_8TeV_ALT.txt hzz4l_4e_8TeV=hzz4l_4eS_8TeV_ALT.txt hzz4l_2e2mu_8TeV=hzz4l_2e2muS_8TeV_ALT.txt  > comb.txt'
        processCmd(cmd)
    else:
        raise RuntimeError,'Unknown sqrts! Please choose 7, 8, or 7p8.'

    os.chdir(curdir)



def makeBinaryWorkspace(dir,muType,mass):

    os.chdir(dir)
    
    if muType == 'fixed':
        cmd = 'text2workspace.py -m '+str(mass)+' comb.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o fixedMu.root'
        processCmd(cmd)
    elif muType == 'float':
        cmd = 'text2workspace.py -m '+str(mass)+' comb.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=muFloating -o floatMu.root'
    else:
        raise RuntimeError,'Unknown muType! Please choose fixed or float(not implemented yet).'

    os.chdir(curdir)
    

def makeParallelToys():
    global opt, args
    parseOptions()
    
    
    ###################
    SEED=1234
    SEEDINCR=12345
    ###################

    totalDirs = int(opt.NTOYS/TOYSPERDIR)
    print '----------------- Generate Parallel Toys -----------------'
    print '>>>>> '+str(totalDirs)+' total directories to be created'
    print '>>>>> '+str(opt.NTOYS)+' total toys to be run'
    print '>>>>> '+str(TOYSPERDIR)+' toys per directory'
    print '>>>>> combining cards for '+str(opt.ENERGY)+' TeV'
    combineCards(opt.SOURCEDIR,opt.ENERGY)
    print '>>>>> making binary workspace'
    makeBinaryWorkspace(opt.SOURCEDIR,opt.MUTYPE,opt.MASS)
    

    for i in xrange(1,totalDirs+1):
        theDir = opt.MODEL+'/'+opt.DIRNAME+str(i)

        cmd = 'mkdir -p '+str(theDir)
        #print cmd
        processCmd(cmd)
        
        cmd = 'cp '+str(opt.SOURCEDIR)+'/* '+str(theDir)+'/'
        #print cmd
        processCmd(cmd)

    return totalDirs

    


def submitToysLSF(totalDirs):
    global opt, args
    parseOptions()

    ###################
    SEED=1234
    SEEDINCR=12345
    MYFITNUIS = 'FALSE'
    if opt.FITNUIS:
        MYFITNUIS = 'TRUE'
    ###################
    CARD = ''
    if opt.MUTYPE == 'fixed':
        CARD = 'fixedMu.root'
    else:
        CARD = 'floatMu.root'

    print '>>>>> submitting toys'
    
    for i in xrange(1,totalDirs+1):
        theDir = curdir+'/'+opt.MODEL+'/'+opt.DIRNAME+str(i)
        cmd = "bsub -q 8nh -o "+str(theDir)+"/lsfLog.txt -J "+str(opt.MODEL)+"_"+str(i)+" run_source.sh "+str(theDir)+" "+str(opt.MUTYPE)+" "+str(MYFITNUIS)+" "+str(TOYSPERDIR)+" "+str(opt.MASS)+" "+str(SEED)+" "+CARD
        #print cmd
        processCmd(cmd)
        SEED=SEED+SEEDINCR





def haddToysCombine():
    global opt, args
    parseOptions()
    
    print '>>>>> hadding toys from combine' 

    os.chdir(opt.MODEL)
    myFile = 'floatMu'
    myNuis = ''
    
    if opt.MUTYPE == 'fixed':
        myFile = 'fixedMu'
    
    if opt.FITNUIS:
        myNuis = 'fitNuis'
    
    if not opt.FITNUIS:
        cmd = 'hadd combineToys.'+myFile+'.'+str(opt.MASS)+'.root submission_*/higgsCombine'+myFile+'.HybridNew.mH*root' 
        processCmd(cmd)
    else:
        cmd = 'hadd combineToys.'+myFile+'.'+str(opt.MASS)+'.root submission_*/higgsCombine'+myFile+'.'+myNuis+'.HybridNew.mH*root' 
        processCmd(cmd)
    
    os.chdir(curdir)

    filename = 'combineToys.'+myFile+myNuis+'.'+str(opt.MASS)+'.root'
    
    print '>>>>> hadded toy file name: '+filename

    return filename


def plotToysCombine(file):
    global opt, args
    parseOptions()


    cmd = 'cp extractSignificanceStats.C '+opt.MODEL
    processCmd(cmd)
    cmd = 'cp tdrstyle.C '+opt.MODEL
    processCmd(cmd)
    os.chdir(opt.MODEL)

    myNuis = ''
    if opt.FITNUIS:
        myNuis = 'fitNuis'

    print '>>>>> extracting toys '
    
    cmd = 'root -q -b '+file+' "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx(\\"'+opt.MODEL+'.qmu.'+opt.MUTYPE+myNuis+'.root\\",'+str(opt.MASS)+',1,\\"x\\")"'
    processCmd(cmd)
    cmd = 'cp '+opt.MODEL+'.qmu.'+opt.MUTYPE+myNuis+'.root qmu.root'
    processCmd(cmd)


    print '>>>>> plotting '+opt.MODEL
    ROOT.gROOT.ProcessLine(".L extractSignificanceStats.C")
    cmd = 'extractSignificanceStats("'+opt.MODEL+'_'+opt.MUTYPE+myNuis+'"); > '+opt.MODEL+'_'+opt.MUTYPE+myNuis+'.log'
    #cmd = 'extractSignificanceStats()'
    ROOT.gROOT.ProcessLine(cmd)

    #print '>>>>> results: '
    #cmd = 'cat '+opt.MODEL+'_'+opt.MUTYPE+myNuis+'.log'
    #output = processCmd(cmd)
    #print output
    os.chdir(curdir)





# run the create_RM_cfg() as main()
if __name__ == "__main__":
    global opt, args
    parseOptions()
    
    print '++++++++++++++++++++++++ MODEL: '+str(opt.MODEL)+' ++++++++++++++++++++++++'
    
    if opt.generateToys:
        cmd = 'mkdir -p '+str(opt.MODEL)
        if os.path.exists(opt.MODEL):
            print opt.MODEL+' exists! Exiting...'
            sys.exit()
        else:
            processCmd(cmd)
                    
        totalDirs = makeParallelToys()
        submitToysLSF(totalDirs)

    elif opt.haddToys:
        if not os.path.exists(opt.MODEL):
            print opt.MODEL+' does not exist!'
            sys.exit()
            
        if opt.TOOL == 'combine':
            haddToysCombine()
        else:
            raise RuntimeError,'not implemented for lands yet'

    elif opt.plot:
        if not os.path.exists(opt.MODEL):
            print opt.MODEL+' does not exist!'
            sys.exit()
            
        if opt.TOOL == 'combine':
            outputFile = haddToysCombine()
            plotToysCombine(outputFile)
        else:
            raise RuntimeError,'not implemented for lands yet'

    else:
        raise RuntimeError,'No method chosen!'
    

        
    sys.exit()
        
