#!/usr/bin/python
import sys, os, pwd, commands
import optparse, shlex, re
import math
from array import array


def parseOptions():
    
    usage = ('usage: %prog [options] \n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-M', '--method', dest='method', type='string', default="",    help='type [ASCLS,PLP,PLPE,PLS,PLSE,ML]')
    parser.add_option('-f', '--massfile',   dest='massfile',  type='string',  default="masses_half.txt", help='mass file [masses_half.txt]')
    parser.add_option('-q', '--scheduler', dest='scheduler', type='string', default="lsf",    help='scheduler [pbs or lsf]')
    parser.add_option('-t', '--tool', type='string', dest='tool', default="combine" ,help='tool [lands or combine]')
    parser.add_option('-o', '--options',   dest='options',  type='string', default="",     help='options [-S for strict]')
    
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.method != "ASCLS" and opt.method != "PLP" and opt.method != "PLPE" and opt.method != "ML"  and opt.method != "PLS" and opt.method != "PLSE" ):
        print opt.method, " is not an optional method"
        print "Please choose [ASCLS,PLP,PLPE,PLS,PLSE,ML]"
        sys.exit()

    if not os.path.exists(opt.massfile):
        print opt.massfile, " does not exist"
        sys.exit()

    if (opt.scheduler != "pbs" and opt.scheduler != "lsf"):
        print opt.scheduler, " is not an option"
        print "Please choose [pbs or lsf]"
        sys.exit()

    
                
                
def processCmd(cmd):
    #    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()



def submit():
    global opt, args
    parseOptions()

    processCmd('mkdir -p errFiles')
    processCmd('mkdir -p outFiles')
    processCmd('mkdir -p results')

    for line in open(opt.massfile,'r'):
        f = line.split()
        if len(f) < 1: continue

        print f[0]
        
        cmd = ""
        if opt.scheduler == "lsf" :
            cmd = "bsub -q 1nh -o "+str(f[0])+"/lsflog_"+opt.method+".txt makeLimits.lsf.sh "+opt.method+" "+str(f[0])+" "+opt.tool+" "+opt.options
        elif opt.scheduler == "pbs" :
            cmd = "qsub makeLimits.pbs.sh -v MASS="+str(f[0])+",TYPE="+opt.method+",OPTIONS="+opt.options+",TOOL="+opt.tool+" -N \""+opt.method+"_"+str(f[0])+"\""
        processCmd(cmd)


    sys.exit()

if __name__ == "__main__":
    submit()



