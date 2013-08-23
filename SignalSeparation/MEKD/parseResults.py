#!/usr/bin/python
#-----------------------------------------------
# Latest update: 06.05.2013
# by Predrag Milenovic, Matt Snowball
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import math


def runAllModels():

    #########################
    outputFile='results.txt'
    plotDir = 'plots'
    models=['gg0-','gg2m+','qq1+','qq1-','qq2m+','gg0h+','2h+','2h-','2b+','PI1+','PI1-','PI2+']
    muTypes=['float','fixed']
    #########################

    data = []
    obsJP = []
    obsSM = []
    cls = []
    expMu1 = []
    exp = []
    
    for model in models:
        for muType in muTypes:
            if muType == 'float':
                fitNuisType = 'fitNuis'
            else:
                fitNuisType = ''
            
            theInput = model+'/'+model+'_'+muType+fitNuisType+'.log'
            for line in open(theInput,'r'):
                f = line.split()
                if len(f) < 1: continue
                if f[0].startswith("#"): continue
                if f[0].lower().startswith("keyword"):
                    
                    if f[1].lower().startswith("data"):
                        if muType == 'float':
                            data.append(float(f[2]))
                    if f[1].lower().startswith("obssm"):
                        if muType == 'float':
                            obsSM.append(float(f[2])*-1)
                    if f[1].lower().startswith("obsjp"):
                        if muType == 'float':
                            obsJP.append(float(f[2]))
                    if f[1].lower().startswith("cls"):
                        if muType == 'float':
                            cls.append(float(f[2]))
                    if f[1].lower().startswith("expsep"):
                        if muType == 'float':
                            exp.append(float(f[2]))
                        elif muType == 'fixed':
                            expMu1.append(float(f[2]))
                        
                                          

    if len(models) != len(data):
        print "len(data) != len(models)"
        sys.exit()
    if len(models) != len(exp):
        print "len(exp) != len(models)"
        sys.exit()
    if len(models) != len(expMu1):
        print "len(expMu1) != len(models)"
        sys.exit()
    if len(models) != len(obsJP):
        print "len(obsJP) != len(models)"
        sys.exit()
    if len(models) != len(obsSM):
        print "len(obsSM) != len(models)"
        sys.exit()


    f = open(outputFile,'w')
    f.write('model:  exp  (mu=1)    obsSM    obsJP    CLs   data\n')
    f.write('----------------------------------------------------\n')
    i=0
    for model in models:
        f.write('{0}:   {1:.2f}  ({2:.2f})   {3:.2f}   {4:.2f}   {5}   {6:.2f}\n'.format(model,exp[i],expMu1[i],obsSM[i],obsJP[i],cls[i],data[i]))
        i+=1
        cmd = 'cp '+model+'/*.eps '+plotDir
        processCmd(cmd)
        cmd = 'cp '+model+'/*.png '+plotDir
        processCmd(cmd)
    
                
#define function for processing of os command
def processCmd(cmd):
    #print cmd
    status, output = commands.getstatusoutput(cmd)
#    if status !=0:
#        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
#        sys.exit()
    return output



if __name__ == "__main__":
    runAllModels()
