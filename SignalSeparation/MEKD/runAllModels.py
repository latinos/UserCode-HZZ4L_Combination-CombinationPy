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
    numberOfToys = 1000000
    models=['gg0-']
    muTypes=['float']
    fitNuisTypes=[0,1]
    #########################

    for model in models:
        print '+++++++++++++ Model: '+model+' +++++++++++++'
        for muType in muTypes:
            print '>>>>> mu = '+muType
            for fitNuisType in fitNuisTypes:
                print '>>>>> fitNuis = '+str(fitNuisTypes)

                if fitNuisType == 1:
                    print '>>>> 1D'
                    cmd = 'python runJCPModels.py --generateToys -M {0} --mu={1} --fitNuis -t {2} -d HCG_1D_{0}'.format(model,muType,numberOfToys)
                    processCmd(cmd)
                    print '>>>> 2D'
                    cmd = 'python runJCPModels.py --generateToys -M {0} --mu={1} --fitNuis -t {2} -d HCG_2D_{0}'.format(model,muType,numberOfToys)
                    processCmd(cmd)
                else:
                    print '>>>> 1D'
                    cmd = 'python runJCPModels.py --generateToys -M {0} --mu={1} -t {2} -d HCG_1D_{0}'.format(model,muType,numberOfToys)
                    processCmd(cmd)
                    print '>>>> 2D'
                    cmd = 'python runJCPModels.py --generateToys -M {0} --mu={1} -t {2} -d HCG_2D_{0}'.format(model,muType,numberOfToys)
                    processCmd(cmd)





#define function for processing of os command
def processCmd(cmd):
    #print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()
    return output



if __name__ == "__main__":
    runAllModels()
