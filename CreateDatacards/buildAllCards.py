#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2012.08.30
# by Matt Snowball
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import math
from array import array

#define function for processing of os command
def processCmd(cmd):
#    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        sys.exit()
        
        

def creationLoop():

    cmd = 'mkdir -p buildLogs/1D'
    processCmd(cmd)
    cmd = 'mkdir -p buildLogs/2D'
    processCmd(cmd)
    
    models=['gg0-']

    print '++++++++++++++++++ Building Cards ++++++++++++++++++'

    for model in models:

        print '>>>> ',model,'2D'

        cmd = 'python makeDCsandWSs.py -b -i SM_inputs_7TeV -a 2D_7TeV_{0} -t templates2D_{0}/7TeV >& buildLogs/2D/{0}_7TeV.txt'.format(model)
        processCmd(cmd)
        cmd = 'python makeDCsandWSs.py -b -i SM_inputs_8TeV -a 2D_8TeV_{0} -t templates2D_{0}/8TeV >& buildLogs/2D/{0}_8TeV.txt'.format(model)
        processCmd(cmd)
        cmd = 'mkdir -p HCG_2D_{0}; cp cards_2D_7TeV_{0}/HCG/126/* HCG_2D_{0};cp cards_2D_8TeV_{0}/HCG/126/* HCG_2D_{0}'.format(model)
        processCmd(cmd)

        print '>>>> ',model,'1D'

        cmd = 'python makeDCsandWSs.py -b -i SM_inputs_7TeV -a 1D_7TeV_{0} -t templates2D_{0}/7TeV --unfold >& buildLogs/1D/{0}_7TeV.txt'.format(model)
        processCmd(cmd)
        cmd = 'python makeDCsandWSs.py -b -i SM_inputs_8TeV -a 1D_8TeV_{0} -t templates2D_{0}/8TeV --unfold >& buildLogs/1D/{0}_8TeV.txt'.format(model)
        processCmd(cmd)
        cmd = 'mkdir -p HCG_1D_{0}; cp cards_1D_7TeV_{0}/HCG/126/* HCG_1D_{0};cp cards_1D_8TeV_{0}/HCG/126/* HCG_1D_{0}'.format(model)
        processCmd(cmd)
        



if __name__ == "__main__":
    creationLoop()
