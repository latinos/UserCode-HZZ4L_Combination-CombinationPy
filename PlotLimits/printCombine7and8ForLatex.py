#! /usr/bin/env python
import os
import re
import math
import ROOT
from array import array

from optparse import OptionParser

file1 = open("test7TeV.txt")
file2 = open("test8TeV.txt")

values1 = []
for line in file1:
    curLine = line.split()
    values1.append( curLine )
values2 = []
for line in file2:
    curLine = line.split()
    values2.append( curLine )

print values1
print values2

combined = []

if len(values1) == len(values2):
        
    print len(values1),",",len(values2)
    for i in xrange(len(values1)):
        rowcombined = []
        print len(values1[i]),",",len(values2[i])
        for j in xrange(len(values1[i])):
            rowcombined.append( float(values1[i][j]) + float(values2[i][j]) )
        combined.append( rowcombined )

    print combined

    massToPrint = [120.,126.,130.,200.,350.,500.]
    nMasses = len(massToPrint)

    print "\n"
    print "\n"
    print "\n"

    print "ZZ background & ",round(combined[0][0],2)," \pm ",round(combined[0][1],2)," & ",round(combined[0][3],2)," \pm ",round(combined[0][4],2)," & ",round(combined[0][6],2)," \pm ",round(combined[0][7],2)," \\\\"
    
    print "Z + X & ",round(combined[1][0],2),"^{ + ",round(combined[1][1],2),"}_{ - ",round(combined[1][2],2),"} & ",round(combined[1][3],2),"^{ + ",round(combined[1][4],2),"}_{ - ",round(combined[1][5],2),"} & ",round(combined[1][6],2),"^{ + ",round(combined[1][7],2),"}_{ - ",round(combined[1][8],2),"} \\\\"
    
    print "\\hline"
    ## print all backgrounds    

    print "All background expected & ",round(combined[2][0],2),"^{ + ",round(combined[2][1],2),"}_{ - ",round(combined[2][2],2),"} & ",round(combined[2][3],2),"^{ + ",round(combined[2][4],2),"}_{ - ",round(combined[2][5],2),"} & ",round(combined[2][6],2),"^{ + ",round(combined[2][7],2),"}_{ - ",round(combined[2][8],2),"} \\\\"    

    print "All background (140 < $m_{4\\ell}$} < 300 GeV) & ",round(combined[3][0],2),"^{ + ",round(combined[3][1],2),"}_{ - ",round(combined[3][2],2),"} & ",round(combined[3][3],2),"^{ + ",round(combined[3][4],2),"}_{ - ",round(combined[3][5],2),"} & ",round(combined[3][6],2),"^{ + ",round(combined[3][7],2),"}_{ - ",round(combined[3][8],2),"} \\\\"
    
    ## print signals
    print "\\hline"    
    for i in range(nMasses):
        y = i+4
        print "$m_H$ = ",int(massToPrint[i]),"GeV & ",round(combined[y][0],2)," \pm ",round(combined[y][1],2)," & ",round(combined[y][3],2)," \pm ",round(combined[y][4],2)," & ",round(combined[y][6],2)," \pm ",round(combined[y][7],2)," \\\\" 

else:
    print "mismatch of event files!!!"
