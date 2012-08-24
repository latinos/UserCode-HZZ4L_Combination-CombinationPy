from glob import glob
import re, os.path

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] mass")
(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    sys.exit(2)
mass = float(args[0])
mprefix = "%g/" % mass
datacards = [ D for D in glob(mprefix+"*.txt") if (mprefix+"comb" not in D) ]

nuisances = {}
for D in datacards:
    basename = os.path.basename(D)
    nuisStart = False
    for l in [ l.strip() for l in open(D, "r") ]:
        if l.startswith("rate"): 
            nuisStart = True; continue
        elif nuisStart:
            cols = l.split()
            if len(cols) > 2:
                nuisName = cols[0]
                if nuisName not in nuisances: nuisances[nuisName] = []
                nuisances[nuisName].append(basename)
nuisorted = nuisances.keys()
nuisorted.sort()
print "| *nuisance* | *cards* |"
for n in nuisorted:
    print "| <b>%s</b>   | %s   |" % (n,  ", ".join(nuisances[n]))
        
