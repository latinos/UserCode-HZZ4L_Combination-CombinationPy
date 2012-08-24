from optparse import OptionParser
from glob import glob
import json
import re

parser = OptionParser(usage="usage: %prog [options] datacard.txt [json]")
parser.add_option("-m", "--metric", dest="metric", default="all",  type="string", help="Metric: exp, obs, all")
parser.add_option("-T", "--td", "--threshold-default", dest="threshDef", default=0.01, type="float", help="Default threshold")
parser.add_option("-t", "--threshold",                 dest="threshCh",  default=[], action="append", type="string", nargs=2, help="Threshold for individual channel or set of channels")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

threshList = []
for lbl, thresh in options.threshCh:
    threshList.append( (re.compile(lbl), float(thresh)) )

metric_exp = lambda map : abs(map['exp']['AllIn']);
metric_obs = lambda map : abs(map['obs']['AllIn']);
metric_all = lambda map : max(metric_exp(map), metric_obs(map))
metrics = { 'exp':metric_exp, 'obs':metric_obs, 'all':metric_all };
metric = metrics[options.metric]

cleanName = args[0].replace(".txt","")+".clean.txt"
datacard = open(args[0], "r")
lines = [ l for l in datacard]
txtdatacards = {}
if lines[0].startswith("Combination of "):
    for X in lines[0].strip().split(", "):
        if "=" in X: txtdatacards[X.split("=")[1]] = True
print txtdatacards

toKeep, toDrop = {}, {}
for J in glob("*.txt.d/Removed1.json"):
    if len(txtdatacards) and J.replace(".d/Removed1.json","") not in txtdatacards: 
        print "skip",J.replace(".d/Removed1.json","")
        continue
    print "Reading json ",J
    report = json.loads(" ".join([l for l in open(J,"r")]))
    if not report: raise RuntimError, "Couldn't load %s" % args[0]
    threshold = options.threshDef
    for relbl, thresh in threshList:
        if re.search(relbl, J): threshold = thresh
    print " ... will apply threshold ",threshold
    for (nuisList, map) in report:
        outcome = metric(map)
        if outcome < threshold:
            print "   ... candidate for exlusion: %s of effect %.2f%%" % (", ".join(nuisList), 100*outcome)
            for L in nuisList: toDrop[L] = True
        else:
            print "   ... must be kept: %s of effect %.2f%%" % (", ".join(nuisList), 100*outcome)
            for L in nuisList: toKeep[L] = True

toReallyDrop = []
for N in toDrop.keys():
    if N not in toKeep: 
        toReallyDrop.append(N)
toReallyDrop.sort()
print "Ok, so the nuisances to drop are:\n"+("\n".join(["   * %s" % N for N in toReallyDrop])) + "\n";

nin, nout = (0,0)
cleanout = open(cleanName, "w")
cleanout.write("# from %s excluding these nuisances: %s" % (args[0], ", ".join(toReallyDrop)))
for line in lines:
    if re.match("^kmax\s+",line): line = "kmax *\n";
    m = re.match("^\s*(\w+)\s+(lnN|shape|gmN)\s+.*", line)
    if m:
        if m.group(1) in toReallyDrop:
            line = "#"+line
            nout+=1
        else:
            nin+=1
    cleanout.write(line)
print "Commented out ",nout," out of ",(nin+nout)," nuisances ===> ",cleanName
            
