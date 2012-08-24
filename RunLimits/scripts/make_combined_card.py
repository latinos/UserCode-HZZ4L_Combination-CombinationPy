#! /usr/bin/env python
import sys, re, os, os.path

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] mass")
parser.add_option("-p", "--production", dest="production",  help="production topologies (regexp)", default="", metavar="PATTERN")
parser.add_option("-d", "--decay",      dest="decay",       help="decay modes (regexp)",           default="", metavar="PATTERN")
parser.add_option("-e", "--energy",     dest="energy",      help="energy (7,8,0=all)",             type="int", default="0", metavar="SQRT(S)")
parser.add_option("-H", "--HPA",        dest="hpa",         help="include only HPAs + a few other high mass modes driving the combination", default=False, action="store_true")
parser.add_option("-v", "--verbose",    dest="verbose",     help="list the datacards that will go into this combination", default=False, action="store_true")
parser.add_option("-n", "--dry-run",    dest="pretend",     help="(use with -v) just list the datacards that will go into this combination", default=False, action="store_true")
parser.add_option("-z", "--zip",        dest="zip",         help="compress output datacard", default=False, action="store_true")
parser.add_option("-o", "--out",        dest="outname",     help="output name",              default=None, metavar="PATTERN")
parser.add_option("-L", "--loose",      dest="loose",       help="turn off strict validation", default=False, action="store_true")
parser.add_option("-x", "--exclude",    dest="exclude",     help="exclude these datacards (regexp)", default=[], action="append")
parser.add_option("-s", "--select",     dest="select",      help="select these datacards (regexp)", default=[], action="append")
parser.add_option("--xnf", "--exclude-nuisance-file", dest="nuisVetoFile", help="exclude the nuisances in this file", default=None)
(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    sys.exit(2)
mass = float(args[0])

class Filter:
    def __init__(self,select,exclude):
        self._select = [re.compile(p) for p in select]
        self._exclude= [re.compile(p) for p in exclude]
    def __call__(self,obj):
        for p in self._exclude:
            if re.search(p, obj): return False
        for p in self._select:
            if re.search(p, obj): return True
        return len(self._select) == 0

class Channel:
    def __init__(self,production,decay,subdecay,energy,massrange,hpa,cards):
        self.production = production
        self.decay      = decay
        self.subdecay   = subdecay
        self.energy     = energy
        self.massrange  = massrange
        self.hpa        = hpa
        self.cards      = cards
    def listCards(self,mass, filter = lambda x : True):
        ret = []
        for C,D in self.cards.items():
            if not filter(D): continue
            CI = "%s%s_%d_%s" % (self.decay, self.subdecay, self.energy, C)
            if '%' in D: ret.append( (CI, D % mass) )
            else:        ret.append( (CI, D) )
        return ret
    def prepareCards(self,mass, filter = lambda x : True):
        for C,D in self.cards.items():
            if not filter(D): continue
            if '%' in D: D = D % mass 
            if not os.path.exists("%g/%s" % (mass,D)):
                raise RuntimeError, "Missing datacard %g/%s" % (mass,D)
    def __str__(self):
        return "%s -> %s -> %s at %d TeV" % (self.production, self.decay, self.subdecay, self.energy)
class CommonChannel(Channel):
    def __init__(self,production,decay,subdecay,energy,massrange,hpa,cards):
        Channel.__init__(self,production,decay,subdecay,energy,massrange,hpa,cards)
    def prepareCards(self,mass, filter = lambda x : True):
        for C,D in self.cards.items():
            if not filter(D): continue
            if '%' in D: D = D % mass 
            if not os.path.exists("%g/%s" % (mass,D)):
                if os.path.exists("%s/%s" % ("common",D)):
                    os.system("cd %g && cp -v -s ../common/%s ." % (mass,D))
                else: raise RuntimeError, "Missing datacard %s" % D

CHANNELS = [
    # ============ H -> bb ==================
    Channel( "VH", "hbb","",    7,(110,135),True,  { "Wln":"vhbb_Wln_7TeV.txt", "Zll":"vhbb_Zll_7TeV.txt", "Znn":"vhbb_Znn_7TeV.txt" } ),
    Channel( "VH", "hbb","",    8,(110,135),True,  { "Wln":"vhbb_Wln_8TeV.txt", "Zll":"vhbb_Zll_8TeV.txt", "Znn":"vhbb_Znn_8TeV.txt" } ),

    CommonChannel( "ttH","hbb","ttH", 7,(110,140),False, { "":"ttH_7TeV.txt" } ),

    # ============ H -> tt ==================
    Channel( "ggH","htt","mm",  7,(110,145),True, { "0":"htt_mm_0_7TeV.txt", "1":"htt_mm_1_7TeV.txt", "2":"htt_mm_2_7TeV.txt", "2":"htt_mm_3_7TeV.txt" } ),
    Channel( "ggH","htt","em",  7,(110,145),True, { "0":"htt_em_0_7TeV.txt", "1":"htt_em_1_7TeV.txt", "2":"htt_em_2_7TeV.txt", "2":"htt_em_3_7TeV.txt" } ),
    Channel( "ggH","htt","et",  7,(110,145),True, { "0":"htt_et_0_7TeV.txt", "1":"htt_et_1_7TeV.txt", "2":"htt_et_2_7TeV.txt", "2":"htt_et_3_7TeV.txt" } ),
    Channel( "ggH","htt","mt",  7,(110,145),True, { "0":"htt_mt_0_7TeV.txt", "1":"htt_mt_1_7TeV.txt", "2":"htt_mt_2_7TeV.txt", "2":"htt_mt_3_7TeV.txt" } ),
    Channel( "ggH","htt","mm",  8,(110,145),True, { "0":"htt_mm_0_8TeV.txt", "1":"htt_mm_1_8TeV.txt", "2":"htt_mm_2_8TeV.txt", "2":"htt_mm_3_8TeV.txt" } ),
    Channel( "ggH","htt","em",  8,(110,145),True, { "0":"htt_em_0_8TeV.txt", "1":"htt_em_1_8TeV.txt", "2":"htt_em_2_8TeV.txt", "2":"htt_em_3_8TeV.txt" } ),
    Channel( "ggH","htt","et",  8,(110,145),True, { "0":"htt_et_0_8TeV.txt", "1":"htt_et_1_8TeV.txt", "2":"htt_et_2_8TeV.txt", "2":"htt_et_3_8TeV.txt" } ),
    Channel( "ggH","htt","mt",  8,(110,145),True, { "0":"htt_mt_0_8TeV.txt", "1":"htt_mt_1_8TeV.txt", "2":"htt_mt_2_8TeV.txt", "3":"htt_mt_3_8TeV.txt" } ),

    Channel( "qqH","htt","mm",  7,(110,145),True, { "5":"htt_mm_5_7TeV.txt" } ),
    Channel( "qqH","htt","em",  7,(110,145),True, { "5":"htt_em_5_7TeV.txt" } ),
    Channel( "qqH","htt","et",  7,(110,145),True, { "5":"htt_et_5_7TeV.txt" } ),
    Channel( "qqH","htt","mt",  7,(110,145),True, { "5":"htt_mt_5_7TeV.txt" } ),
    Channel( "qqH","htt","mm",  8,(110,145),True, { "5":"htt_mm_5_8TeV.txt" } ),
    Channel( "qqH","htt","em",  8,(110,145),True, { "5":"htt_em_5_8TeV.txt" } ),
    Channel( "qqH","htt","et",  8,(110,145),True, { "5":"htt_et_5_8TeV.txt" } ),
    Channel( "qqH","htt","mt",  8,(110,145),True, { "5":"htt_mt_5_8TeV.txt" } ),
 
    Channel( "VH", "htt","2lt", 7,(110,145),False, { "":"vhtt_0_7TeV.txt", } ),
    Channel( "VH", "htt","2l2t",7,(110,150),False, { "":"vhtt_1_7TeV.txt" } ),

    # ============ H -> gg ==================
    CommonChannel( "ggH", "hgg","",7,(110,150),True, { "inc": "hgg_inc_7TeV.txt" } ),
    CommonChannel( "ggH", "hgg","",8,(110,150),True, { "inc": "hgg_inc_8TeV.txt" } ),
    CommonChannel( "qqH", "hgg","",7,(110,150),True, { "vbf": "hgg_vbf_7TeV.txt" } ),
    CommonChannel( "qqH", "hgg","",8,(110,150),True, { "vbf": "hgg_vbf_8TeV.txt" } ),

    # ============ H -> WW ==================
    Channel( "ggH","hww","2l2v",7,(110,600),True,  { "of0j": "hwwof_0j_shape_7TeV.txt", "sf0j": "hwwsf_0j_shape_7TeV.txt", "of1j": "hwwof_1j_shape_7TeV.txt", "sf1j": "hwwsf_1j_shape_7TeV.txt" } ),
    Channel( "ggH","hww","2l2v",8,(110,600),True,  { "of0j": "hwwof_0j_cut_8TeV.txt", "sf0j": "hwwsf_0j_cut_8TeV.txt", "of1j": "hwwof_1j_cut_8TeV.txt", "sf1j": "hwwsf_1j_cut_8TeV.txt" } ),
    Channel( "qqH","hww","2l2v",7,(110,600),True,  { "2j":   "hww_2j_cut_7TeV.txt" } ),
    Channel( "qqH","hww","2l2v",8,(110,600),True,  { "sf2j": "hwwsf_2j_cut_8TeV.txt", "of2j": "hwwof_2j_cut_8TeV.txt" } ),
    Channel( "ggH","hww","lvjj",7,(170,600),True,  { "":"hwwlvjj_shape_7TeV.txt"} ),
    #Channel( "ggH","hww","lvjj",8,(180,200),True,  { "":"hwwlvjj_shape_8TeV.txt"} ),
    #Channel( "ggH","hww","lvjj",8,(300,600),True,  { "":"hwwlvjj_shape_8TeV.txt"} ),
    Channel( "ggH","hww","lvjj",8,(180,600),True,  { "":"hwwlvjj_shape_8TeV_5p1fbinv.txt"} ),
    Channel( "VH", "hww","3l",  7,(110,200),False, { "3l": "vh3l_cut_7TeV.txt" }),
    Channel( "VH", "hww","jj2l",7,(118,190),False, { "sf2ljj": "vhwwlnlnjj_OF_7TeV.txt", 'sf2ljj': "vhwwlnlnjj_SF_7TeV.txt" } ),

    # ============ H -> ZZ ==================
    Channel( "ggH","hzz","4l",  7,(110,600),True,  { "4mu":"hzz4l_4muS_7TeV.txt", "4e":"hzz4l_4eS_7TeV.txt", "2e2mu":"hzz4l_2e2muS_7TeV.txt" } ),
    Channel( "ggH","hzz","4l",  8,(110,600),True,  { "4mu":"hzz4l_4muS_8TeV.txt", "4e":"hzz4l_4eS_8TeV.txt", "2e2mu":"hzz4l_2e2muS_8TeV.txt" } ),

    Channel( "ggH","hzz","2l2q",  7,(130,164),False,  { "ee0b":"hzz2l2q_ee0b_7TeV.txt", "ee1b":"hzz2l2q_ee1b_7TeV.txt", "ee2b":"hzz2l2q_ee2b_7TeV.txt",
                                                        "mm0b":"hzz2l2q_mm0b_7TeV.txt", "mm1b":"hzz2l2q_mm1b_7TeV.txt", "mm2b":"hzz2l2q_mm2b_7TeV.txt" } ),
    Channel( "ggH","hzz","2l2q",  7,(200,600),True,  { "ee0b":"hzz2l2q_ee0b_7TeV.txt", "ee1b":"hzz2l2q_ee1b_7TeV.txt", "ee2b":"hzz2l2q_ee2b_7TeV.txt",
                                                       "mm0b":"hzz2l2q_mm0b_7TeV.txt", "mm1b":"hzz2l2q_mm1b_7TeV.txt", "mm2b":"hzz2l2q_mm2b_7TeV.txt" } ),

    Channel( "ggH","hzz","2l2v",  7,(200,600),True,  { "ee0j":"hzz2l2v_%g_7TeV_eeeq0jets.dat", "ee1j":"hzz2l2v_%g_7TeV_eeeq1jets.dat", "ee2j":"hzz2l2v_%g_7TeV_eegeq2jets.dat", 
                                                       "mm0j":"hzz2l2v_%g_7TeV_mumueq0jets.dat", "mm1j":"hzz2l2v_%g_7TeV_mumueq1jets.dat", "mm2j":"hzz2l2v_%g_7TeV_mumugeq2jets.dat" } ),
    Channel( "qqH","hzz","2l2v",  7,(200,600),True,  { "eevbf":"hzz2l2v_%g_7TeV_eevbf.dat", "mmvbf":"hzz2l2v_%g_7TeV_mumuvbf.dat" } ),
    Channel( "ggH","hzz","2l2v",  8,(200,600),True,  { "ee0j":"hzz2l2v_%g_8TeV_eeeq0jets.dat", "ee1j":"hzz2l2v_%g_8TeV_eeeq1jets.dat", "ee2j":"hzz2l2v_%g_8TeV_eegeq2jets.dat", 
                                                       "mm0j":"hzz2l2v_%g_8TeV_mumueq0jets.dat", "mm1j":"hzz2l2v_%g_8TeV_mumueq1jets.dat", "mm2j":"hzz2l2v_%g_8TeV_mumugeq2jets.dat" } ),
    Channel( "qqH","hzz","2l2v",  8,(200,600),True,  { "eevbf":"hzz2l2v_%g_8TeV_eevbf.dat", "mmvbf":"hzz2l2v_%g_8TeV_mumuvbf.dat" } ),

    Channel( "ggH","hzz","2l2t",  7,(200,600),True,  { "":"hzz2l2t.shape_7TeV.txt" } ),
    Channel( "ggH","hzz","2l2t",  8,(200,600),True,  { "":"hzz2l2t.shape_8TeV.txt" } ),
]

filter = Filter(options.select, options.exclude)

toCombine = []
if options.verbose: print "Considering the following analyses: "
for C in CHANNELS:
    if options.production != "" and re.search(options.production, C.production) == None: continue
    if options.decay != "" and re.search(options.decay, C.decay) == None: continue
    if options.energy != 0 and options.energy != C.energy: continue
    if options.hpa and not C.hpa: continue
    if not(C.massrange[0] <= mass and mass <= C.massrange[1]): continue
    if options.verbose: print "   in mode ",C
    if not options.pretend: C.prepareCards(mass, filter)
    for SC, D in C.listCards(mass, filter):
        if options.verbose:     print "      channel",SC,"\t datacard",D
        toCombine.append((SC,D))

if not(options.pretend) and len(toCombine):
    outname = options.outname
    if outname == None: 
        outname = "comb"
        if options.energy != 0: outname += str(options.energy)
        if options.production != "": outname += ("_%s" % options.production)
        if options.decay != "": outname += ("_%s" % options.decay)
        if options.hpa: outname += "_HPA"
        if options.nuisVetoFile: outname += "_clean" 
        outname += ".txt"
    pipe = ""
    if options.zip:
        pipe = " | gzip"; outname += ".gz"
    command =  "cd %g && " % mass
    command += "combineCards.py -S "
    if options.loose: command += " --X-no-jmax "
    if options.nuisVetoFile: command += " --xn-file=%s " % options.nuisVetoFile
    command += " ".join(["%s=%s" %(C,D) for (C,D) in toCombine])
    command += pipe
    command += " > "+outname
    os.system(command)
    print "Output saved to",outname
