TString channelFromName(TString who, bool withSpaces=false, bool noLumi=true) {
    TString name = who, space, lumi;
    int e = (who.Contains("_7") ? 7 : (who.Contains("_8") ? 8 : 0));
    if (who.Contains("comb"))  { name = "Combined"; space = "", lumi=(e == 0 ? "5.0+1.6 fb^{-1}" : (e == 7 ? "5.0 fb^{-1}" : "1.6 fb^{-1}")); }
    if (who.Contains("hzz"))  { 
        name = withSpaces ? "H #rightarrow ZZ" : "H #rightarrow ZZ"; 
        space = "              "; lumi= "4.7 fb^{-1}"; 
    }
    if (who.Contains("combp") || who.Contains("hires"))  { 
        name = "H #rightarrow ZZ + #gamma#gamma"; 
    }
    if (who.Contains("combl") || who.Contains("lowres"))  { 
        name = "H #rightarrow bb + #tau#tau + WW" ; 
    }
    if (who.Contains("fermions"))  name = "H #rightarrow bb + #tau#tau "; 
    if (who.Contains("bosons"))    name = "H #rightarrow ZZ + WW + #gamma#gamma" ; 

    if (who.Contains("combs"))  { 
        name = isSquareCanvas ? "VBF+VH excl." : "VBF+VH exclusive" ; 
        space = "   "; lumi= "4.6-4.8 fb^{-1}"; 
        if (withSpaces) { space = ""; lumi = ""; }
    }
    if (who.Contains("hww"))   { name = "H #rightarrow WW"; space = "           "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hww2l"))   { name = "H #rightarrow WW #rightarrow 2l 2#nu"; space = ""; lumi="4.6 fb^{-1}"; }
    if (who.Contains("vhww3l"))   { name = "WH #rightarrow 3l 3#nu"; space = "     "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hbb"))  { name = "H #rightarrow bb"; space = "             "; lumi="4.7 fb^{-1}"; }
    if (who.Contains("hgg"))   { name = "H #rightarrow #gamma#gamma"; space = "              "; lumi="4.8 fb^{-1}"; }
    if (who.Contains("htt"))   { name = "H #rightarrow #tau#tau"; space = "              "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("vhtt"))   { name = "VH #rightarrow #tau_{h} 2l"; space = "         "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("httm"))   { name = "H #rightarrow #tau#tau #rightarrow #mu#mu"; space = "    "; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hzz4l")) { name = "H #rightarrow ZZ #rightarrow 4l"; space = "     "; lumi="4.7 fb^{-1}"; }
    if (who.Contains("hzz2l2q"))  { name = "H #rightarrow ZZ #rightarrow 2l 2q"; space=""; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hzz2l2nu")) { name = "H #rightarrow ZZ #rightarrow 2l 2#nu"; space=""; lumi="4.6 fb^{-1}"; }
    if (who.Contains("hzz2l2t")) { name = "H #rightarrow ZZ #rightarrow 2l 2#tau"; space=""; lumi="4.6 fb^{-1}"; }
    if (who.Contains("ggH")) {  name = "Inclusive"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("qqH")) {  name = "VBF modes"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("VH"))  {  name = "VH modes"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("ttH")) {  name = "ttH modes"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("VH_hbb"))   { name = "H #rightarrow bb (VH)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("ttH_hbb"))   { name = "H #rightarrow bb (ttH)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("ggH_hgg"))   { name = "H #rightarrow #gamma#gamma (inc.)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("qqH_hgg"))   { name = "H #rightarrow #gamma#gamma (VBF)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("ggH_hww"))   { name = "H #rightarrow WW (0/1j)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("qqH_hww"))   { name = "H #rightarrow WW (2j)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("VH_hww"))    { name = "H #rightarrow WW (VH)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("ggH_htt"))   { name = "H #rightarrow #tau#tau (inc.)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("qqH_htt"))   { name = "H #rightarrow #tau#tau (VBF)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("VH_htt"))    { name = "H #rightarrow #tau#tau (VH)"; space = "   "; lumi = "? fb^{-1}"; }
    if (who.Contains("0")) name += " (old)";

    if (lessSpacesInLegends) {
        int nsp = space.Length()-12; space = "";
        for (int i = 0; i < nsp; ++i) space += " ";
    }
    if (isMultiEnergy) {
        if (who.Contains("_7")) { name = "#sqrt{s} = 7 TeV"; lumi  = "5.0 fb^{-1}"; }
        else if (who.Contains("_8")) { name = "#sqrt{s} = 8 TeV"; lumi  = "1.6 fb^{-1}"; }
    } else if (isMultiHPA) {
        if (who.Contains("_HPA")) { name = "HPA only"; }
        if (who.Contains("_HPA_clean")) { name = "HPA cleaned"; }
    }

    if (noLumi) return name;
    if (withSpaces) return (name == "Combined" || lumi == "") ? name : (name + space + "  ("+lumi+")");
    if (name == "Combined" && justLumiForCombined) return lumiSymbol+" = "+lumi;
    return name + ", "+lumiSymbol+" = "+lumi;
}

int colorFromName(TString who, bool dark=false) {
    if (isMultiEnergy) {
        if (who.Contains("_7")) return 4;
        if (who.Contains("_8")) return 206;
        return 1;
    } else if (isMultiHPA) {
        if (who.Contains("_HPA_clean")) return 214;
        if (who.Contains("_HPA")) return 99;
        return 1;
    }

    if (who.Contains("hires")) return 62;
    if (who.Contains("lowres")) return 98;
    if (who.Contains("comb"))  return (dark ? 19 : 1);

    if (who.Contains("fermions"))  return 98; 
    if (who.Contains("bosons"))    return 62; 

    // split by prod. and decay
    if (who.Contains("VH_hbb"))  return 201;
    if (who.Contains("ttH_hbb"))  return 28;
    if (who.Contains("qqH_hgg"))  return  62;
    if (who.Contains("ggH_hgg"))  return 209;
    if (who.Contains("ggH_hww"))  return 4;
    if (who.Contains("qqH_hww"))  return 213;
    if (who.Contains("VH_hww"))  return 62;
    if (who.Contains("ggH_htt"))  return 221;
    if (who.Contains("qqH_htt"))  return 52;
    if (who.Contains("VH_htt"))  return 223;
    // split by decay
    if (who.Contains("hww"))   return (dark ? 4 : 4);
    if (who.Contains("hgg"))   return (dark ? 209 : 209);
    if (who.Contains("hbb"))  return (dark ? 67 : 67 );
    if (who.Contains("htt"))   return (dark ? 221 : 221);
    if (who.Contains("hzz4l"))  return (dark ? 2 : 2);
    if (who.Contains("hzz2l2q"))  return (dark ? 93 : 93);
    if (who.Contains("hzz2l2nu")) return (dark ? 28 : 28);
    if (who.Contains("hzz2l2t"))  return (dark ? 223 : 223);
    if (who.Contains("hzz")) return (dark ? 205 : 2);
    return 39;
}
int lineStyleFromName(TString who) {
    if (noLineStyles) return 1;
    if (who.Contains("combzz"))   return 1;
    if (who.Contains("comb"))     return 1;
    if (who.Contains("vhbb"))     return 3;
    if (who.Contains("hww"))      return 7;
    if (who.Contains("hgg"))      return 5;
    if (who.Contains("htt"))      return 2;
    if (who.Contains("hzz4l"))    return 1;
    if (who.Contains("hzz2l2q"))  return 5;
    if (who.Contains("hzz2l2nu")) return 2;
    if (who.Contains("hzz2l2t"))  return 3;
    return 39;
}

