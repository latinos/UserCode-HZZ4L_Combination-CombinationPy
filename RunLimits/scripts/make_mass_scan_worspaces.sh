cd 125;
#text2workspace.py comb_hzz.txt.gz -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=115,135 -o FloatMass_comb_hzz.root
text2workspace.py comb_hgg.txt.gz -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=120,130 -o FloatMass_comb_hgg.root
text2workspace.py comb_ggH_hgg.txt.gz -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=120,130 -o FloatMass_comb_ggH_hgg.root
text2workspace.py comb_qqH_hgg.txt.gz -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=115,135 -o FloatMass_comb_qqH_hgg.root
#text2workspace.py comb_hires.txt.gz -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=120,130 -o FloatMass_comb_hires.root
#text2workspace.py comb_hires.txt.gz -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO higgsMassRange=120,130 --PO verbose \
#        --PO 'map=hgg.*inc.*/(V|gg|qq|tt)H:r_ggInc[1,0,5]' \
#        --PO 'map=hgg.*vbf.*/(V|gg|qq|tt)H:r_ggVBF[1,0,10]' \
#        --PO 'map=hzz.*/([WZ]|gg|qq|tt)H:r_zz[1,0,10]' \
#        -o FloatMassMu_comb_hires.root 
