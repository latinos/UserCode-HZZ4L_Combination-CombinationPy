rm cards_sm7_2D_base_12fb/ -rf; python makeDCsandWSs_mass.py -i SM_inputs_7TeV -a sm7_2D_base_12fb -b -d 1 -e 0 -u 0 > log2
cd cards_sm7_2D_base_12fb/HCG/126/;
cd -
rm cards_sm8_2D_base_20fb/ -rf; python makeDCsandWSs_mass.py -i SM_inputs_8TeV/ -a sm8_2D_base_20fb -b -d 1 -e 0 -u 0 > log3
cd cards_sm8_2D_base_20fb/HCG/126/;
cp ../../../cards_sm7_2D_base_12fb/HCG/126/h* .
combineCards.py -S hzz4l_7_4mu=hzz4l_4muS_7TeV.txt hzz4l_8_4mu=hzz4l_4muS_8TeV.txt hzz4l_7_4e=hzz4l_4eS_7TeV.txt hzz4l_8_4e=hzz4l_4eS_8TeV.txt hzz4l_7_2e2mu=hzz4l_2e2muS_7TeV.txt hzz4l_8_2e2mu=hzz4l_2e2muS_8TeV.txt  > comb.txt
text2workspace.py comb.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=115,135 -o FloatMass_comb_hzz.root
cd -
#exit


rm cards_sm7_2D_ebe_12fb/ -rf; python makeDCsandWSs_mass.py -i SM_inputs_7TeV -a sm7_2D_ebe_12fb -b -d 0 -e 1 -u 0 > log2
cd cards_sm7_2D_ebe_12fb/HCG/126/;
cd -
rm cards_sm8_2D_ebe_20fb/ -rf; python makeDCsandWSs_mass.py -i SM_inputs_8TeV/ -a sm8_2D_ebe_20fb -b -d 0 -e 1 -u 0 > log3
cd cards_sm8_2D_ebe_20fb/HCG/126/;
cp ../../../cards_sm7_2D_ebe_12fb/HCG/126/h* .
combineCards.py -S hzz4l_7_4mu=hzz4l_4muS_7TeV.txt hzz4l_8_4mu=hzz4l_4muS_8TeV.txt hzz4l_7_4e=hzz4l_4eS_7TeV.txt hzz4l_8_4e=hzz4l_4eS_8TeV.txt hzz4l_7_2e2mu=hzz4l_2e2muS_7TeV.txt hzz4l_8_2e2mu=hzz4l_2e2muS_8TeV.txt  > comb.txt
text2workspace.py comb.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=115,135 -o FloatMass_comb_hzz.root
cd -

rm cards_sm7_3D_base_12fb/ -rf; python makeDCsandWSs_mass.py -i SM_inputs_7TeV -a sm7_3D_base_12fb -b -d 1 -e 1 -u 0 > log2
cd cards_sm7_3D_base_12fb/HCG/126/;
cd -
rm cards_sm8_3D_base_20fb/ -rf; python makeDCsandWSs_mass.py -i SM_inputs_8TeV/ -a sm8_3D_base_20fb -b -d 1 -e 1 -u 0 > log3
cd cards_sm8_3D_base_20fb/HCG/126/;
cp ../../../cards_sm7_3D_base_12fb/HCG/126/h* .
combineCards.py -S hzz4l_7_4mu=hzz4l_4muS_7TeV.txt hzz4l_8_4mu=hzz4l_4muS_8TeV.txt hzz4l_7_4e=hzz4l_4eS_7TeV.txt hzz4l_8_4e=hzz4l_4eS_8TeV.txt hzz4l_7_2e2mu=hzz4l_2e2muS_7TeV.txt hzz4l_8_2e2mu=hzz4l_2e2muS_8TeV.txt  > comb.txt
text2workspace.py comb.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=115,135 -o FloatMass_comb_hzz.root
cd -

#exit

rm cards_sm7_1D_12fb/ -rf; python makeDCsandWSs_mass.py -i SM_inputs_7TeV -a sm7_1D_12fb -b -d 0 -e 0 -u 0 > log2
cd cards_sm7_1D_12fb/HCG/126/;
cd -
rm cards_sm8_1D_20fb/ -rf; python makeDCsandWSs_mass.py -i SM_inputs_8TeV/ -a sm8_1D_20fb -b -d 0 -e 0 -u 0 > log3
cd cards_sm8_1D_20fb/HCG/126/;
cp ../../../cards_sm7_1D_12fb/HCG/126/h* .
combineCards.py -S hzz4l_7_4mu=hzz4l_4muS_7TeV.txt hzz4l_8_4mu=hzz4l_4muS_8TeV.txt hzz4l_7_4e=hzz4l_4eS_7TeV.txt hzz4l_8_4e=hzz4l_4eS_8TeV.txt hzz4l_7_2e2mu=hzz4l_2e2muS_7TeV.txt hzz4l_8_2e2mu=hzz4l_2e2muS_8TeV.txt  > comb.txt
text2workspace.py comb.txt -P HiggsAnalysis.CombinedLimit.PhysicsModel:floatingHiggsMass   --PO higgsMassRange=115,135 -o FloatMass_comb_hzz.root
cd -
exit

