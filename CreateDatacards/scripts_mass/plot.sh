


for post in 2D_base 1D   2D_ebe 3D_base  
do 


	echo  $post
	#post="2D_ebe"
	curdir=`pwd` 
	dir=cards_sm8_${post}_20fb/HCG/126
	cd $dir

  #hadd -f -k moriond_${post}_PLP_126.root higgsCombine.TXT_PLP.ProfileLikelihood.mH1* 
  #hadd -f -k moriond_${post}_PLPBE_126.root higgsCombine.TXT_PLPBE.ProfileLikelihood.mH1* 
  #cp moriond_* ~/www/posthcp
  

  	cp higgsCombineHZZ_MASS_SCAN_1D_FAST.MultiDimFit.mH126.root ~/www/posthcp/higgsCombineHZZ_MASS_SCAN_1D-${post}-moriond_obs_fast.root 
	echo $dir
	hadd -f higgsCombineHZZ_MASS_SCAN_1D.root higgsCombineHZZ_MASS_SCAN_1D.*Mul* 
	cp higgsCombineHZZ_MASS_SCAN_1D.root ~/www/posthcp/higgsCombineHZZ_MASS_SCAN_1D-${post}-moriond_obs.root 
	#echo hadd 1 

	#hadd -f higgsCombineHZZ_MASS_SCAN_2D.root higgsCombineHZZ_MASS_SCAN_2D.*Mul* 
	#cp higgsCombineHZZ_MASS_SCAN_2D.root ~/work/cmshcg/trunk/hcp2012/results/ 
	#cd ~/work/cmshcg/trunk/hcp2012/
	#root -l -b -q makeBands.cxx
	#cp results/bands.root ~/www/posthcp/bands_mass_${post}_moriond_obs.root 

	cd $curdir
done
