for post in 2D_ebe 3D_base 1D 2D_base 

do
curdir=`pwd` 
dir="cards_sm8_${post}_20fb/HCG"
cp runbatch.sh make_MASS_SCAN.sh lxbatch_runner.sh $dir
cd $dir/126
ln -s FloatMass_comb_hzz.root FloatMass_comb_hzz_exp.root 
cd ..
sh runbatch.sh make_MASS_SCAN.sh --1d --fast 126 comb_hzz.root
for ((first=0;first<10;first++))
do
echo $first
sh runbatch.sh make_MASS_SCAN.sh --1d  -K $first 126 comb_hzz.root -v 1   
#sh runbatch.sh make_MASS_SCAN.sh --2d --less -k2d $first 126 comb_hzz.root  -v 1  --setPhysicsModelParameterRanges MH=120,130:r=0,3

#sh runbatch.sh make_MASS_SCAN.sh --2d --less 126 comb_hzz.root  
#sh runbatch.sh make_MASS_SCAN.sh --1d  -K $first 126 comb_hzz_exp.root --expectSignal=1 --expectSignalMass=126 -t -1 -v 1  
done
cd $curdir
done 

