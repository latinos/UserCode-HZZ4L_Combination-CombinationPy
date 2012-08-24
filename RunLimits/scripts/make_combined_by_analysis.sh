#MASSES=$(./massrange 110 145 | grep -v '\.5')
MASSES=$(./massrange 145 600 | grep -v '\.5')
for D in hww hzz hgg hbb htt; do    
    for E in 0 7 8; do 
        for M in $MASSES; do
            #if [[ "$M" != "128" ]]; then continue; fi;
            python scripts/make_combined_card.py $M -e $E -d $D  --zip
            echo "Done $M";
        done
        ~gpetrucc/sh/bann -w 180 "Done $(echo $D | tr '[a-z]' '[A-Z]') $E TeV"
    done
    for M in $MASSES; do
        test -f $M/comb_${D}.txt.gz && bash runbatch.sh make_all_first.sh $M comb7_${D} comb8_${D} comb_${D}
    done
done
