SCALEFACTORS=$(cat scalefactors.txt)
SEQUENCES=$(cat sequences.txt)
# To toggle local vs dirac job, check line 7 and 11/12 and 16
for seq in ${SEQUENCES[@]}; do
    mkdir -p seq_"$seq"_"$1"
    cd seq_"$seq"_"$1"
    for SFa in ${SCALEFACTORS[@]}; do
        for SFb in ${SCALEFACTORS[@]}; do
            for SFc in ${SCALEFACTORS[@]}; do
                DIRNAME=dpt2_scaled_"$SFa"_"$SFb"_"$SFc"_"$1"
                mkdir -p "$DIRNAME"
                # cp ../SUBMISSIONSCRIPT.sh ./"$DIRNAME"/"$DIRNAME.sh"
                cp ../"$1".inpt ./"$DIRNAME"/"$1".inpt
                python ../diaphite_creator_scaled.py --seq "$seq" --out_file ./"$DIRNAME"/"positions.data" --nx 4 --ny 4 --nz 3 --sfa "$SFa" --sfb "$SFb" --sfc "$SFc" 
                cd "$DIRNAME"
                lmp_serial -in "$1".inpt
                # qsub "$DIRNAME".sh  
                cd ../
            done
        done  
    done
    cd ../
done

