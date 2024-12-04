#!/bin/bash
sampleName=Data2
smode=FALSE
bulks=(cangelosi seqc) 
slength=${#bulks[@]}

job0='job'
echo "#!/bin/bash" > $job0
echo "#SBATCH -c 1" >> $job0
echo "#SBATCH --mem-per-cpu=100G"  >> $job0
echo "#SBATCH -t 2-10:00:00" >> $job0
echo "#SBATCH --output=%x.%j.out" >> $job0
echo "source ~/.bashrc"  >> $job0
echo "module load singularity " >> $job0

for (( i=0; i<${slength}; i++ ));
do
        
        bulk=${bulks[$i]}
        bulk_file=${bulk}_mtx.txt
        jobk=${job0}_${bulk}
        cp $job0 $jobk
        outDir=data/intermediate/RNA/CIBERSORTx/eachGroupDownTo10kCells_${sampleName}_${bulk}
        mkdir -p  ${outDir}
        echo "singularity exec --bind /cm/shared --bind data/intermediate/RNA/CIBERSORTx:/src/data --bind ${outDir}:/src/outdir --bind ${outDir} fractions_latest.sif /src/CIBERSORTxFractions --username yuw1@chop.edu --token your_token_from_cibersortx_website --single_cell TRUE --refsample seurat_rna4Cibersortx_eachGroupDownTo10kCells_${sampleName}.txt --mixture $bulk_file --rmbatchSmode ${smode} --verbose TRUE --perm 100
" >> $jobk
        sbatch $jobk
done

