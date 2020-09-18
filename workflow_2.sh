############ Run the containment screening test
bash contamination_screen.sh

## Download candidate genomes (To be edited later)
mkdir -p /groups/lorolab/Tamer/
cd /groups/lorolab/Tamer/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/205/GCF_000091205.1_ASM9120v1/GCF_000091205.1_ASM9120v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/645/GCF_000149645.2_ASM14964v3/GCF_000149645.2_ASM14964v3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/355/GCF_000151355.1_v2.0/GCF_000151355.1_v2.0_genomic.fna.gz

gunzip *gz

################## kmer containment test
# cDBG for genomes @ k75 and min_ab 1
KMER_SIZE=75
MIN_ABUNDANCE=1
MAX_RAM_MB=50000
THREADS=5
GENOMES_DIR="/groups/lorolab/Tamer"

mkdir genomes_cDBGs && cd genomes_cDBGs

for GENOME in ${GENOMES_DIR}/*;
do
    echo "Processing $GENOME"
    basename "$GENOME"
    GENOME_BASENAME="$(basename -- $GENOME)"
    mkdir -p cDBG_${GENOME_BASENAME};
    cd cDBG_${GENOME_BASENAME};
    OUTPUT_PREFIX=cDBG_k${KMER_SIZE}_${GENOME_BASENAME}
    CMD="bcalm -kmer-size ${KMER_SIZE} -nb-cores ${THREADS} -max-memory ${MAX_RAM_MB} -abundance-min ${MIN_ABUNDANCE} -out ${OUTPUT_PREFIX} -in ${GENOME}"
    clusterize -d -nosub -n ${THREADS} ${CMD} > ${GENOME_BASENAME}_bcalm.qsub
    qsub ${GENOME_BASENAME}_bcalm.qsub
    cd ..
done

cd ..

# -------------------------------------------------------------------------------------------
# Samples cDBG creation @ k75 and min_ab 2

KMER_SIZE=75
MIN_ABUNDANCE=2
MAX_RAM_MB=50000
THREADS=5
SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
SAMPLES_DIR="/groups/lorolab/Astrangia/Astrangia2019"
cDBG_partitioner="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/build/cDBG_partitioner"
IDX_PREFIX="/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/idx_multiSpecies_10_k75.fa"


mkdir -p samples_cDBGs && cd samples_cDBGs

for SAMPLE in $SAMPLES;
do
    mkdir -p ${SAMPLE};
    cd ${SAMPLE};
    OUTPUT_PREFIX=cDBG_k${KMER_SIZE}_${SAMPLE}
    ls -1 ${SAMPLES_DIR}/${SAMPLE}*fastq.gz > reads_${SAMPLE}
    CMD1="bcalm -kmer-size ${KMER_SIZE} -nb-cores ${THREADS} -max-memory ${MAX_RAM_MB} -abundance-min ${MIN_ABUNDANCE} -out ${OUTPUT_PREFIX} -in reads_${SAMPLE}"
    clusterize -d -nosub -n ${THREADS} "${CMD1}" > ${SAMPLE}_bcalm.qsub
    qsub ${SAMPLE}_bcalm.qsub
    rm -rf *h5
    cd ..
done

# -------------------------------------------------------------------------------------------
# kmer counting 

cd /groups/lorolab/mr-eyes/oveview_exp
SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
SAMPLES_DIR="/groups/lorolab/mr-eyes/final_experiment/samples_cDBGs"
OUTPUT_DIR="/groups/lorolab/mr-eyes/oveview_exp/samples_kfs"
kmerCount="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/kmerCount.py"


for SAMPLE in $SAMPLES;
do 
    CMD="python $kmerCount $SAMPLES_DIR/$SAMPLE/cDBG_k75_$SAMPLE.unitigs.fa $OUTPUT_DIR";
    clusterize -d -n 2 "${CMD}"; 
done


GENOMES_KFS=/groups/lorolab/mr-eyes/oveview_exp/genomes_kf
GENOMES_DIR=/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/*unitigs.fasta
for GENOME in $GENOMES;
do 
    CMD="python $kmerCount $GENOME $GENOMES_KFS";
    clusterize -d -n 2 "${CMD}"; 
done


# --------------------------------------------------
# Contamnination study

SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
OUTPUT_DIR="/groups/lorolab/mr-eyes/oveview_exp/samples_kfs"

commonKmers="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/build/commonKmers"
GENOMES_KFS=/groups/lorolab/mr-eyes/oveview_exp/genomes_kf/*mqf

OUT_TSV="contamination_report.tsv"

touch $OUT_TSV
echo -e "ref\tsample\tcommon_kmers\tsamples_kmers\tcontainment_percentage" >> $OUT_TSV

for GENOME in $GENOMES_KFS;
do  
    for SAMPLE in $SAMPLES;
    do 
        CMD="$commonKmers $GENOME $OUTPUT_DIR/idx_cDBG_k75_$SAMPLE.unitigs.fa.mqf >> $OUT_TSV"
        clusterize -d -n 1 "${CMD}"; 
    done 
done

touch detailed_Contamination_report.tsv
echo -e "ref\tsample\tcommon_kmers\tsamples_kmers\tcontainment_percentage" >> detailed_Contamination_report.tsv
cat contamStats*tsv >> detailed_Contamination_report.tsv
python /groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/transform.py detailed_Contamination_report.tsv


## ----------- genomes_cDBG Indexing ----------------------
# Move genomes of interest into /groups/lorolab/Tamer/10genomes
mv *gz 10genomes/
cd 10genomes/ && gunzip *gz


cd genomes_cDBGs
python /groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/unitigsTokProcessorFormat.py multiSpecies_9_k75 10genomes/*fa
clusterize -d -nosub -n 8 python /groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/indexing.py multiSpecies_9_k75.fa multiSpecies_9_k75.fa.names 21 > indexing_multiSpecies.qsub
qsub indexing_multiSpecies.qsub
IDX_PREFIX="/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/idx_multiSpecies_10_k75.fa"

# -------------------------------------------------------------------------------------------


# cDBG Partitioning
SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
SAMPLES_DIR="/groups/lorolab/mr-eyes/final_experiment/samples_cDBGs"
cDBG_partitioner="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/build/cDBG_partitioner"
IDX_PREFIX="/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/idx_multiSpecies_10_k75.fa"
KMER_SIZE=75

for SAMPLE in $SAMPLES;
do
    cd ${SAMPLE};
    OUTPUT_PREFIX=cDBG_k${KMER_SIZE}_${SAMPLE}
    CMD="/usr/bin/time -v ${cDBG_partitioner} ${OUTPUT_PREFIX}.unitigs.fa ${IDX_PREFIX}"
    clusterize -d -nosub -n 1 "${CMD}" > ${SAMPLE}_partitioning.qsub
    qsub ${SAMPLE}_partitioning.qsub
    cd ..
done

cd ..

# -------------------------------------------------------------------------------------------

# Combined index (Genomes index + Samples partitions)

SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"

declare -A groupNames
groupNames[1]=cDBG_k75_apoculataassemblyscaffolds_chromosome_level
groupNames[2]=cDBG_k75_con_Cfortranscriptome
groupNames[3]=cDBG_k75_con_Scargenome_assembly
groupNames[4]=cDBG_k75_GCA_0005073051_ASM50730v1_genomic
groupNames[5]=cDBG_k75_GCA_0005120851_Reti_assembly10_genomic
groupNames[6]=cDBG_k75_GCA_0017500851_Fracy1_genomic
groupNames[7]=cDBG_k75_GCA_0020494552_P_umbilicalis_v1_genomic
groupNames[8]=cDBG_k75_GCF_0001435352_ASM14353v4_genomic
groupNames[9]=cDBG_k75_GCF_0003727251_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v10_genomic
groupNames[10]=cDBG_k75_Ost699_genome_LATEST

contigsFasta=$(pwd)/allSamples_contigs.fa
contigsNames=$(pwd)/allSamples_contigs.fa.names

touch ${contigsFasta}
touch ${contigsNames}

CONTIGS_COUNTER=1

for SAMPLE in $SAMPLES;
  do
    echo "Processing $SAMPLE"
    for GENOME_ID in 1 2 3 4 5 6 7 8 9 10;
      do
          ((CONTIGS_COUNTER++))
          echo "Processing Genome ${GENOME_ID}"
          originalCDBG=${groupNames[$GENOME_ID]}
          cat ${SAMPLE}/genome_${GENOME_ID}_partition.fa | awk -v gid=$GENOME_ID -v sample=$SAMPLE 'BEGIN{OFS="\n";}!/^>/{print ">"gid"."sample"."NR/2,$0}' >> ${contigsFasta}
          cat ${SAMPLE}/genome_${GENOME_ID}_partition.fa | awk -v gid=$GENOME_ID -v sample=$SAMPLE -v seqName=$originalCDBG 'BEGIN{OFS="\t";}!/^>/{print ">"gid"."sample"."NR/2,seqName}' >> ${contigsNames}
    done;
done;

# Merging the contigs with the original cDBGs


multiSpecisCDBG_fasta=/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/multiSpecies_10_k75.fa
multiSpecisCDBG_names=/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/multiSpecies_10_k75.fa.names

cat ${contigsFasta} ${multiSpecisCDBG_fasta} > allSamples_with_10genomes.fa
cat ${contigsNames} ${multiSpecisCDBG_names} > allSamples_with_10genomes.fa.names


#------------------------------------------------------------------------
# Combined indexing

INDEXING=/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/indexing.py

clusterize -d -n 6 /usr/bin/time -v python ${INDEXING} allSamples_with_10genomes.fa allSamples_with_10genomes.fa.names 21

FULL_IDX_PREFIX=/groups/lorolab/mr-eyes/final_experiment/idx_allSamples_with_10genomes.fa


# ----------------------------------------------------------------------------------

# Reads partitioning

SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
SAMPLES_DIR="/groups/lorolab/Astrangia/Astrangia2019"
fastqPartitioner="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/build/fastq_partitioner"
FULL_IDX_PREFIX=/groups/lorolab/mr-eyes/final_experiment/idx_allSamples_with_10genomes.fa


mkdir reads_partitions && cd reads_partitions

for SAMPLE in $SAMPLES;
do
    mkdir ${SAMPLE} && cd ${SAMPLE}
    R1=${SAMPLES_DIR}/${SAMPLE}_R1_001.fastq.gz;
    R2=${SAMPLES_DIR}/${SAMPLE}_R2_001.fastq.gz;
    CMD="${fastqPartitioner} ${FULL_IDX_PREFIX} ${R1} ${R2}"
    clusterize -d -nosub -n 1 ${CMD} > ${SAMPLE}_partitioning.qsub
    qsub ${SAMPLE}_partitioning.qsub
    cd ..
done

cd ..

# -------------------------------------------------------------------------------------

# merging samples (STILL DEVELOPMENT)

genome_${GENOME_ID}_readsPartition_R1.fastq
genome_${GENOME_ID}_readsPartition_R2.fastq

for GENOME_ID in 1 2 3 4 5 6 7 8 9 10;
do
    touch decontaminated_genome_${GENOME_ID}_R1.fastq
    touch decontaminated_genome_${GENOME_ID}_R2.fastq
done

for SAMPLE in $SAMPLES;
  do
    echo "Processing $SAMPLE"
    for GENOME_ID in 1 2 3 4 5 6 7 8 9 10;
      do
          echo "Processing Genome ${GENOME_ID}"
          originalCDBG=${groupNames[$GENOME_ID]}
          cat ${SAMPLE}/genome_${GENOME_ID}_readsPartition_R1.fastq | awk 'BEGIN{OFS="\n";}!/^>/{print ">"NR/2,$0}' >> decontaminated_genome_${GENOME_ID}_R1.fastq
          cat ${SAMPLE}/genome_${GENOME_ID}_readsPartition_R1.fastq | awk 'BEGIN{OFS="\n";}!/^>/{print ">"NR/2,$0}' >> decontaminated_genome_${GENOME_ID}_R1.fastq
    done;
done;


# -----------------------------------------------

# Generate lengths summary stats TSV

SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
SAMPLES_ORIGINAL_DIR="/groups/lorolab/Astrangia/Astrangia2019"
SAMPLES_OUTPUT_DIR="/groups/lorolab/mr-eyes/final_experiment/reads_partitions"

SUMMARY_TSV="partitioning_summary.tsv"
touch ${SUMMARY_TSV}
printf "sample\ttotal" >> ${SUMMARY_TSV}

for i in 1 2 3 4 5 6 7 8 9 10
do  
    printf "\tgenome_${i}" >> ${SUMMARY_TSV}
done

printf "\tunmapped\n" >> ${SUMMARY_TSV}


for SAMPLE in $SAMPLES
do 
    TSV_LINE="${SAMPLE}\t"
    echo "Processing ${SAMPLE}"
    ORIGINAL_GENOME_LINES=$(zcat ${SAMPLES_ORIGINAL_DIR}/${SAMPLE}_R1_001.fastq.gz | wc -l)
    GENOMES_SEQS=$((ORIGINAL_GENOME_LINES / 4))
    TSV_LINE+="${GENOMES_SEQS}"
    for GENOME_ID in 1 2 3 4 5 6 7 8 9 10
    do      
            echo "Processing ${SAMPLE}/${GENOME_ID}"
            PARTITION_LINES=$(cat ${SAMPLES_OUTPUT_DIR}/${SAMPLE}/genome_${GENOME_ID}_readsPartition_R1.fastq | wc -l)
            PARTITION_SEQS=$((PARTITION_LINES / 4))
            TSV_LINE+="\t${PARTITION_SEQS}"
    done
    UNMAPPED_LINES=$(cat ${SAMPLES_OUTPUT_DIR}/${SAMPLE}/unmapped_partition.fa_R1.fastq | wc -l)
    UNMAPPED_SEQS=$((UNMAPPED_LINES / 4))
    printf "${TSV_LINE}\t${UNMAPPED_SEQS}\n" >> ${SUMMARY_TSV}
done



# ------------------------------------------------------------------------------

# TRIMMING

adap="$CONDA_PREFIX/share/trimmomatic-0.39-1/adapters"
SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
SAMPLES_OUTPUT_DIR="/groups/lorolab/mr-eyes/final_experiment/reads_partitions"

mkdir -p trimmed_reads && cd trimmed_reads

for SAMPLE in $SAMPLES;
    do
        mkdir -p $SAMPLE && cd $SAMPLE
        echo "Processing $SAMPLE"
        for GENOME_ID in 1 2 3 4 5 6 7 8 9 10;
        do
            R1=${SAMPLES_OUTPUT_DIR}/${SAMPLE}/genome_${GENOME_ID}_readsPartition_R1.fastq
            R2=${SAMPLES_OUTPUT_DIR}/${SAMPLE}/genome_${GENOME_ID}_readsPartition_R2.fastq

            OP_R1_SE=trimmed_genome_${GENOME_ID}_readsPartition_R1_SE.fastq
            OP_R2_SE=trimmed_genome_${GENOME_ID}_readsPartition_R2_SE.fastq
            OP_R1_PE=trimmed_genome_${GENOME_ID}_readsPartition_R1_PE.fastq
            OP_R2_PE=trimmed_genome_${GENOME_ID}_readsPartition_R2_PE.fastq

            CMD="trimmomatic PE -threads 4 -phred33 ${R1} ${R2} ${OP_R1_PE} ${OP_R1_SE} ${OP_R2_PE} ${OP_R2_SE} ILLUMINACLIP:${adap}/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:25"
            clusterize -d -nosub -n 4 ${CMD} > ${SAMPLE}_trimmomatic.qsub
            qsub ${SAMPLE}_trimmomatic.qsub
        done;
        cd ..
done;

# -----------------------------------------------------------
# Compressing fastq files to save disk space


SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
TRIMMED_SAMPLES_DIR="/groups/lorolab/mr-eyes/final_experiment/trimmed_reads"

for SAMPLE in $SAMPLES;
do
    echo "Gzipping $SAMPLE";
    clusterize -d -n 1 gzip "$SAMPLE/*fastq"
done


# -----------------------------------------------------------
# Assembly (Trinity)

SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
TRIMMED_SAMPLES_DIR="/groups/lorolab/mr-eyes/final_experiment/trimmed_reads"
# MERGED_SE_FILES=${TRIMMED_SAMPLES_DIR}/all_trimmed_SE.fastq
THREADS=32

TRINITY_CMD="Trinity --seqType fq --CPU 32--max_memory 400G --left  --right "


for SAMPLE in $SAMPLES;
do
    cd $SAMPLE
    
    for GENOME_ID in 1; # Just the apoculata
    do
        
    done
    cd ..
done
