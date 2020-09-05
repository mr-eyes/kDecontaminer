fastqPartitioner="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/build/fastq_partitioner"

CURRENT=$(pwd)


# cDBG for genomes @ k75


KMER_SIZE=75
MIN_ABUNDANCE=2
MAX_RAM_MB=50000
THREADS=5
GENOMES_DIR="/groups/lorolab/Tamer/7genomes"

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

# ----------- genomes_cDBG Indexing ----------------------

cd genomes_cDBGs
python /groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/unitigsTokProcessorFormat.py multiSpecies_7_k75 */*fa
clusterize -d -nosub -n 1 python /groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/indexing.py multiSpecies_7_k75.fa multiSpecies_7_k75.fa.names 21 > indexing_multiSpecies.qsub
qsub indexing_multiSpecies.qsub
IDX_PREFIX="/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/idx_multiSpecies_7_k75.fa"

# -------------------------------------------------------------------------------------------

# Samples cDBG creation & partitioning

KMER_SIZE=75
MIN_ABUNDANCE=2
MAX_RAM_MB=50000
THREADS=5
SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
SAMPLES_DIR="/groups/lorolab/Astrangia/Astrangia2019"
cDBG_partitioner="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/build/cDBG_partitioner"
IDX_PREFIX="/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/idx_multiSpecies_7_k75.fa"


mkdir -p samples_cDBGs && cd samples_cDBGs

for SAMPLE in $SAMPLES;
do
    mkdir -p ${SAMPLE};
    cd ${SAMPLE};
    OUTPUT_PREFIX=cDBG_k${KMER_SIZE}_${SAMPLE}
    ls -1 ${SAMPLES_DIR}/${SAMPLE}*fastq.gz > reads_${SAMPLE}
    CMD1="bcalm -kmer-size ${KMER_SIZE} -nb-cores ${THREADS} -max-memory ${MAX_RAM_MB} -abundance-min ${MIN_ABUNDANCE} -out ${OUTPUT_PREFIX} -in reads_${SAMPLE}"
    CMD2="/usr/bin/time -v ${cDBG_partitioner} ${OUTPUT_PREFIX}.unitigs.fa ${IDX_PREFIX}"
    clusterize -d -nosub -n ${THREADS} "${CMD1} && ${CMD2}" > ${SAMPLE}_bcalm.qsub
    qsub ${SAMPLE}_bcalm.qsub
    rm -rf *h5
    cd ..
done

cd ..


# -------------------------------------------------------------------------------------------

# Combined index (Genomes index + Samples partitions)

SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"


declare -A groupNames
groupNames[1]=cDBG_k75_apoculataassemblyscaffolds_chromosome_level
groupNames[2]=cDBG_k75_Cadsp1_AssemblyScaffolds
groupNames[3]=cDBG_k75_con_Cfortranscriptome
groupNames[4]=cDBG_k75_con_Scargenome_assembly
groupNames[5]=cDBG_k75_GCA_0003970851_Porphyridium_purpureum_genomic
groupNames[6]=cDBG_k75_GCA_0005073051_ASM50730v1_genomic
groupNames[7]=cDBG_k75_GCA_0005120851_Reti_assembly10_genomic


contigsFasta=$(pwd)/allSamples_contigs.fa
contigsNames=$(pwd)/allSamples_contigs.fa.names

touch ${contigsFasta}
touch ${contigsNames}

CONTIGS_COUNTER=1

for SAMPLE in $SAMPLES;
  do
    echo "Processing $SAMPLE"
    for GENOME_ID in 1 2 3 4 5 6 7;
      do
          echo "Processing Genome ${GENOME_ID}"
          originalCDBG=${groupNames[$GENOME_ID]}

          grep -v "^>" ${SAMPLE}/genome_${GENOME_ID}_partition.fa | while read -r seq;
            do
              echo -e ">${CONTIGS_COUNTER}\n$seq" >> ${contigsFasta}
              echo -e "$CONTIGS_COUNTER\t${originalCDBG}" >> ${contigsNames};
              ((CONTIGS_COUNTER=CONTIGS_COUNTER+1))
          done;
    done;
done;

# Merging the contigs with the original cDBGs

multiSpecisCDBG_fasta=/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/multiSpecies_7_k75.fa
multiSpecisCDBG_names=/groups/lorolab/mr-eyes/final_experiment/genomes_cDBGs/multiSpecies_7_k75.fa.names

cat ${contigsFasta} ${multiSpecisCDBG_fasta} > allSamples_with_7Genomes.fa
cat ${contigsNames} ${multiSpecisCDBG_names} > allSamples_with_7Genomes.fa.names

INDEXING=/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/indexing.py


/usr/time/bin -v python ${INDEXING} allSamplesContigs_with_multiSpeciesGenomes.fa allSamplesContigs_with_multiSpeciesGenomes.fa.names 21


# Reads partitioning

SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"

mkdir reads_final_partitions && cd reads_final_partitions

for SAMPLE in $SAMPLES;
do
    mkdir ${SAMPLE};
    cd ${SAMPLE};
    R1=${SAMPLES_DIR}/${SAMPLE}_R1_001.fastq.gz;
    R2=${SAMPLES_DIR}/${SAMPLE}_R2_001.fastq.gz;
    CMD="${fastqPartitioner} ${IDX_PREFIX} ${R1} ${R2}"
    clusterize -d -nosub -n 1 ${CMD} > ${SAMPLE}_partitioning.qsub
    qsub ${SAMPLE}_partitioning.qsub
    cd ..
done

cd ..