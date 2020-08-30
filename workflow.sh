IDX_PREFIX="/groups/lorolab/mr-eyes/partitioning_cDBG_exp/idx_withContigs_MiltiSpecies.fa"
#fastqPartitioner="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/build/fastq_partitioner"

CURRENT=$(pwd)


# -------------------------------------------------------------------------------------------


# cDBG of all reads at k=75
cDBG_KMER_SIZE=75
MIN_ABUNDANCE=1
MAX_RAM_MB=60000
THREADS=16
SAMPLES_DIR="/groups/lorolab/Astrangia/Astrangia2019"
SAMPLES_CDBG_DIR=samples_cDBG_k75
mkdir ${SAMPLES_CDBG_DIR};
cd ${SAMPLES_CDBG_DIR};
OUTPUT_PREFIX=cDBG_k${cDBG_KMER_SIZE}_samples
ls -1 ${SAMPLES_DIR}/*fastq.gz > reads_paths
CMD="/usr/bin/time -v bcalm -kmer-size ${cDBG_KMER_SIZE} -nb-cores ${THREADS} -max-memory ${MAX_RAM_MB} -abundance-min ${MIN_ABUNDANCE} -out ${OUTPUT_PREFIX} -in reads_paths"
clusterize -d -nosub -n ${THREADS} ${CMD} > bcalm_allSamples.qsub
qsub bcalm_allSamples.qsub
cd ..

# -------------------------------------------------------------------------------------------

# cDBG Partitioning
cDBG_partitioner="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/build/cDBG_partitioner"
mkdir contigs_partitions && cd contigs_partitions
/usr/bin/time -v ${cDBG_partitioner} ${CURRENT}/${SAMPLES_CDBG_DIR}/${OUTPUT_PREFIX}.unitigs.fa ${IDX_PREFIX}

# -------------------------------------------------------------------------------------------

# Indexing

declare -A groupNames
groupNames[1]=cDBG_k21_apoculataassemblyscaffolds_chromosome_level
groupNames[1]=cDBG_k21_Cadsp1_AssemblyScaffolds
groupNames[3]=cDBG_k21_con_Cfortranscriptome
groupNames[4]=cDBG_k21_con_Scargenome_assembly
groupNames[5]=cDBG_k21_GCA_0005073051_ASM50730v1_genomic
groupNames[6]=cDBG_k21_GCA_0005120851_Reti_assembly10_genomic
groupNames[7]=cDBG_k21_GCF_0003727251_Emiliana_huxleyi_CCMP1516_main_genome_assembly_v10_genomic


contigsFasta=allSamples_contigs.fa
contigsNames=allSamples_contigs.fa.names

touch ${contigsFasta}
touch ${contigsNames}

for GENOME_ID in 1 2 3 4 5 6 7;
do
    originalCDBG=${groupNames[$GENOME_ID]}
    cat contigs_partitions/genome_${GENOME_ID}_partition.fa >> ${contigsFasta}
    grep ">" contigs_partitions/genome_${GENOME_ID}_partition.fa | cut -c2- | awk -F' ' -v groupname=${originalCDBG} '{print $0"\t"groupname}' >> ${contigsNames};
done;

# Merging the contigs with the original cDBGs

multiSpecisCDBG_fasta=/groups/lorolab/mr-eyes/oveview_exp/multiSpecies_indexing/multiSpecies7.fa
multiSpecisCDBG_names=/groups/lorolab/mr-eyes/oveview_exp/multiSpecies_indexing/multiSpecies7.fa.names

cat ${contigsFasta} ${multiSpecisCDBG_fasta} > allSamplesContigs_with_multiSpeciesGenomes.fa
cat ${contigsNames} ${multiSpecisCDBG_names} > allSamplesContigs_with_multiSpeciesGenomes.fa.names

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