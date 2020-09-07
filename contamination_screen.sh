# UNITE (https://unite.ut.ee/) is a database centered on the eukaryotic nuclear ribosomal ITS region. All eukaryotic ITS sequences from the International Nucleotide Sequence Database Collaboration (INSDC) are clustered to approximately the species level (97-100% similarity in steps of 0.5%), and all such species hypotheses are given a DOI to facilitate unambiguous scientific communication and data assembly.
# We are using the ITS sequences of the UNITE database to identify possible contaminants in RNAseq sample.
# Currently we are using v8.2 (2020-02-04) for All eukaryotes (12666 RefS + 91074 RepS)
# When using this resource, please cite it as follows: Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2020): UNITE general FASTA release for eukaryotes 2. Version 04.02.2020. UNITE Community. https://doi.org/10.15156/BIO/786371


work_dir=$(pwd)
mkdir -p its_iden && cd its_iden
wget https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz
tar xvzf 1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz

## blast blast+/2.7.1
trinity_tran=$work_dir/Ast2019_allAssembly_trinity.fasta
its_ref=$work_dir/its_iden/sh_general_release_s_all_04.02.2020/sh_general_release_dynamic_s_all_04.02.2020.fasta
makeblastdb -in $its_ref -input_type fasta -dbtype nucl

blastn -query $trinity_tran \
-db $its_ref \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
-dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
-max_target_seqs 10  -out Ast2019_allEuComp ## -perc_identity 90 -qcov_hsp_perc 50

sort -k1,1 -k11,11g -k5,5nr Ast2019_allEuComp | sort -u -k1,1 --merge | sort -k2,2 > Ast2019_allEuComp_best

echo "qseqid ID GenBank_accession SH_accession curation kingdom phylum class order family genus species pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" | tr " " "\t" > Ast2019_allEuComp_best_ext
cat Ast2019_allEuComp_best | tr "|" "\t" | tr ";" "\t" | sort -k6,6 -k11,11 -k12,12 >> Ast2019_allEuComp_best_ext
tail -n+2 Ast2019_allEuComp_best_ext | awk -F"\t" '{print $6}' | sort | uniq -c | sort -k1,1nr > kingdom.stats
tail -n+2 Ast2019_allEuComp_best_ext | awk -F"\t" '{print $7}' | sort | uniq -c | sort -k1,1nr > phylum.stats
tail -n+2 Ast2019_allEuComp_best_ext | awk -F"\t" '{print $9}' | sort | uniq -c | sort -k1,1nr > order.stats
tail -n+2 Ast2019_allEuComp_best_ext | awk -F"\t" '{print $11}' | sort | uniq -c | sort -k1,1nr > genus.stats
tail -n+2 Ast2019_allEuComp_best_ext | awk -F"\t" '{print $12}' | sort | uniq -c | sort -k1,1nr > species.stats

head -n1 Ast2019_allEuComp_best_ext > Ast2019_allEuComp_best_ext_best
tail -n+2  Ast2019_allEuComp_best_ext | awk -F"\t" '{if($13>95 && $14>200)print}' >> Ast2019_allEuComp_best_ext_best
tail -n+2  Ast2019_allEuComp_best_ext | awk -F"\t" '{if($13>90 && $14>100)print}' | grep -i p__Cercozoa  >> Ast2019_allEuComp_best_ext_best

## edit the file to choose the best enteries and save in Ast2019_allEuComp_best_ext_best_manual
cat Ast2019_allEuComp_best_ext_best_manual | tr '\t' ';' > Ast2019_allEuComp_best_ext_best_manual.csv  ## 19 enteries plus header line 

#-------
#-------
#### Test possible contaminants in published transcriptomes
## 1. Astrangia_poculata_transcriptome
mkdir -p refs
wget http://sites.bu.edu/davieslab/files/2020/05/Astrangia_poculata_transcriptome_2020v1.zip
unzip Astrangia_poculata_transcriptome_2020v1.zip
mv Astrangia_poculata_transcriptome_2020v1/ast_MPCC2017_transcriptome.fasta refs/coral_Apoc.transcriptome.fasta
rm -rf __MACOSX  Astrangia_poculata_transcriptome_2020v1

coral_Apoc=$work_dir/candidates/coral_Apoc.transcriptome.fasta
blastn -query $coral_Apoc \
-db $its_ref \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
-dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
-max_target_seqs 10  -out coral_Apoc_allEuComp ## -perc_identity 90 -qcov_hsp_perc 50

sort -k1,1 -k11,11g -k5,5nr coral_Apoc_allEuComp | sort -u -k1,1 --merge | sort -k2,2 > coral_Apoc_allEuComp_best
grep -v p__Cnidaria coral_Apoc_allEuComp_best > coral_Apoc_allEuComp_contaminants
#-------
## 2. coral Orbicella faveolata (Scleractinia-Merulinidae) 
wget https://dfzljdn9uc3pi.cloudfront.net/2016/1616/1/SI_4a_coral_assembly.fasta.zip
unzip SI_4a_coral_assembly.fasta.zip
mv "SI 2 coral assembly.fasta" refs/coral_Ofav.transcriptome.fasta
rm -rf __MACOSX SI_4a_coral_assembly.fasta.zip

coral_Ofav=$work_dir/Tamer/coral_Ofav.transcriptome.fasta
blastn -query $coral_Ofav \
-db $its_ref \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
-dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
-max_target_seqs 10  -out coral_Ofav_allEuComp ## -perc_identity 90 -qcov_hsp_perc 50

sort -k1,1 -k11,11g -k5,5nr coral_Ofav_allEuComp | sort -u -k1,1 --merge | sort -k2,2 > coral_Ofav_allEuComp_best
grep -v p__Cnidaria coral_Ofav_allEuComp_best > coral_Ofav_allEuComp_contaminants
#-------
## 3. ACROPORA TENUIS
wget -O aten_july2014.zip https://www.dropbox.com/s/7sv415goixnnhvy/aten_july2014.zip?dl=0
unzip aten_july2014.zip
mv aten_july2014/aten.fasta refs/coral_Aten.transcriptome.fasta
rm -rf __MACOSX aten_july2014 aten_july2014.zip

coral_Aten=$work_dir/Tamer/coral_Aten.transcriptome.fasta
blastn -query $coral_Aten \
-db $its_ref \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
-dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
-max_target_seqs 10  -out coral_Aten_allEuComp ## -perc_identity 90 -qcov_hsp_perc 50

sort -k1,1 -k11,11g -k5,5nr coral_Aten_allEuComp | sort -u -k1,1 --merge | sort -k2,2 > coral_Aten_allEuComp_best
grep -v p__Cnidaria coral_Aten_allEuComp_best > coral_Aten_allEuComp_contaminants
#-------
## 4. MONTIPORA AEQUITUBERCULATA
wget -O maequituberculata_transcriptome_july2014.zip https://www.dropbox.com/s/mnkqwz7oi3l226c/maequituberculata_transcriptome_july2014.zip?dl=0
unzip maequituberculata_transcriptome_july2014.zip -d maequituberculata_transcriptome
mv maequituberculata_transcriptome/maeq_coral.fasta refs/coral_Maeq.transcriptome.fasta
rm -rf maequituberculata_transcriptome maequituberculata_transcriptome_july2014.zip

coral_Maeq=$work_dir/Tamer/coral_Maeq.transcriptome.fasta
blastn -query $coral_Maeq \
-db $its_ref \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
-dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
-max_target_seqs 10  -out coral_Maeq_allEuComp ## -perc_identity 90 -qcov_hsp_perc 50

sort -k1,1 -k11,11g -k5,5nr coral_Maeq_allEuComp | sort -u -k1,1 --merge | sort -k2,2 > coral_Maeq_allEuComp_best
grep -v p__Cnidaria coral_Maeq_allEuComp_best > coral_Maeq_allEuComp_contaminants
#-------
## 5. PORITES ASTREOIDES
wget -O pastreoides_2014.zip https://www.dropbox.com/s/dzvr219crsbgt0n/pastreoides_2014.zip?dl=0
unzip pastreoides_2014.zip
mv pastreoides_2014/pastreoides_may2014/past.fasta refs/coral_Past.transcriptome.fasta
rm -rf __MACOSX pastreoides_2014 pastreoides_2014.zip

coral_Past=$work_dir/Tamer/coral_Past.transcriptome.fasta
blastn -query $coral_Past \
-db $its_ref \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qcovs qcovhsp" \
-dust 'yes' -best_hit_overhang 0.25 -best_hit_score_edge 0.25 \
-max_target_seqs 10  -out coral_Past_allEuComp ## -perc_identity 90 -qcov_hsp_perc 50

sort -k1,1 -k11,11g -k5,5nr coral_Past_allEuComp | sort -u -k1,1 --merge | sort -k2,2 > coral_Past_allEuComp_best
grep -v p__Cnidaria coral_Past_allEuComp_best > coral_Past_allEuComp_contaminants
#-------




