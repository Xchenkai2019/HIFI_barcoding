#1. extract high quality ccs reads
ccs subreads.bam subreads.ccs.bam --min-passes 3 --min-rq 0.99   ## with 99% accuracy and 3 subreads required to generate CCS for a ZMW
bam2fastx -a -A -o subreads.ccs.bam.fa subreads.ccs.bam ## convert bam format to fasta format

#2. use hmmscan to alignment ccs reads to Rfam rDNA database and extract rDNA sequences
hmmscan all.rDNA.hmm subreads.ccs.bam.fa > subreads.ccs.bam.fa.domain.nhmmer
more subreads.ccs.bam.fa.domain.nhmmer|grep 'SSU'|perl -ne 'if ($!~/^\#/){print "$";}'|sort -k1,1 -k5,5g -k4,4nr > subreads.ccs.bam.fa.domain.nhmmer.SSU
more subreads.ccs.bam.fa.domain.nhmmer|grep 'LSU'|perl -ne 'if ($!~/^\#/){print "$";}'|sort -k1,1 -k5,5g -k4,4nr > subreads.ccs.bam.fa.domain.nhmmer.LSU
more subreads.ccs.bam.fa.domain.nhmmer|grep '5.8S_rRNA_euk|5S_rRNA'|perl -ne 'if ($!~/^\#/){print "$";}'|sort -k1,1 -k5,5g -k4,4nr > subreads.ccs.bam.fa.domain.nhmmer.5.8S
perl get.m8.first.line.pl subreads.ccs.bam.fa.domain.nhmmer.SSU > subreads.ccs.bam.fa.domain.nhmmer.SSU.firstline
perl get.m8.first.line.pl subreads.ccs.bam.fa.domain.nhmmer.LSU > subreads.ccs.bam.fa.domain.nhmmer.LSU.firstline
perl get.m8.first.line.pl subreads.ccs.bam.fa.domain.nhmmer.5.8S > subreads.ccs.bam.fa.domain.nhmmer.5.8S.firstline
more subreads.ccs.bam.fa.domain.nhmmer.SSU.firstline|awk '$9=="+"' |awk '{print $1"\t"$12"\t"$13}' >subreads.ccs.bam.fa.domain.nhmmer.SSU.firstline.plusregion
more subreads.ccs.bam.fa.domain.nhmmer.SSU.firstline|awk '$9=="-"' |awk '{print $1"\t"$13"\t"$12}' > subreads.ccs.bam.fa.domain.nhmmer.SSU.firstline.minusregion
more subreads.ccs.bam.fa.domain.nhmmer.LSU.firstline|awk '$9=="+"' |awk '{print $1"\t"$12"\t"$13}' >subreads.ccs.bam.fa.domain.nhmmer.LSU.firstline.plusregion
more subreads.ccs.bam.fa.domain.nhmmer.LSU.firstline|awk '$9=="-"' |awk '{print $1"\t"$13"\t"$12}' > subreads.ccs.bam.fa.domain.nhmmer.LSU.firstline.minusregion
more subreads.ccs.bam.fa.domain.nhmmer.5.8S.firstline|awk '$9=="+"' |awk '{print $1"\t"$12"\t"$13}' >subreads.ccs.bam.fa.domain.nhmmer.5.8S.firstline.plusregion
more subreads.ccs.bam.fa.domain.nhmmer.5.8S.firstline|awk '$9=="-"' |awk '{print $1"\t"$13"\t"$12}' > subreads.ccs.bam.fa.domain.nhmmer.5.8S.firstline.minusregion
seqkit subseq subreads.ccs.bam.fa subreads.ccs.bam.fa.domain.nhmmer.SSU.firstline.plusregion >SSU.fasta
seqkit subseq subreads.ccs.bam.fa subreads.ccs.bam.fa.domain.nhmmer.SSU.firstline.minusregion|perl reverse_seq.pl - >>SSU.fasta
seqkit subseq subreads.ccs.bam.fa subreads.ccs.bam.fa.domain.nhmmer.LSU.firstline.plusregion >LSU.fasta
seqkit subseq subreads.ccs.bam.fa subreads.ccs.bam.fa.domain.nhmmer.LSU.firstline.minusregion|perl reverse_seq.pl - >>LSU.fasta
seqkit subseq subreads.ccs.bam.fa subreads.ccs.bam.fa.domain.nhmmer.5.8S.firstline.plusregion >5.8S.fasta
seqkit subseq subreads.ccs.bam.fa subreads.ccs.bam.fa.domain.nhmmer.5.8S.firstline.minusregion|perl reverse_seq.pl - >>5.8S.fasta

#3. use hmmscan to alignment ccs reads to CO1 HMM (PF00115) and extract CO1 sequences
hmmscan PF00115.hmm subreads.ccs.bam.fa > CO1.domain.nhmmer


#4. taxonomy assignments
vsearch --cluster_size SSU.fasta --strand both --id 0.97 --centroids SSU.centroids97 --biomout SSU.biomout --otutabout SSU.otutabout97 --sizeout --threads 15 --uc SSU.uc.97 --strand both
blastn -query SSU.centroids97 -db SILVA_138.1_SSURef_NR99_tax_silva.fasta  -out SSU.silva -outfmt 6 -evalue 0.01 -num_threads 24
blastn -query SSU.centroids97 -db pr2_version_4.14.0_SSU_UTAX.fasta  -out SSU.pr2 -outfmt 6 -evalue 0.01 -num_threads 12
less SSU.v97.silva |cut -f 2 |perl try.pl silva_NR99.taxonomy - |cut -d ' ' -f 2- |paste SSU.silva - >x
less x  |grep -v metagenome |sort -u -k 1,1  >SSU.silva.tophit
less SSU.v97.pr2 |sort -u -k 1,1 >SSU.pr2.tophit
less SSU.silva.tophit |grep -v Eukaryota >xxx
less SSU.silva.tophit |grep -v Eukaryota |cut -f 13 |perl -e 'while(<>){chomp;my @a=split/\;/;my @b;my $n = 1;for($i = 0;$i < 7;$i++){if ($a[$i]){$b[$i] = $a[$i];}else{if ($n == 1){$b[$i] = $b[$i-1]._X;$n++}else{$b[$i] = $b[$i-1]._X;$n++}}}my $str=join(";",@b);print "$str\n";}'|sed 's/_X_X_X_X_X/_XXXXX/g'|sed 's/_X_X_X_X/_XXXX/g' |sed 's/_X_X_X/_XXX/g'|sed 's/_X_X/_XX/g'|perl -ne 'chomp;@ar=split(/\;/,$_);if($ar[-1]=~/.*\_X.*/){print "$_ sp.\n"}else{print "$_\n"}'|paste xxx - |perl -ne 'chomp;@ar=split(/\t/,$_);if(@a=split(/\;/,$ar[-1])){print "$ar[0]\t$ar[1]\t$ar[2]\tK__$a[0]\tP__$a[1]\tC__$a[2]\tO__$a[3]\tF__$a[4]\tG__$a[5]\tS__$a[6]\n"}'  >SSU.silva.prokaryate.taxonomy
less SSU.centroids97 |grep '>'|sed 's/>//g'  >SSU.centroids97.id
less SSU.silva.prokaryate.taxonomy  |cut -f 1 |grep -w -vFf - SSU.centroids97.id  |perl /data/chenkai/software/perl/try.pl  SSU.pr2.tophit - 2>>err |sed 's/:nucl//g'|sed 's/:plas//g' |sed 's/d:Stramenopiles,p:Opalozoa/d:Stramenopiles,p:Ochrophytes/g'|sed 's/d:Stramenopiles,p:Sagenista/d:Stramenopiles,p:Ochrophytes/g'|sed 's/d:Amoebozoa,p:Lobosa/d:Amoebozoa,p:Discosea/g' |sed 's/d:Hacrobia,p:Centroheliozoa/d:Hacrobia,p:Ochrophyta/g'|sed 's/d:Opisthokonta,p:Mesomycetozoa/d:Opisthokonta,p:Opisthokonta_X/g'|sed 's/d:Amoebozoa,p:Conosa/d:Amoebozoa,p:Amoebozoa_X/g'|sed 's/:nucl//g'|sed 's/d:Alveolata,p:Chrompodellids/d:Alveolata,p:Alveolata_X/g'|sed 's/d:Apusozoa,p:Apusomonadidae/d:Apusozoa,p:Apusozoa_X/g'|sed 's/d:Apusozoa,p:Hilomonadea/d:Apusozoa,p:Apusozoa_X/g'|sed 's/d:Hacrobia,p:Telonemia/d:Hacrobia,p:Hacrobia_X/g'|sed 's/d:Stramenopiles,p:Pseudofungi/d:Stramenopiles,p:Oomycota/g'|sed 's/d:Stramenopiles,p:Opalozoa/d:Stramenopiles,p:Ochrophytes/g'|sed 's/d:Stramenopiles,p:Sagenista/d:Stramenopiles,p:Ochrophytes/g'|sed 's/d:Amoebozoa,p:Lobosa/d:Amoebozoa,p:Discosea/g' |sed 's/d:Hacrobia,p:Centroheliozoa/d:Hacrobia,p:Ochrophyta/g'|sed 's/d:Opisthokonta,p:Mesomycetozoa/d:Opisthokonta,p:Opisthokonta_X/g'|sed 's/d:Amoebozoa,p:Conosa/d:Amoebozoa,p:Amoebozoa_X/g'|sed 's/:nucl//g'|sed 's/d:Alveolata,p:Chrompodellids/d:Alveolata,p:Alveolata_X/g'|sed 's/d:Apusozoa,p:Apusomonadidae/d:Apusozoa,p:Apusozoa_X/g'|sed 's/d:Apusozoa,p:Hilomonadea/d:Apusozoa,p:Apusozoa_X/g'|sed 's/d:Hacrobia,p:Telonemia/d:Hacrobia,p:Hacrobia_X/g'|sed 's/d:Stramenopiles,p:Pseudofungi/d:Stramenopiles,p:Oomycota/g'|sed 's/d:Hacrobia,p:Katablepharidophyta/d:Hacrobia,p:Ochrophyta/g'|sed 's/p:Metamonada/p:Fornicata/g'|sed 's/p:Ochrophytes/p:Ochrophyta/g'|sed 's/p:Discoba,c:Euglenida/p:Euglenozoa,c:Euglenida/g'|sed 's/p:Discoba,c:Kinetoplastea/p:Euglenozoa,c:Kinetoplastea/g'|sed 's/p:Discoba,c:Heterolobosea/p:Heterolobosea,c:Heterolobosea/g'|sed 's/p:Discoba/p:Discoba_X/g' |perl -ne 'chomp;@ar=split(/\t/,$_);if(@a=split(/\,/,$ar[1])){if(@b=split(/\;/,$a[0])){print "$ar[0]\t$b[0]\t$ar[2]\t$b[1]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\n"}}' |cut -f 1-12 |sed 's/tax=//g'|sed 's/k:/K__/g' |sed 's/d:/D__/g' |sed 's/p:/P__/g'|sed 's/c:/C__/g' |sed 's/o:/O__/g'|sed 's/f:/F__/g' |sed 's/g:/G__/g'|sed 's/s:/S__/g'  > SSU.pr2.Eukaryota.taxonomy
less  SSU.pr2.Eukaryota.taxonomy |cut -f 1-4,6-|cat - SSU.silva.prokaryate.taxonomy |cut -f 1,4-  >SSU.v97.taxonomy

#5.phylogeny tree
mafft --thread 24   SSU.centroids97 > SSU.centroids97.mafft
Gblocks SSU.centroids97.mafft -t=d -b2=0.85 -b3=8 -b4=2 -b5=a -e=.gb
sed -i 's/ //g' SSU.centroids97.mafft.gb
raxmlHPC-PTHREADS -x 1234567890 -p 1234567890 -f a -# 1000 -T 36 -m GTRGAMMA -s   SSU.centroids97.mafft.gb  -n  SSU.nwk

