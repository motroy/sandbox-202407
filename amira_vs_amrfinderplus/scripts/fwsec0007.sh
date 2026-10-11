#!/usr/bin/env bash
# SAMN08768108 (E. coli O103:H2 FWSEC0007): ONT SRR7896269 + MiSeq SRR7896254, truth = GCF_005037795.1 (complete)
# usage: scripts/fwsec0007.sh {download|myloasm|polish|amrfinder|amira}   (run from amira_vs_amrfinderplus/)
set -euo pipefail
D=data/FWSEC0007; R=results/FWSEC0007; T=${THREADS:-4}
ONT=$D/SRR7896269_1.fastq.gz; R1=$D/SRR7896254_1.fastq.gz; R2=$D/SRR7896254_2.fastq.gz
mkdir -p $D $R
case "$1" in
download)
  B=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789
  [ -s $ONT ] || curl -fsSL -o $ONT $B/009/SRR7896269/SRR7896269_1.fastq.gz
  [ -s $R1 ] || curl -fsSL -o $R1 $B/004/SRR7896254/SRR7896254_1.fastq.gz
  [ -s $R2 ] || curl -fsSL -o $R2 $B/004/SRR7896254/SRR7896254_2.fastq.gz
  [ -s $D/truth.fna ] || curl -fsSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/037/795/GCF_005037795.1_ASM503779v2/GCF_005037795.1_ASM503779v2_genomic.fna.gz | gzip -dc > $D/truth.fna ;;
myloasm)
  # defaults. (Tried --quality-value-cutoff 85 to keep more R9 reads: worse, 359 contigs/9.5 Mb/N50 43 kb vs 273/8.0 Mb/66 kb)
  pixi run -e asm myloasm $ONT -o $R/myloasm -t $T --clean-dir
  ls $R/myloasm ;;
polish)
  A=$(ls $R/myloasm/assembly_primary.fa); W=$R/polish; mkdir -p $W; cp $A $W/myloasm.fasta
  pixi run -e asm bwa index $W/myloasm.fasta
  pixi run -e asm bwa mem -t $T -a $W/myloasm.fasta $R1 > $W/a1.sam
  pixi run -e asm bwa mem -t $T -a $W/myloasm.fasta $R2 > $W/a2.sam
  pixi run -e asm polypolish filter --in1 $W/a1.sam --in2 $W/a2.sam --out1 $W/f1.sam --out2 $W/f2.sam
  pixi run -e asm polypolish polish $W/myloasm.fasta $W/f1.sam $W/f2.sam > $W/polypolish.fasta
  rm -f $W/a1.sam $W/a2.sam $W/f1.sam $W/f2.sam
  pixi run -e asm pypolca run -a $W/polypolish.fasta -1 $R1 -2 $R2 -t $T -o $W/pypolca --careful -f
  cp $W/pypolca/*_corrected.fasta $W/polypolish_pypolca.fasta ;;
amrfinder)
  mkdir -p $R/amrfinder
  for n in truth:$D/truth.fna myloasm:$R/polish/myloasm.fasta polypolish:$R/polish/polypolish.fasta polished:$R/polish/polypolish_pypolca.fasta; do
    pixi run -e amrfinder amrfinder -n ${n#*:} -O Escherichia --plus --threads $T -o $R/amrfinder/${n%%:*}.tsv
  done ;;
amira)
  pixi run -e amira amira --reads $ONT --output $R/amira --species Escherichia_coli \
    --panRG-path data/Escherichia.coli.panidx.zip --pandora-path "$PWD/bin/pandora" --cores $T ;;
esac
