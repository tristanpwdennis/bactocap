for f in *
do
	bedtools coverage -d -a ../ancillary/annotation/aln.Myco.4columns.bed -b $f/$f".rg.bam" > $f/$f".total.txt"
	bedtools coverage -d -a ../ancillary/annotation/complement.aln.Myco.bed -b $f/$f".rg.bam" > $f/$f".complement.total.txt"
done
