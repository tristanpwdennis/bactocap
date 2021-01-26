for f in *total.txt
do
	for i in `seq 0 5 200`
	do
		echo `awk -v num=$i '$6 < num {print$4}' $f | uniq | wc -l ` `basename $f .total.txt` $i >> covstats.txt
	done
done
