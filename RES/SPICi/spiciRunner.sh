for f in ./../GOSIM/$1/forSPICi/*.tsv
do
	echo "$f"
	t_start=`date +%s%3N`
	/usr/bin/time --format='%M' spici -i $f -o "./LINE_BY_LINE_CLUSTERS/$1/$(basename $f)" -d 0.5 -s 2 -g 0.5
	t_end=`date +%s%3N`
	echo `expr $t_end - $t_start` # ms
done
