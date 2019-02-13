for f in ./*.tsv
do
	echo "$f"
	spici -i $f -o "./CLUSTERS/$(basename $f)" -d 0.5 -s 2 -g 0.5
	echo ""
done
