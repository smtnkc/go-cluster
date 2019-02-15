for f in ./INPUT/*.tsv
do
	echo "$f"
	spici -i $f -o "./LINE_BY_LINE_CLUSTERS/$(basename $f)" -d 0.5 -s 2 -g 0.5
	echo ""
done
