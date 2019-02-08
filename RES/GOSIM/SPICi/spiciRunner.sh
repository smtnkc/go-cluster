for D in `find . -maxdepth 1 -mindepth 1 -type d`
do
	if [ ! -d $D/CLUSTERS ] 
	then
    	mkdir -p $D/CLUSTERS
	fi
	
	for f in $D/*.tsv
	do
		echo "$(basename $f)"
		spici -i $f -o "$D/CLUSTERS/$(basename $f)"
		echo ""
	done
done
