process CHECK_STRANDNESS {
	input:
		tuple val(sample_name), path(reads)
    		path(reference_cdna)
    		path(annotation)

	output: 
		env STRANDNESS
    require([
	FOR_ALL("f", ITER("*.fastq"), 
		{ f -> 
			IF_THEN(
				NOT(
					EQUAL(
						NUM("\$(( \$(wc -l $f | cut -d' ' -f1)/4*4 ))"), 
						NUM("\$(wc -l $f | cut -d' ' -f1)")
					)
				), 
				"exit 1"
			)
		}
	), 
	"""for file in *.fastq; do 
		if ! awk 'BEGIN{c=0;} {if (c % 4 == 0 && substr(\$0,1,1) != "@" || substr(\$0,1,1) == "@ ") {exit 1} if (c % 4 == 1) {len = length(\$0)} if (c % 4 == 2 && substr(\$0,1,1) != "+") {exit 1} if (c % 4 == 3 && len != length(\$0)) {exit 1} c = c+1}' \$file; then 
			exit 1; 
		fi;
	done"""
])
	promise(['if [ "$STRANDNESS" = "error" ]; then exit 1; fi'])

	shell:
	'''
	check_strandedness -g !{annotation} -r1 !{reads[0]} -r2 !{reads[1]} --transcripts !{reference_cdna} > result.txt
	result=$( tail -n 1 result.txt )
	if [[ $result == *"likely unstranded"* ]]; then
  		STRANDNESS="unstranded"
	elif [[ $result == *"likely RF/fr-firststrand"* ]]; then
  		STRANDNESS="firststrand"
	elif [[ $result == *"likely FR/fr-secondstrand"* ]]; then
  		STRANDNESS="secondstrand"
	else
		STRANDNESS="error"
	fi
    #    STRANDNESS="firststrand"
     '''
}

