#!/usr/bin/env bash
#
#generate chunks for panel merging based on current legends
function find_min(){
pos_array=($@)
a_min=${pos_array[0]}

# Loop through all elements in the array
for i in "${pos_array[@]}"
do
    # Update min if applicable
    if [[ "$i" -lt "$a_min" ]]; then
        a_min="$i"
    fi
done
echo ${a_min}

}

function find_max(){
pos_array=($@)
a_max=${pos_array[0]}

# Loop through all elements in the array
for i in "${pos_array[@]}"
do
    # Update max if applicable
    if [[ "$i" -gt "$a_max" ]]; then
        a_max="$i"
    fi
done

echo ${a_max}
}

###############################################################
chunk_mode=$1
infile=$2
chunk_size=$3
chr=$4
outfile=$5

case ${chunk_mode} in
	vcf )
	# bcftools view 
	start_chr=`bcftools view -H ${infile} | head -1 | awk '{print $2}'`
	end_chr=`bcftools view -H ${infile} | tail -n1 | awk '{print $2}'`
	;;
    range )
    region=$2
    start_chr=`echo ${region} | cut -f 1 -d "-"`
    end_chr=`echo ${region} | cut -f 2 -d "-"`
    ;;
	*)
	# s2=`zcat ${leg_2} | head -2 | tail -n1 | awk '{print $2}'`
	# starts=(${s1} ${s2})

	# e1=`zcat ${leg_1} | tail -n1 | awk '{print $2}'`
	# e2=`zcat ${leg_2} | tail -n1 | awk '{print $2}'`
	# ends=(${e1} ${e2})

	# start_chr=`find_min "${starts[@]}"`
	# end_chr=`find_max "${ends[@]}"`

	;;
esac

start_pos=${start_chr}
end_pos=0
# than create the chunk file to use to submit the job array
while [ ${start_pos} -lt ${end_chr} ]
do
end_pos=$[start_pos + chunk_size]
if [ ${end_pos} -ge ${end_chr} ];then
end_pos=${end_chr}
fi
echo "${chr} ${start_pos} ${end_pos}"
start_pos=$[end_pos + 1 ] 
done > ${outfile}