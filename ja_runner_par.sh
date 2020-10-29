#!/usr/local/bin/bash

# This is the runner file run by bsub
# Arguments: runner.sh filelist
# Environment variables: LSB_JOBINDEX
#Example command:
#mkdir -p LOGS;size=`wc -l result.list|cut -f 1 -d " "`;bsub -J "p_check[1-${size}]" -o "LOGS/%J_p_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh ~/Work/bash_scripts/check_pvals_test.sh result.list
#mkdir -p LOGS;size=`wc -l match_file.list|cut -f 1 -d " "`;bsub -J "sh_dens[1-${size}]" -o "LOGS/%J_sh_dens.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh ~/Work/bash_scripts/sharingDensity.sh match_file.list map.list
echo "${@}"
while getopts ":dstql" opt; do
  case $opt in
    d)
      echo $opt
      echo "Double column list mode triggered!" >&2
      file=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $1}'`
      file2=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $2}'`
      echo ${file}
      echo ${file2}
      script=$2
      $script ${file} ${file2} "${@:4}"
      ;;
    t)
      echo $opt
      echo "Triple column list mode triggered!" >&2
      file1=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $1}'`
      file2=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $2}'`
      file3=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $3}'`
      echo ${file1}
      echo ${file2}
      echo ${file3}
      file=${file1}\:${file2}\:${file3}
      script=$2
      $script ${file} "${@:4}"
      ;;
    q)
      echo $opt
      echo "Quadruple column list mode triggered!" >&2
      file1=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $1}'`
      file2=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $2}'`
      file3=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $3}'`
      file4=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $4}'`
      echo ${file1}
      echo ${file2}
      echo ${file3}
      echo ${file4}
      file=${file1}\:${file2}\:${file3}\:${file4}
      script=$2
      # $script ${file} "${@:4}"
      $script ${file1} ${file2} ${file3} ${file4} "${@:4}"
      ;;
    s)
      echo $opt
      echo "Single list mode triggered!!" >&2
      file=`sed -n "${LSB_JOBINDEX}p" $3`
      # echo ${file}
      script=$2
      # $script ${file} $4 $5 $6 $7 $8
      $script ${file} "${@:4}"
      ;;
    l)
      echo $opt
      echo "Script list mode triggered!!" >&2
      file=`sed -n "${LSB_JOBINDEX}p" $2`
      # echo ${file}
      script=$file
      # $script ${file} $4 $5 $6 $7 $8
      $script "${@:3}"
      ;;
    *)
      echo $opt
    ;;
  esac

done
