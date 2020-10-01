#!/usr/bin/bash

DIR=$1

wd="$( cd "$(dirname "$( readlink -f ${BASH_SOURCE[0]} )")" >/dev/null 2>&1 && pwd)"

rep_head="${wd}/detectISlisting.tex"
rep_css="${wd}/detectISlisting.css"

singimg=${wd}"../../../utils/detectIS.simg"

array=($(ls -d $DIR*.md))

 
for name in "${array[@]}"
do
	name=${name%".md"}
	singularity exec $singimg bash -c "pandoc ${name}.md --listings -H ${rep_head} -o ${name}.pdf" 
	singularity exec $singimg bash -c "pandoc --css=${rep_css} --mathjax  --to=html5 ${name}.md -o ${name}.html"
done 
