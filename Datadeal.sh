#!/bin/bash

ls fst/*Backgrond.windowed.weir.fst | while read line
do 
	id=${line#*/}
	id=${id%%.*}
	awk '{print $6}' fst/${id}.Backgrond.windowed.weir.fst | sed '1d' > File/${id}.back
	awk '{print $6}' fst/${id}.windowed.weir.fst | sed '1d' > File/${id}.gene
	#paste File/${id}.back File/${id}.gene > File/${id}.fst
done

