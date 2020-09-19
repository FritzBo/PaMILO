#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

rm $2 2> /dev/null
touch $2
while IFS= read -r line; do
	if grep -q "obj:" <<< "$line"; then
		declare -i i=1
		my_array=($(printf "$line\n" | tr " " "\n"))
		for substr in "${my_array[@]}"; do
			if [[ $substr == "z"* ]]; then
				printf " obj$i:\r\n  z$i\r\n" >> $2
				i=$((i+1))
			fi
		done
	elif grep -q "imize" <<< "$line"; then
		line=$(echo $line | sed 's/\r//g')
		printf "$line multi-objectives\r\n" >> $2
	else
		echo "$line" >> $2
	fi
done < $1

