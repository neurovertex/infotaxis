#!/bin/bash

output=found.gif
if [[ $# > 0 ]]; then
	output=$1
fi

lastframe=pictures/$( ls pictures | tail -n 1 )
convert -delay 10 pictures/iteration*.png -delay 400 $lastframe $output && \
echo "Done! Wrote gif to $output"
