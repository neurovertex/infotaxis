#!/bin/bash

output=found.gif
if [[ $# > 0 ]]; then
	output="$1"
fi

dir="pictures"
if [[ $# > 1 ]]; then
	dir="$2"
fi

lastframe=$( ls $dir/*.png | tail -n 1 )
firstsize=$( file "$dir/iteration0000.png" | cut -d" " -f 5,7 )
lastsize=$( file "$lastframe" | cut -d" " -f 5,7 )
[[ "$firstsize" != "$lastsize" ]] && echo "Error : image sizes differ ($firstsize / $lastsize). Resize images first" && exit 1
convert -delay 10 $dir/iteration*.png -delay 400 $lastframe $output && \
echo "Done! Wrote gif to $output"
