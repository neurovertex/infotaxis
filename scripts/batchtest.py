#!/bin/python3
import json
import os
import subprocess
from sys import argv
from math import sqrt

if argv[1] == 'test':
	prog = "./bin/main"
	args = " --size=100:100 -a 0.1 --windvel=1 --windangle=-135 --diff=1 --resolution=10 -q -n "
	logarg = " --log=batch/log_%d_%d.txt"
	startarg = " --start=%d:%d"
	sourcearg = " --source=%d:%d"
	source =90
	nbiter = int(argv[2]) if len(argv) > 2 else 1
	incr = int(argv[3]) if len(argv) > 3 else 10
	fastarg = " --fast" if len(argv) > 4 and argv[4] == "fast" else ""
	filename = argv[5] if len(argv) > 5 else "batch/result.json"

	obj = {}

	print("Starting batch tests. %d tests per distance" % nbiter)
	for dist in range(10, 81, incr):
		arr = []
		for i in range(0, nbiter):
			arglist=prog + args + sourcearg % (source, source) + \
					startarg % ((source-dist), (source-dist)) + \
					fastarg + logarg % (dist, i)
			print(arglist)
			count=int(subprocess.check_output(arglist, shell=True, universal_newlines=True))
			print("Count : \"%s\"" % count)
			arr.append(count)
		obj[dist] = arr

	json=json.dumps(obj, sort_keys=True, indent=4, separators=(',', ': '))
	f=open(filename, mode='w')
	f.write(json)
	f.write('\n')
	f.close()
elif argv[1] == 'process':
	infile = argv[2] if len(argv) > 2 else "batch/result.json"
	outfile = argv[3] if len(argv) > 3 else "batch/data.txt"
	dotfile = argv[4] if len(argv) > 4 else None
	histfile = argv[5] if len(argv) > 5 else None

	f=open(infile, mode='r')
	obj = json.load(f)
	f.close()
	
	dot = open(dotfile, mode='w') if dotfile != None else None
	hist = open(histfile, mode='w') if histfile != None else None
	
	histo = {}

	f=open(outfile, mode='w')
	f.write("dist average median deviation success rate\n")
	for dist in sorted(obj):
		filtered=sorted([i for i in obj[dist] if i <= 2048])
		success=len(filtered) * 100. / len(obj[dist])
		avg=sum(filtered)/len(filtered)
		stddev=sqrt(sum([(i-avg)**2 for i in filtered])/len(filtered))
		f.write("%d %d %d %f %f\n" % (int(dist), avg, filtered[int(len(filtered)/2)], stddev, success))
		for i in filtered:
			v = int(i * 10 / int(dist))
			histo[v] = histo[v] + 1 if v in histo.keys() else 1
			if dot != None:
				dot.write("%s %d\n" % (dist, i))
	
	if hist != None:
		for i in range(0, max(histo.keys())+1):
			hist.write("%d %d\n" % (i, histo[i] if i in histo.keys() else 0))
	
	f.close()
	if dot != None:
		dot.close()
	if hist != None:
		hist.close()
	
	print("Success. Wrote to "+ outfile)
else:
	print("Command "+ argv[1] +" unrecognized")