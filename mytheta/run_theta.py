# !/bin/python
import subprocess
import glob
import copy
import re
import sys

from optparse import OptionParser

parser = OptionParser()

parser.add_option('--input', metavar='F', type='string', action='store',
				  dest='input',
				  help='input files to be globbed')

(options, args) = parser.parse_args()

argv = []

files = glob.glob(options.input)
print files

fname = 'driver.py'

for file in files:

	#So we can keep the report
	if file.endswith('.root'):
		dirname = file[:-5]
	else:
		dirname = file

	#Create the driver file that theta expects
	with open(fname, 'w') as fout:
		fout.write(
			"import sys\n"+
			"sys.path.append('/uscms_data/d3/dfehling/theta')\n"+
			"sys.path.append('/home/dfehling/work/theta/')\n"+
			"from base_theta import *\n\n"+
			"had_xs_model('"+str(file)+"', 245)")

	#Commands to be executed
	commands = [
		'rm analysis.py',
		'rm -rf analysis/',
		'ln -s ' + './driver.py ' + ' ./analysis.py',
		'./utils2/theta-auto.py',
		# 'rm -rf analysis_' + dirname
		'mv histos-' + file + ' analysis/',
		'mv analysis analysis_' + dirname
		]

	#Execute the commands
	for s in commands :
		print 'executing ' + s
		subprocess.call( [s], shell=True )

