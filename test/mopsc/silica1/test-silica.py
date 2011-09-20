#!/usr/bin/python2.6
# Copyright William Menz (wjm34) 2011

# Load part data file. Returns primary and collision particle size lists.
def LoadPart (inpfile):
	# Require the regexp module
	import re
	
	data = open(inpfile, 'rb')
	headers = data.readline()
	items = headers.split(",")
	
	m0_id = 4
	n_si_id = 0
	n_o_id = 0
	sl_ld = 0
	dcol_id = 0
	dpri_id = 0
	i = 0
	
	for item in items:
		n_si_m = re.search("Si atoms",item)
		n_o_m = re.search("O atoms",item)
		sl_m = re.search("Sintering Level",item)
		dcol_m = re.search("Collision Diameter",item)
		dpri_m = re.search("primary diameter",item)
		
		if n_si_m:
			n_si_id = i
		if n_o_m:
			n_o_id = i
		if sl_m:
			sl_id = i
		if dcol_m:
			dcol_id = i
		if dpri_m:
			dpri_id = i
		i=i+1
	del i
	print dcol_id
	# Initialise data storage arrays
	m0 = []
	n_si = []
	n_o = []
	sl = []
	dcol = []
	dpri = []
	
	# Split header file by comma
	for csvline in data:
		line = csvline.split(",")
		m0_val = float(line[m0_id].strip())
		m0.append(m0_val)
		
		n_si_val = float(line[n_si_id].strip())
		n_si.append(n_si_val)
		
		n_o_val = float(line[n_o_id].strip())
		n_o.append(n_o_val)
		
		sl_val = float(line[sl_id].strip())
		sl.append(sl_val)
		
		dcol_val = float(line[dcol_id].strip())
		dcol.append(dcol_val)
		
		dpri_val = float(line[dpri_id].strip())
		dpri.append(dpri_val)
	
	# Close file
	data.close()
	
	return (m0, n_si, n_o, sl, dcol, dpri)

#######################################################################
# Main body of script
#######################################################################
import sys
# Load data
try:
	part = "silica-part.csv"
	data_part = LoadPart(part)
except:
	print "Required files could not be found or loaded"
	sys.exit(2)

# Average values at final time found at end of part file
dcol	= 1.0e9*data_part[4][len(data_part[4])-1]
dpri	= data_part[5][len(data_part[4])-1]
m0		= data_part[0][len(data_part[0])-1]
n_si	= data_part[1][len(data_part[1])-1]
n_o		= data_part[2][len(data_part[2])-1]
si_to_o = n_si/n_o
sl		= data_part[3][len(data_part[3])-1]

# Define test values and error lists (absolute)
# Order is dcol dpri m0 si_to_o n_sl
data_labs = ['dcol','dpri','m0','si_to_o','sl']
data_this = [dcol,dpri,m0,si_to_o,sl]
data_comp = [74.5,7.27,2.15e14,0.545,0.103]
# Errors list represents 1/2 distance to 99.9% CI
# Found for running the test case 5 times
# i.e. min bound = mean - err, max = mean + err for full 99.9% CI.
data_errs = [2.65,1.12,1.5e13,0.027,0.022]
data_fail = [0,0,0,0,0]

# Check if the test results exceed error margins
for lab, this, comp, err in zip(data_labs, data_this, data_comp, data_errs):
	if (this >= comp-err and this <= comp+err):
		print("[PASS]: {0} returned {1}".format(lab,this))
	else:
		print("[FAIL]: {0} returned {1}, expected in range {2}--{3}".format(lab,this,comp-err,comp+err))
		data_fail[data_labs.index(lab)] = 1

# Exit script and inform of failure
fail = max(data_fail)
if fail > 0:
	print("************")
	print("TEST FAILURE")
	print("************")
	sys.exit(2)
else:
	print("All tests passed :D")
	sys.exit(0)
