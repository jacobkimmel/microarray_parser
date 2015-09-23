
'''
Control - Exp:
negative numbers = higher expression

Highlight gene if:
2 out of three animals demonstrate a similar directionality 

If both are increased in expression, label "Upregulated" with an average of the two genes in one direction presented, as well as an average of all three genes
'''

'''
USAGE:

#Reports all genes, regardless of sample variability
	python microarray_parser.py -a samples.csv output.csv 

#Only reports genes with lower variability
	python microarray_parser.py -l samples.csv output.csv
'''

import csv
import sys
import numpy as np
import operator

def CommonDirection(x, y, z):
	negatives = []
	positives = []
	for var in (x, y, z):
		if var < 0:
			negatives.append(var)
		else:
			positives.append(var)
	if len(negatives) > 1:
		avg = np.mean(negatives)
		included_vars = negatives
	else:
		avg = np.mean(positives)
		included_vars = positives
	included_vars = [str(x) for x in included_vars]
	return avg, ' '.join(included_vars)

def calculate_variability(x, y, z):
	a = np.array( [ x, y, z] )
	SD = np.std(a)
	return SD

raw_file = open(sys.argv[2], 'rb')

with raw_file as f:
	reader = csv.reader(f, delimiter=',' , skipinitialspace = True)
	next(reader)
	Row_num,Row_name,ACCNUM,SYMBOL,GO,Exp1,Exp2,Exp3 = zip(*reader)

Exp1 = [float(x) for x in Exp1]
Exp2 = [float(x) for x in Exp2]
Exp3 = [float(x) for x in Exp3]
GO = list(GO)
GO[:] = [s.replace(",", ";") for s in GO]

avg_expression = []
mice_used = []
standard_deviations = []
i = 0
while i < len(Exp1):
	exp, mice = CommonDirection(Exp1[i], Exp2[i], Exp3[i])
	SD = calculate_variability(Exp1[i], Exp2[i], Exp3[i])
	avg_expression.append(exp)
	mice_used.append(mice)
	standard_deviations.append(SD)
	i += 1

genes_with_averages = []

i = 0
while i < len(SYMBOL):
	line = (SYMBOL[i], GO[i], Exp1[i], Exp2[i], Exp3[i], avg_expression[i], mice_used[i], standard_deviations[i])
	genes_with_averages.append(line)
	i += 1

genes_with_lowvar = []
i = 0

percentile80th_SD = np.percentile( np.array(standard_deviations), 95)
print percentile80th_SD

while i < len(SYMBOL):
	if standard_deviations[i] < percentile80th_SD: 
		line = (SYMBOL[i], GO[i], Exp1[i], Exp2[i], Exp3[i], avg_expression[i], mice_used[i], standard_deviations[i])
		genes_with_lowvar.append(line)
	else:
		None
		#print "Skipped line ", i, "with SD = ", standard_deviations[i]
	i += 1

#Sort list based on average relative expression from selected samples
#Sorts with negative numbers (highest exp) at the top
genes_with_averages = sorted(genes_with_averages, key = operator.itemgetter(5))
genes_with_lowvar = sorted(genes_with_lowvar, key = operator.itemgetter(5))

#for line in genes_with_averages:
#	print line

#for line in standard_deviations:
#	print line

for line in genes_with_lowvar:
	print line

# Use CSV library to write to CSV
# write_all_genes() writes every calculated gene with average expression to a csv
# write_low_variability writes only genes with the selected low variability criteria

def write_all_genes():
	out_file = open( sys.argv[3], 'wb' )
	with out_file as f:
		selectric = csv.writer(out_file, dialect='excel', quoting=csv.QUOTE_NONE)
		selectric.writerow( ["SYMBOL", "GO", "Exp1", "Exp2", "Exp3", "SelectiveAVG", "SelectedMice", "SD"] )
		for row in genes_with_averages:
			selectric.writerow(row)

def write_low_variability_genes():
	out_file = open( sys.argv[3], 'wb' )
	with out_file as f:
		selectric = csv.writer(out_file, dialect='excel', quoting=csv.QUOTE_NONE)
		selectric.writerow( ["SYMBOL", "GO", "Exp1", "Exp2", "Exp3", "SelectiveAVG", "SelectedMice", "SD"] )
		for row in genes_with_lowvar:
			selectric.writerow(row)


variability = sys.argv[1]
unacceptable_flag_error = "Unacceptable flags entered. \n Please enter -a for all genes, -l for low variability"

if variability == "-a":
	write_all_genes()
elif variability == "-l":
	write_low_variability_genes()
else:
	print unacceptable_flag_error


