#from collections import defaultdict
from collections import OrderedDict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os, sys, argparse, csv, json, yaml


""" Command line parsers """
parser = argparse.ArgumentParser()
parser.add_argument("files",nargs='+',help="JSON files to merge")
parser.add_argument("-c", "--config", help="Config YAML file")
args = parser.parse_args()

if args.config is not None:
    with open(args.config, 'r') as ymlfile:
        cfg = yaml.safe_load(ymlfile)

def main(files):
	""" This function merges a series of python dicts stored as json files supplied as arguments
		into a large dict, one at a time, possibly with Nones for missing values. 
	"""
	
	dict1 = parseDict(files[0])
	for i in range(1, len(files)):
		dict2 = parseDict(files[i])

		# add new values for the keys that exist in both dicts
		output = dict((k, newval(dict1[k], dict2.get(k))) for k in dict1)

		# update dict for keys that exist in the second dict but not the first
		output.update((k, [None]*i + [dict2[k]]) for k in dict2 if k not in dict1)

		# rewind and roll again
		dict1 = output

	# sort the dict entities based on the sum of all their occurrances
	odict = OrderedDict(sorted(dict1.items(), reverse=True, key = lambda e: sum(filter(None, e[1]))))

	records = (genSeqRec(seq,val) for seq,val in odict.iteritems())
	SeqIO.write(records, "consolidated_seqs.fasta", "fasta")

	rows = ([x] + y for x,y in odict.items())
	headers = ["Seq"]
	if args.config is not None:
		for s,e in OrderedDict(sorted(cfg["Samples"].items())).items():
			headers += [e]
	else:
		for i in range(0, len(files)):
			headers += ["Sample" + str(i + 1)]

	with open('consolidated_seqs.csv', 'w') as fp:
		w = csv.writer(fp)
		w.writerow(headers)
		w.writerows(rows)


	with open('consolidated_seqs.json', 'w') as fp:
		json.dump(odict, fp, indent=4)
		

def newval(v1, v2):
	if isinstance(v1,list):
		return v1 + [v2]
	else: 
		return [v1, v2]

def parseDict(f):
	data = {}
	with open(f, 'r') as fp:
		data = json.load(fp)
	return data

def genSeqRec(seq, vals):
	try:
		genSeqRec.counter += 1
	except AttributeError:
		genSeqRec.counter = 1

	seq_id = "ab_seq" + str(genSeqRec.counter)
 
	seq_desc = "Counts:" + "/".join(map(lambda x: '0' if x is None else str(x), vals))
	return SeqRecord(Seq(seq), id = seq_id , description = seq_desc)

if __name__ == "__main__":
	main(args.files)
