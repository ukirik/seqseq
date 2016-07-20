from __future__ import print_function
from Bio import Seq, SeqIO
from itertools import islice, izip
from pprint import pprint
from collections import namedtuple, Counter
import os, sys, gzip, bz2, argparse, json, yaml

""" Command line parsers """
parser = argparse.ArgumentParser()
parser.add_argument("fprime", help="5 prime read")
parser.add_argument("tprime", nargs='?', help="3 prime read", default=None)
parser.add_argument("-o", "--output", help="Output file")
parser.add_argument("-c", "--config", help="Configuration file in YAML syntax")
parser.add_argument("-a", "--append_summary", help="A file to append the individual read statistics")
parser.add_argument("-n", "--nbr_of_seqs", type=int, help="number of sequences to process")
parser.add_argument("-t", "--no_trim", action="store_true")
parser.add_argument("-s", "--save_seqs", choices=['lowq','stop', 'all'], \
                        help="whether or not to save filtered sequences")

args = parser.parse_args()

if args.config is not None:
    with open(args.config, 'r') as ymlfile:
        cfg = yaml.safe_load(ymlfile)


""" Variable region definitions 
    TODO: Consider moving them as class members instead
    TODO: These should be user supplied parameters
"""
Region = namedtuple('Region', ['name', 'dir', 'start', 'stop'])
hm_regions = []
if args.config is not None:
    for name,vals in cfg["Regions"].items():
        hm_regions.append(Region(name, vals[0], vals[1], vals[2]))
else:
    hm_regions.append(Region("cdrh1", 5,  0*3, 10*3))
    hm_regions.append(Region("cdrh2", 5, 23*3, 33*3))
    hm_regions.append(Region("cdrh3", 5, 69*3, 87*3))
    hm_regions.append(Region("cdrl3", 3, 33*3, 43*3))

""" String representation of the linker between VH and VL chains
    Note that the parts of VH after CDRH3 and half of VL, containing 
    CDRL1 and CDRL2 are a part of this string as per default QC trimming
"""
linker_str = cfg["Stitching"]["linker"] if args.config is not None else "~"


""" Various counters used by methods """
nseqs = 0
starSeqs = 0
lowQSeq = 0

magic_dict2 = {
    "\x1f\x8b\x08": lambda f : gzip.open(f,'rb'),
    "\x42\x5a\x68": lambda f : bz2.open(f,'rb'),
    }

def filter_on_average(rec, q = 20):
    """ Returns True if the supplied sequence record passes the quality threshold.
        Default threshold is 20.
    """
    qScores = rec.letter_annotations["phred_quality"]
    if (sum(qScores)/len(qScores)) < q:
        return False
    
    return True
        
def filter_on_cdr(rec, direction = 5, q = 20, combined = False):
    """ Returns True if the supplied sequence record passes the given quality threshold,
        specifically on CDR regions. Default threshold is 20. 
        Ideally the start and stop positions should be read in from the user as exp parameters
    """
    try:
        if combined:
            x = None
            for r in filter(lambda r: r.dir == direction, hm_regions):
                x = rec[r.start : r.stop] if x is None else x + rec[r.start : r.stop]
            return filter_on_average(x, q)
        else:
            passQC = True
            for r in filter(lambda r: r.dir == direction, hm_regions):
                passQC = passQC and filter_on_average(rec[r.start : r.stop],q)
            return passQC
    except ZeroDivisionError:
        pprint(hm_regions)
        userinput = raw_input("Press any button to continue")

def filter_quality(rec, direction = 5, q = 20, method='ave'):
    """ Applies the appropriate quality filtering method 
        TODO: Implement conditional storage and printing of filtered sequences
    """
    if method == 'cdr':
        result = filter_on_cdr(rec, direction, combined = False)
    elif method == 'combi':
        result = filter_on_cdr(rec, direction, combined = True)
    else:
        result = filter_on_average(rec)

    if result == False: 
        global lowQSeq
        lowQSeq += 1

    return result

def filter_seq(rec, direction = 5):
    """ Applies all desired filters to SeqRecord 
        All settings from the config file should be applied here!
    """
    
    hasStopCodon = '*' in str(rec.seq.translate())
    
    if hasStopCodon:
        global starSeqs
        starSeqs += 1
        return False
    else:
        if args.config is not None:
            qthreshold = cfg["Quality"]["threshold"]
            qmethod = cfg["Quality"]["filtering"]
            return filter_quality(rec, direction, q = qthreshold, method = qmethod)
        else:
            return filter_quality(rec, direction, method = 'cdr')
    
def check_filetype(filename):
    """ Check magic bytes to figure out the filetype """
    max_len = max(len(x) for x in magic_dict2)
    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict2.items():
        if file_start.startswith(magic):
            return filetype(filename)
    return filename

def single_read(read1, direction = 5, nbrofitems = 10**8, fileout = None):
    """ Procedure for reading both sequences and stitching them together 
        Unless specified, it will read 10^8 sequences from the supplied read
    """
    seqFreqs = Counter()

    # TODO: Enfore trimming parameters (or rather YAML config file)
    if cfg is not None:
        trim5 = cfg["Trim"]["fwdread"]
        trim3 = cfg["Trim"]["revread"]
    else:
        trim5 = [27,None]
        trim3 = [21, 150]

    for rec in islice(read1, nbrofitems):
        if(direction == 5):
            rec = rec[trim5[0] : trim5[1]]                          # Trim the primer variable sequence
        else:   
            rec = rec[trim3[0] : trim3[1]].reverse_complement()  # Trim the low Q half of the 3' read, the primer AND take rev complement
            
        aaSeq = rec.seq.translate()
        if filter_seq(rec, direction) :
            seqFreqs.update({ str(aaSeq) : 1 }) 
        
        global nseqs 
        nseqs += 1

    if args.no_trim is not True:
        """ Trim out sequences that occur just once """
        seqFreqs = seqFreqs - Counter(k for k in seqFreqs.keys())

    if fileout is not None:
        fout = open(fileout, "w")
        sys.stdout = fout
        jsonf = os.path.join(os.path.split(fileout), "seqdata.json")
        with open(jsonf, 'w') as fp:
            json.dump(seqFreqs, fp, indent=4)


    pprint(seqFreqs.most_common(100), width = 120)

    if fileout is not None:
        sys.stdout = sys.__stdout__
        fout.close()

def paired_read(read1, read2, nbrofitems = 10**8, fileout = None):
    """ Procedure for reading both sequences and stitching them together 
        Unless specified, it will read 10^8 sequences from the supplied reads
    """
    seqFreqs = Counter()

    # TODO: Enfore trimming parameters (or rather YAML config file)
    if args.config is not None:
        trim5 = cfg["Trim"]["fwdread"]
        trim3 = cfg["Trim"]["revread"]
    else:
        trim5 = [27,None]
        trim3 = [21, 150]

    for rec1, rec2 in islice(izip(read1, read2), nbrofitems):

        rec1 = rec1[trim5[0] : trim5[1]]                         # Trim the primer variable sequence
        rec2 = rec2[trim3[0] : trim3[1]].reverse_complement()    # Trim the low Q half of the 3' read, the primer AND take rev complement

        global nseqs 
        nseqs += 1

        if filter_seq(rec1, direction=5) and filter_seq(rec2, direction=3):
            aa1 = rec1.seq.translate()
            aa2 = rec2.seq.translate()

            # Stitch the strings together
            if args.config is not None:
                i = str(aa1).rfind(cfg["Stitching"]["f_anchor"])
                j = str(aa2).find(cfg["Stitching"]["r_anchor"])
                
                # Check whether or not stitching is done in the expected place
                # TODO: this should be done in a more graceful way
                if i < len(str(aa1)) * 0.75:
                    print("Warning: linker anchor on VH side not found where it was expected (i = {})".format(i))
                    print("read1: {} (i = {})".format(str(aa1), i))

                if j > len(str(aa2)) * 0.25:
                    print("Warning: linker anchor on VL side not found where it was expected (j = {})".format(j))
                    print("read2: {} (j = {})".format(str(aa2),j))
                    
            else:
                i = None
                j = None

            aakey = str(aa1)[:i] + linker_str + str(aa2)[j:]
            seqFreqs.update({ aakey : 1 })  

    if args.append_summary is not None:
        """ Export read stats before trimming sequences that occur just once """ 
        filtseqs = sum(seqFreqs.values())
        dist_seqs = len(list(seqFreqs))

        promille_seqs = 0
        for k,v in islice(seqFreqs.most_common(), 1000):
            if v > filtseqs / 1000:
                promille_seqs +=1 
            else:
                break

        with open(args.append_summary, 'a') as statfile:
            print(os.path.dirname(fileout), nseqs, lowQSeq, starSeqs, filtseqs, dist_seqs, promille_seqs, sep="\t", file=statfile)

    if args.no_trim is not True:
        """ Trim out sequences that occur just once """
        seqFreqs = seqFreqs - Counter(k for k in seqFreqs.keys())

    if fileout is not None:
        fout = open(fileout, "w")
        sys.stdout = fout

        outdir = os.path.dirname(fileout)
        jsonf = os.path.join(outdir, "seqdata_paired.json")

        with open(jsonf, 'w') as fp:
            json.dump(seqFreqs, fp, indent=4)

    pprint(seqFreqs.most_common(100), width = 240)
    
    if fileout is not None:
        sys.stdout = sys.__stdout__
        fout.close()

def main():
    """ Parses the paired FASTQ files and does a freq table """
    
    read1 = None
    read2 = None

    file1 = args.fprime
    file2 = args.tprime if args.tprime is not None else None
    f = args.output
    n = args.nbr_of_seqs

    
    read1 = SeqIO.parse(check_filetype(file1), 'fastq')

    if file2 is not None:
        read2 = SeqIO.parse(check_filetype(file2), 'fastq')
        paired_read(read1,read2, nbrofitems=n, fileout=f)
    else:
        single_read(read1, direction=5, nbrofitems=n, fileout=f)

if __name__ == "__main__":
    main()
