import h5py
import numpy as np
import os
import io
import pandas as pd
import subprocess
import argparse
#from make_pentaCY import make_penta

parser = argparse.ArgumentParser(description="Produce five channel representation files from sorted bams")

parser.add_argument('-l', '--positions_file', type=str, 
            help='Optional file of reference position to pass to samtools mpileup.')

parser.add_argument('-f', '--reference', type=str, help="The indexed reference genome each bam was aligned to.")

parser.add_argument('-i', '--bam_infile', type = str, help = "Sorted bam file")

parser.add_argument('-q', '--qual_min', type = int, help = "Phred quality minimum to count a base (default is 50)", default=50)

parser.add_argument('-b', '--buffer', type = float, help = "How many elements in the list before writing to file", default=1e6)

parser.add_argument('-o', '--h5_outfile', type = str, help = "Name of HDF5 outfile -- will default to prefix.h5, where infiles are name {prefix}.sorted.bam")

args = parser.parse_args()

if args.h5_outfile:
    h5_outfile = args.h5_outfile
else:
    h5_outfile = args.bam_infile.replace(".sorted.bam", "h5")

##################
# CORE FUNCTIONS #
##################

def parse_line(pileup):
    mp = pileup.strip().split("\t")
    chrom, pos, ref = mp[0:3] #site data
    pop_bam = mp[3:]
    idx = list(range(0,len(pop_bam), 3))
    site_dict = {"chrom": chrom, "ref": ref.upper(), "pos": pos, "pop_bam": pop_bam, "idx":idx}
    return site_dict


def init_penta(bam_in, h5_out):
    with h5py.File(h5_out, 'w') as f:

        #channels = f.create_dataset('channels', (0, 0), maxshape=(None,None), dtype='f8', chunks=True,  compression="gzip", compression_opts=9)
        #makers = f.create_dataset('markers', (0, 0), maxshape=(None,2), dtype = 'S100', chunks=True,  compression="gzip", compression_opts=9)
 
        channels = f.create_dataset('channels', (0, 0), maxshape=(None,None), dtype='f8', chunks=True,  compression="lzf")
        makers = f.create_dataset('markers', (0, 0), maxshape=(None,2), dtype = 'S100', chunks=True,  compression="lzf")
        

        metadata = {'id': bam_in
                   }
        f.attrs.update(metadata)
        
def add_pentas(h5_out, penta_arr, marker_arr):
    try:
        with h5py.File(h5_out, 'a') as f:
            channels = f['channels']
            ch_end = channels.shape[0]
            
            shp = channels.shape[0] + penta_arr.shape[0]
            channels.resize((shp, 5))    
            channels[ch_end:shp] = penta_arr 
            
            markers = f['markers']
            mrk_end = markers.shape[0]
            mrk_shp = markers.shape[0] + marker_arr.shape[0]
            markers.resize((mrk_shp, 2))
            markers[mrk_end:mrk_shp] = marker_arr
            
    except FileNotFoundError:
        print("{0} not found".format(h5_out))


def make_penta(seq_str, qual_str, depth_str, ref, qual_min = 50):
    seq_str = seq_str.upper()
   
    nuc_dict = {"A":0, "T":1, "G":2, "C":3, "*":4, "N":5}
    nucs = {"A": 0, "T": 0, "G":0, "C":0, "*":0, "N":0}
    inserts = {"A": 0, "T": 0, "G":0, "C":0, "*":0, "N":0} 
    
    if depth_str == "0" and seq_str == "*" and qual_str == "*":
        seq_channel = [[0,0,0,0,0], [0,0,0,0,0]]
        return seq_channel
    
    else:
        i = 0
        q = 0
        while i < len(seq_str):
            if seq_str[i] == "$": 
                i += 1
            elif seq_str[i] == "^": 
                i += 2
                
            elif seq_str[i] in [".", ","]:
                if ord(qual_str[q]) > qual_min:
                    nucs[ref] += 1
                i += 1
                q += 1

            elif seq_str[i] in nucs:
                if ord(qual_str[q]) > qual_min:
                    nucs[seq_str[i]] += 1
                i += 1
                q += 1

                
            elif seq_str[i] == "+":
                i += 1
                j = 0
                insert_str = ""
                while seq_str[i].isnumeric():
                    insert_str += seq_str[i]
                    i += 1
                    j += 1
                insert_int = int(insert_str)
                insert_seq = seq_str[i:i + insert_int]
                i += insert_int

                for s in range(len(insert_seq)):
                    inserts[insert_seq[s]] += 1

            elif seq_str[i] == "-":
                i += 1
                j = 0
                gap_str = ""
                while seq_str[i].isnumeric():
                    gap_str += seq_str[i]
                    i += 1
                    j += 1
                gap_int = int(gap_str)
                i += gap_int
                
        seq_channel = [list(nucs.values())[0:5], list(inserts.values())[0:5]]

        return seq_channel



#if no reference flag is passed, run mpileup without it
if args.reference:
    if args.positions_file:
        cmd = "samtools mpileup -A -aa -f {0} -l {1} {2}".format(args.reference, args.positions_file, args.bam_infile).split()
    else:
        cmd = "samtools mpileup -A -aa -f {0} {1}".format(args.reference, args.bam_infile).split()
else:
    if args.positions_file: 
        cmd = "samtools mpileup -A -aa -l {0} {1}".format(args.positions_file, args.bam_infile).split()
    else:
        cmd = "samtools mpileup -A -aa {1}".format(args.bam_infile).split()


penta_all = []
mrks_all = []

init_penta(args.bam_infile, h5_outfile)
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):

    pile_up = parse_line(line)

    c1 = make_penta(seq_str = pile_up["pop_bam"][1], 
                    qual_str = pile_up["pop_bam"][2], 
                    depth_str = pile_up["pop_bam"][0], 
                    ref = pile_up["ref"],
                    qual_min=args.qual_min
                    )

    penta_all.append(c1[0])
    penta_all.append(c1[1])
    chrom = pile_up["chrom"]
    pos = pile_up["pos"]
    mrk_init = [chrom, pos]
    mrk_2 = [chrom, "i"+pos]
    mrks_all.append(mrk_init)
    mrks_all.append(mrk_2)
    if len(penta_all) >= args.buffer:
        add_pentas(h5_outfile, np.array(penta_all), np.squeeze(np.array(mrks_all, dtype = "S100")))
        penta_all = []
        mrks_all = []

add_pentas(h5_outfile, np.array(penta_all), np.squeeze(np.array(mrks_all, dtype = "S100")))

