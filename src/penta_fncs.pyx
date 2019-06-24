import numpy as np


def make_penta(seq_str, qual_str, depth_str, ref, qual_min = 50):
    seq_str = seq_str.upper()
    
    nuc_dict = {"A":0, "T":1, "G":2, "C":3, "*":4}
    nucs = {"A": 0, "T": 0, "G":0, "C":0, "*":0}
    inserts = {"A": 0, "T": 0, "G":0, "C":0, "*":0}
    
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
                
        seq_channel = [list(nucs.values()), list(inserts.values())]
        return seq_channel

