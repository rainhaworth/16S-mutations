from Bio import SeqIO
import numpy as np
import time
import argparse
import os
import csv
from collections import defaultdict

from utils.find_gapped import find_gapped

# create gapped kmers from universal kmers
def get_univ_kmers(fn='./universal-residues-ecoli.txt', k=6, maxsize=100):
    with open(fn, 'r') as f:
        residues = f.readlines()
    
    residues = [(x[0], int(x[1:])) for x in residues]
    kmers = [residues[i:i+k] for i in range(len(residues)-k+1)]
    
    drop = []
    for i in range(len(kmers)):
        offset = kmers[i][0][1]
        if kmers[i][-1][1] - offset <= maxsize:
            kmers[i] = tuple((x[0],x[1]-offset) for x in kmers[i])
        else:
            drop.append(i)
    for i in reversed(drop):
        kmers.pop(i)
    return kmers

# search sequences for candidate universal kmers
def search_for_kmers(kmers, nmax=1000, threshold=0.5, fn='./SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta'):
    kmer_counts = np.zeros(len(kmers), dtype=int)
    kmers_of_interest = set()
    total_seqs = 0
    seqs_none_found = 0

    # make numpy arrays for input
    kmer_chars = np.array([[x[0] for x in kmer] for kmer in kmers])
    kmer_offs = np.array([[x[1] for x in kmer] for kmer in kmers])

    # output data structures
    seq_idx_to_kmers = [[]] * nmax
    kmer_sets_to_seqs = defaultdict(list)

    for i in SeqIO.parse(fn, 'fasta'):
        s = str(i.seq).upper().replace('T','U')
        #if 'eukaryota' in d.lower(): continue

        # find kmers
        found = find_gapped(s, kmer_chars, kmer_offs)
        kmer_counts += found
        if np.sum(found) == 0:
            print('none found:', total_seqs, d)
            seqs_none_found += 1

        # populate output data structures
        kmer_hits = np.argwhere(found)[:,0]
        seq_idx_to_kmers[total_seqs] = kmer_hits
        kmer_sets_to_seqs[tuple(kmer_hits)].append(total_seqs)

        # iterate
        total_seqs += 1
        #print(total_seqs,end='\r')
        if total_seqs == nmax: break

    for i,v in enumerate(kmer_counts):
        if v/total_seqs > threshold:
            kmers_of_interest.add(kmers[i])

    print()
    print(f"k-mers used for analysis: {len(kmers)}")
    print(f"Observed frequency threshold: {threshold}")
    print(f"Sequences analyzed: {total_seqs}")
    print(f"Total k-mers found: {sum(kmer_counts)}")
    print(f"Number of biomarker k-mers identified: {len(kmers_of_interest)}")
    print(f"Sequences containing no kmers: {seqs_none_found}")
    
    return kmer_counts, kmer_counts / total_seqs, seq_idx_to_kmers, kmer_sets_to_seqs

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-k', type=int, default=8)
    p.add_argument('-n', type=int, default=1000)
    p.add_argument('-t', type=float, default=0.9)
    p.add_argument('-i', type=str, default='./SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta')
    p.add_argument('-o', type=str, default='./results')
    p.add_argument('--save', action='store_true')
    args = p.parse_args()

    st = time.time()
    km = get_univ_kmers(k=args.k)
    kc, kf, s2k, k2s = search_for_kmers(km, args.n, args.t, args.i)
    best = np.argmax(kf)

    print('runtime:', time.time() - st)
    
    if args.save:
        print('saving to', args.o)
        
        fn_end = '-' + str(args.k) + '-' + str(args.n) + '.csv'
        fn_k = os.path.join(args.o, 'k' + fn_end)
        fn_s2k = os.path.join(args.o, 's2k' + fn_end)
        fn_k2s = os.path.join(args.o, 'k2s' + fn_end)

        try: os.mkdir(args.o)
        except FileExistsError: pass

        
        with open(fn_k, 'w') as f:
            cw = csv.writer(f, delimiter=';')
            out = [[str(km[i]), kc[i]] for i in range(len(kc))]
            cw.writerows(out)
        
        with open(fn_s2k, 'w') as f:
            cw = csv.writer(f)
            cw.writerows(s2k)

        with open(fn_k2s, 'w') as f:
            cw = csv.writer(f, delimiter=';')
            for k,v in k2s.items():
                cw.writerow([str(k)] + v)
            

    print('freqs:', kf)
    print('best:', kf[best], km[best])
