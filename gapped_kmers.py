from Bio import SeqIO
from collections import defaultdict
import numpy as np
from utils.find_gapped import find_gapped
import time

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
def search_for_kmers(kmers, threshold=0.5, fn='./SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta'):
    kmer_counts = np.zeros(len(kmers), dtype=int)
    kmers_of_interest = set()
    total_seqs = 0

    # make numpy arrays
    kmer_chars = np.array([[x[0] for x in kmer] for kmer in kmers])
    kmer_offs = np.array([[x[1] for x in kmer] for kmer in kmers])

    for h,i in enumerate(SeqIO.parse(fn, 'fasta')):
        d,s,L = str(i.description), str(i.seq).upper().replace('T','U'), len(i.seq)
        if 'eukaryota' in d.lower(): continue

        total_seqs += 1
        kmer_counts += find_gapped(s, kmer_chars, kmer_offs)

        print(total_seqs,end='\r')
        if total_seqs == 100: break

    for i,v in enumerate(kmer_counts):
        if v/total_seqs > threshold:
            kmers_of_interest.add(kmers[i])

    print()
    print(f"k-mers used for analysis: {len(kmers)}")
    print(f"Observed frequency threshold: {threshold}")
    print(f"Sequences analyzed: {total_seqs}")
    print(f"Total k-mers found: {sum(kmer_counts)}")
    print(f"Number of biomarker k-mers identified: {len(kmers_of_interest)}")
    
    return kmers_of_interest, kmer_counts, kmer_counts / total_seqs

if __name__ == '__main__':
    st = time.time()
    km = get_univ_kmers(k=8)
    koi, kc, kf = search_for_kmers(km)
    print(kf)
    print('time:', time.time() - st)