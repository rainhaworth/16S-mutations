from Bio import SeqIO
import numpy as np
import time
import argparse
import os

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

    # make numpy arrays
    kmer_chars = np.array([[x[0] for x in kmer] for kmer in kmers])
    kmer_offs = np.array([[x[1] for x in kmer] for kmer in kmers])

    for i in SeqIO.parse(fn, 'fasta'):
        d, s = str(i.description), str(i.seq).upper().replace('T','U')
        if 'eukaryota' in d.lower(): continue

        total_seqs += 1
        found = find_gapped(s, kmer_chars, kmer_offs)
        kmer_counts += found
        if np.sum(found) == 0:
            print('none found:', d)
            seqs_none_found += 1

        print(total_seqs,end='\r')
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
    
    return kmers_of_interest, kmer_counts, kmer_counts / total_seqs

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-k', type=int, default=8)
    p.add_argument('-n', type=int, default=1000)
    p.add_argument('-t', type=float, default=0.9)
    p.add_argument('-i', type=str, default='./SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta')
    p.add_argument('-o', type=str, default='./')
    args = p.parse_args()

    st = time.time()
    km = get_univ_kmers(k=args.k)
    koi, kc, kf = search_for_kmers(km, args.n, args.t, args.i)
    best = np.argmax(kf)

    print(km)
    print(kf)
    np.save(os.path.join(args.o, 'kf-'+str(args.k)+'-'+str(args.n)+'.o'), kf)

    print('best:', kf[best], km[best])
    print('time:', time.time() - st)