from Bio import SeqIO
from collections import defaultdict, Counter

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
    kmer_counts = defaultdict(int)
    kmers_of_interest = set()
    total_seqs = 0

    for h,i in enumerate(SeqIO.parse(fn, 'fasta')):
        d,s,L = str(i.description), str(i.seq).upper().replace('T','U'), len(i.seq)
        if 'eukaryota' in d.lower(): continue

        total_seqs += 1
        for k in kmers:
            klen = k[-1][1] + 1
            for i in range(0, len(s)+1-klen):
                subseq = s[i:i+klen]
                bp_found = [subseq[bp[1]] == bp[0] for bp in k]
                if min(bp_found) == True:
                    kmer_counts[k] += 1
                    break

        print(total_seqs,end='\r')
        if total_seqs == 100000: break

    for k,v in sorted(kmer_counts.items(),key=lambda x:x[1]):
        if v/total_seqs > threshold:
            kmers_of_interest.add(k)

    print()
    print(f"k-mers used for analysis: {len(kmers)}")
    print(f"Observed frequency threshold: {threshold}")
    print(f"Sequences analyzed: {total_seqs}")
    print(f"Total k-mers observed: {len(kmer_counts)}")
    print(f"Total k-mers possible: {4**klen}")
    print(f"Number of biomarker k-mers identified: {len(kmers_of_interest)}")
    
    return kmers_of_interest, kmer_counts

if __name__ == '__main__':
    km = get_univ_kmers()
    koi, kc = search_for_kmers(km)
    print(kc)