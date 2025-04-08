from Bio import SeqIO

f = 'SILVA_138.2_SSURef_NR99_tax_silva.fasta'
seqs_before_filtering = 0
seqs_after_filtering = 0

with open('SILVA_138.2_SSURef_NR99_tax_silva_filtered.fasta','w') as out:
  for i in SeqIO.parse(f,'fasta'):
    seqs_before_filtering += 1
    d,s = str(i.description),str(i.seq).upper().replace('T','U')
    if 'eukaryota' in d.lower():
      continue
    elif set(s) != set('ACGU'):
      continue
    else:
      out.write(f'>{d}\n{s}\n')
      seqs_after_filtering += 1

print(f'Number of seqs before filtering: {seqs_before_filtering}')
print(f'Number of seqs after filtering: {seqs_after_filtering}')
  
