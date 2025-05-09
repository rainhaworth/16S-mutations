from Bio import SeqIO

f = 'SILVA_138.2_SSURef_NR99_tax_silva.fasta'
f_metadata = set()

for line in open('SILVA_138.2_SSURef_Nr99.full_metadata', 'r'):
  elems = line.strip().split('\t')
  try:
    assert 'metagenome' not in elems[10]
    assert float(elems[37]) == 100.0
    f_metadata.add('.'.join(elems[:3]))
  except:
    continue

print(len(f_metadata), 'IDs from metadata')

seqs_before_filtering = 0
seqs_after_filtering = 0

with open('SILVA_138.2_SSURef_NR99_tax_silva_filtered_strict_2.fasta','w') as out:
  for i in SeqIO.parse(f,'fasta'):
    seqs_before_filtering += 1
    id,d,s = str(i.id),str(i.description),str(i.seq).upper().replace('T','U')
    if id not in f_metadata: # filter out sequences not labeled as 16S rRNA in SILVA metadata
      continue
    elif 'eukaryota' in d.lower(): # filter out eukaryotic sequences e.g., 18S rRNA or mitochondria/plastid 16S rRNA
      continue
    elif set(s) != set('ACGU'): # filter out sequences with low quality base calls
      continue
    else:
      out.write(f'>{d}\n{s}\n')
      seqs_after_filtering += 1

print(f'Number of seqs before filtering: {seqs_before_filtering}')
print(f'Number of seqs after filtering: {seqs_after_filtering}')
