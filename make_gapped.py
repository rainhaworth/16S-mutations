# create gapped kmers from universal kmers
def get_univ_kmers(fn = './universal-residues-ecoli.txt', k=6, maxsize=100):
    with open(fn, 'r') as f:
        residues = f.readlines()
    
    residues = [(x[0], int(x[1:])) for x in residues]
    kmers = [residues[i:i+k] for i in range(len(residues)-k+1)]
    
    drop = []
    for i in range(len(kmers)):
        offset = kmers[i][0][1]
        if kmers[i][-1][1] - offset <= maxsize:
            kmers[i] = [(x[0],x[1]-offset) for x in kmers[i]]
        else:
            drop.append(i)
    for i in reversed(drop):
        kmers.pop(i)
    return kmers

if __name__ == '__main__':
    km = get_univ_kmers()
    print(len(km))
    print(km)
