from Bio import SeqIO
import statistics

# parę słów wstępu
# wg publikacji z bakty "Protein sequences of unknown function: product denoted as hypothetical
# protein, putative protein, uncharacterized protein or conserved predicted protein." więc takie są
# założenia tutaj

def read_file(filename):
    with open(filename, 'rt') as file:
        records = [rec for rec in SeqIO.parse(file, 'genbank')]
    return records

def find_genes(records):
    return [f for r in records for f in r.features if f.type == 'CDS']

def find_function(genes):
    unknown = 0
    known = 0
    for feature in genes:
        product = feature.qualifiers.get('product')[0]
        if product in ('hypothetical protein', 'putative protein', 'uncharacterized protein', 'predicted protein'):
            unknown += 1
        else:
            known += 1
    return unknown, known

def median_len(genes):
    length = [len(gene) for gene in genes]
    return statistics.median(length)

def main():
    records = read_file('assembly.gbff')
    genes = find_genes(records)
    unknown, known = find_function(genes)
    print(f'liczba genow o nieznanej funkcji: {unknown}\nliczba genow o znanej funkcji: {known}')
    print(f'mediana dlugosci genu: {median_len(genes)}')

if __name__ == '__main__':
    main()
