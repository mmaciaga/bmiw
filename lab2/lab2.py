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
    light = read_file('assembly.gbff')
    genes_l = find_genes(light)
    unknown_l, known_l = find_function(genes_l)
    print('-' * 5, 'LIGHT', '-' * 5)
    print(f'liczba genow o nieznanej funkcji: {unknown_l}\nliczba genow o znanej funkcji: {known_l}')
    print(f'mediana dlugosci genu: {median_len(genes_l)}')

    full = read_file('assembly_full.gbff')
    genes_f = find_genes(full)
    unknown_f, known_f = find_function(genes_f)
    print('-' * 5, 'FULL', '-' * 5)
    print(f'liczba genow o nieznanej funkcji: {unknown_f}\nliczba genow o znanej funkcji: {known_f}')
    print(f'mediana dlugosci genu: {median_len(genes_f)}')

if __name__ == '__main__':
    main()
