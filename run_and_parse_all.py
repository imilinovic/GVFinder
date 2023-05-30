import os
import subprocess

from Levenshtein import distance as levenshtein_distance

EXAMPLES_DIR = 'data/'

def main():
    try:
        os.remove("tmp.txt")
    except OSError:
        pass

    for dir in os.listdir(EXAMPLES_DIR):
        if not dir.startswith('J'):
            continue
        
        filename = EXAMPLES_DIR + dir
        print(filename)
        subprocess.run(["build/src/GVFinder"], input=bytes(' '.join([filename, "tmp.txt"]), 'utf-8'))
        
    with open("tmp.txt", 'r') as f:
        data = f.readlines()

    occurences = {}
    most_common = []

    for line in data:
        striped_line = line.strip()
        if striped_line.startswith('data'):
            filename = striped_line.split('/')[-1]
        elif len(striped_line) > 200: # ignoriranje jako losih sampleova
            if striped_line not in occurences:
                occurences[striped_line] = [filename]
            else:
                occurences[striped_line].append(filename)

    for genome, occurence in occurences.items():
        most_common.append((len(occurence), genome))
    most_common = sorted(most_common, reverse=True)

    with open("summary_analized.txt", 'w') as f:
        for cnt, genome in most_common:
            f.write(f"Sekvenca: {genome}\n")
            f.write(f"Broj pojavljivanja: {cnt}\n")
            f.write(f"Pronađena u podacima: {', '.join(occurences[genome])}\n\n")
    
    with open("summary_levenshtein.txt", 'w') as f:
        f.write("Tablica Levenshteinovih udaljenosti između pronađenih sekvenci (i-ta sekvenca u tablici je i-ta najčešća sekvenca iz datoteke jelen_summary_analized.txt):\n")
        f.write("Ako je udaljenost veća ili jednaka 10 je prikazana kao -\n")
        for _, genome1 in most_common:
            f.write(' '.join([str(levenshtein_distance(genome1, genome2)) if levenshtein_distance(genome1, genome2) < 10 else '-' for _, genome2 in most_common]) + '\n')

if __name__ == '__main__':
    main()