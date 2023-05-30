# GVFinder

Genome variant finder from data from sequencing

## Building

Clone the spoa library inside of the root of GVFinder by following CMake (3.12+) instructions from <https://github.com/rvaser/spoa>

Create */data* folder in root and put all your samples there \
Build with following commands:
> *mkdir build* \
 *cd build* \
 *cmake -DCMAKE_BUILD_TYPE=Release ..* \
 *make*

## Running

Run the program with `./build/src/GVFinder` \
You will then be prompted to enter the path (from root) to the sequences you want to find the variant from (in .fastq format) and the output file. \

### Example

>Unesite datoteku s podacima: data/J29_B_CE_IonXpress_005.fastq \
Unesite datoteku za ispis rezultata: output.txt

Variants in the output file are sorted by cluster size (descending). The variants will also be printed to standard output (stdout). \
*tests/* were used to test matching with known samples 29 and 30 but are currently deprecated.

## Summary generation

Summaries were generated with *run_and_parse_all.py*. Run the script without any arguments. By default the script looks for all files in *data/* directory starting with the letter J. It generates files *summary_analizator.txt* and *summary_levenshtein.txt*.
