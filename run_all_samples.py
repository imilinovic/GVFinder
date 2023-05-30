import os
import subprocess

EXAMPLES_DIR = 'data/'

def main():
    for dir in os.listdir(EXAMPLES_DIR):
        if not dir.startswith('J'):
            continue
        
        filename = EXAMPLES_DIR + dir
        print(filename)
        subprocess.run(["build/src/GVFinder"], input=bytes(' '.join([filename, "jelen_summary.txt"]), 'utf-8'))
        
if __name__ == '__main__':
    main()