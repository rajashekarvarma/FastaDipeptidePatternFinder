# FastaDipeptidePatternFinder

## Scan and find di peptide patterns (P-S/T/Y or S/T/Y-P) from Fasta protein sequences. FastaPatternFinder is written in python3 and requires python3 to run.

### prerequisites

     Python modules required
        re
        Pandas
        BioPython
        argparse
        
### Usage

#### python FastaPatternFinder.py -i <input_filename.fa> -o <output_filename.csv>
Defining output file(currently only csv is supported) is an option if the output file name is not provided the output will be written to output.csv in the working directory.

#### python FastaPatternFinder.py -h           ### for help
