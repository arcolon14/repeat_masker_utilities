#!/usr/bin/env python3
import sys, os, argparse, datetime, gzip

DATE = datetime.datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
DESC = 'Extract repeat proportions from the RepeatMasker output table.'

def parse_args():
    '''Command line arguments'''
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-t', '--rm-table', required=True,
                   help='(str) Path to the table (*.tbl) output from Repeat Masker.')
    p.add_argument('-o', '--outdir', required=False, default='.',
                   help='(str) Path to output directory.')
    p.add_argument('-b', '--basename', required=False, default=f'RepeatMasker_{DATE}',
                   help='(str) Basename for the output files  [default=RepeatMasker_YYYYMMDD]')
    # Check input arguments
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: '{args.outdir}' not found.")
    if not os.path.exists(args.rm_table):
        sys.exit(f"Error: '{args.rm_table}' not found.")
    return args

def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

class RepProp:
    '''Store the repeat proportions for a given sequence.'''
    def __init__(self, seq_id):
        self.seq_id       = seq_id
        self.total_length = 0
        self.total_masked = 0
        self.interspersed = 0
        self.dna          = 0
        self.ltr          = 0
        self.line         = 0
        self.sine         = 0
        self.small_rna    = 0
        self.unclassified = 0
    def __str__(self):
        row = f'''Sequence ID:   {self.seq_id}
Total length:  {self.total_length:,} bp
Total masked:  {self.total_masked:,} bp
Interspersed:  {self.interspersed:,} bp
DNA:           {self.dna:,} bp
LTR:           {self.ltr:,} bp
LINE:          {self.line:,} bp
SINE:          {self.sine:,} bp
Small RNA:     {self.small_rna:,} bp
Unclassified:  {self.unclassified:,} bp'''
        return row

def parse_rm_table(rm_table, sequence_id):
    '''Parse the repeat masker output table.'''
    rep_props = RepProp(sequence_id)
    # Open the table file
    with gzip.open(rm_table, 'rt') if rm_table.endswith('.gz') else open(rm_table) as fh:
        for line in fh:
            line = line.strip('\n')
            # Skip empty line and table borders
            if len(line) == 0:
                continue
            if line[0] in {'-', '='}:
                continue
            line = line.lstrip(' ')
            fields = line.split()
            # For the total length
            if line.startswith('total length'):
                rep_props.total_length = int(fields[2])
            # For the masked bases
            elif line.startswith('bases masked'):
                rep_props.total_masked = int(fields[2])
            # Interspersed
            elif line.startswith('Total interspersed'):
                rep_props.interspersed = int(fields[3])
            # DNA elements
            elif line.startswith('DNA transposons'):
                rep_props.dna = int(fields[3])
            # LTRs
            elif line.startswith('LTR elements'):
                rep_props.ltr = int(fields[3])
            # LINEs
            elif line.startswith('LINEs'):
                rep_props.line = int(fields[2])
            # SINEs
            elif line.startswith('SINEs'):
                rep_props.sine = int(fields[2])
            # Small RNA
            elif line.startswith('Small RNA'):
                rep_props.small_rna = int(fields[3])
            # Unclassified
            elif line.startswith('Unclassified'):
                rep_props.unclassified = int(fields[2])
            else:
                continue
    return rep_props


def print_table(rep_props, basename, outdir='.'):
    '''Print summary table'''
    assert isinstance(rep_props, RepProp)
    out_f = f'{outdir}/{basename}.repeat_proportions.tsv'
    with open(out_f, 'w') as fh:
        # Prepare Header
        header = ['SequenceID', 'TotalSeqLen_BP', 'Masked_BP', 'Masked_Prop',
                  'Interspersed_BP', 'Interspersed_Prop', 'DNA_BP', 'DNA_Prop',
                  'LTR_BP', 'LTR_Prop', 'LINE_BP', 'LINE_Prop', 'SINE_BP', 
                  'SINE_Prop', 'smallRNA_BP', 'smallRNA_Prop', 'Unclassified_BP', 
                  'Unclassified_Prop']
        header = '\t'.join(header)
        fh.write(f'{header}\n')
        # Prepare data row
        total = rep_props.total_length
        row = [rep_props.seq_id,                             # Sequence ID
               f'{total}',                                   # Total length
               f'{rep_props.total_masked}',                  # Masked total
               f'{(rep_props.total_masked/total):0.08g}',    # Masked proportion
               f'{rep_props.interspersed}',                  # Interspersed total
               f'{(rep_props.interspersed/total):0.08g}',    # Interspersed proportion
               f'{rep_props.dna}',                           # DNA total
               f'{(rep_props.dna/total):0.08g}',             # DNA proportion
               f'{rep_props.ltr}',                           # LTR total
               f'{(rep_props.ltr/total):0.08g}',             # LTR proportion
               f'{rep_props.line}',                          # LINE total
               f'{(rep_props.line/total):0.08g}',            # LINE proportion
               f'{rep_props.sine}',                          # SINE total
               f'{(rep_props.sine/total):0.08g}',            # SINE proportion
               f'{rep_props.small_rna}',                     # LINE total
               f'{(rep_props.small_rna/total):0.08g}',       # LINE proportion
               f'{rep_props.unclassified}',                  # LINE total
               f'{(rep_props.unclassified/total):0.08g}']    # LINE proportion
        row = '\t'.join(row)
        fh.write(f'{row}\n')
    # Report to log
    print(rep_props)

def main():
    print(f'{PROG} started on {now()}\n')
    args = parse_args()
    # Parse the repeat table
    repeat_proportions = parse_rm_table(args.rm_table, args.basename)
    # Generate output
    print_table(repeat_proportions, args.basename, args.outdir)
    # Finish
    print(f'\n{PROG} finished on {now()}')

# Run Code
if __name__ == '__main__':
    main()
