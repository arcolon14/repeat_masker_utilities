#!/usr/bin/env python3
import sys, os, argparse, datetime, gzip

#
# Globals
#
DATE = datetime.datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
DESC = 'Parse and merge the output from RepeatMasker.'

def parse_args():
    '''Command line arguments'''
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-c', '--cross-match', required=True,
                   help='(str) Path to the cross match (*.out) table from Repeat Masker.')
    p.add_argument('-d', '--divsum', required=True,
                   help='(str) Path to the divergence summary (*.divsum) table from Repeat Masker.')
    p.add_argument('-o', '--outdir', required=False, default='.',
                   help='(str) Path to output directory.')
    p.add_argument('-b', '--basename', required=False, default=f'RepeatMasker_{DATE}',
                   help='(str) Basename for the output files  [default=RepeatMasker_YYYYMMDD]')
    p.add_argument('-m', '--min-length', required=False, default=None,
                   help='(int) Minimum length of well characterized bases (wellCharLen) needed to keep a sequence  [default=None]')
    # Check input arguments
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: '{args.outdir}' not found.")
    if not os.path.exists(args.cross_match):
        sys.exit(f"Error: '{args.cross_match}' not found.")
    if not os.path.exists(args.divsum):
        sys.exit(f"Error: '{args.divsum}' not found.")
    if args.min_length is not None:
        # If not None, a numeric value must be provided
        min_len = str(args.min_length)
        if min_len.isnumeric():
            args.min_length = int(min_len)
        else:
            sys.exit(f'"Error: min len ({min_len}) must be numeric.')
    else:
        # If set to default None, set to a very small negative number
        args.min_length = -1_000_000_000
    return args

def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

class RepeatDivsum:
    '''Divergence summary statistics for a given repeat.'''
    def __init__(self, rep_class, rep_id, abs_len, well_char_len, kimura):
        self.r_class   = rep_class.split('/')[0]
        self.r_family  = rep_class
        self.r_id      = rep_id
        self.abs_len   = abs_len
        self.wchar_len = well_char_len
        self.kimura    = kimura/100
    def __str__(self):
        return f'{self.r_class} {self.r_family} {self.r_id} {self.abs_len} {self.wchar_len} {self.kimura:0.06g}'

class RepeatAnnot:
    '''Annotation metadata for a given repeat.'''
    def __init__(self, rep_name, rep_class, chromosome, start, end, strand, sw_score):
        self.name     = rep_name
        self.r_class  = rep_class.split('/')[0]
        self.r_family = rep_class
        self.chrom    = chromosome
        self.start    = start
        self.end      = end
        self.strand   = strand
        self.score    = sw_score
    def __str__(self):
        return f'{self.name} {self.r_class} {self.r_family} {self.chrom} {self.start} {self.end} {self.score}'

def parse_divsum_table(divsum_f, min_len=0):
    '''Parse the divergence summary table.'''
    fh = open(divsum_f, 'r')
    if divsum_f.endswith('.gz'):
        fh = gzip.open(divsum_f, 'rt')
    print('\nParsing divergence summary table...')
    divsum = dict()
    records = 0
    discard = 0
    # Parse file
    for i, line in enumerate(fh):
        line = line.strip('\n')
        # Skip comments and empty lines
        if line.startswith('#'):
            continue
        if len(line) == 0:
            continue
        # Split into different fields.
        # The target fields have 5 entries
        fields = line.split('\t')
        if len(fields) != 5:
            continue
        # Skip headers or placeholders
        if fields[0] in {'Class', 'ARTEFACT', 'Simple_repeat', '-----'}:
            continue
        # Work with the desired records
        records += 1
        repeat_summary = RepeatDivsum(rep_class=fields[0],
                                       rep_id=fields[1],
                                       abs_len=int(fields[2]),
                                       well_char_len=int(fields[3]),
                                       kimura=float(fields[4]))
        # Remove entries if too small
        if repeat_summary.wchar_len < min_len:
            discard += 1
            continue
        # save for export
        divsum[repeat_summary.r_id] = repeat_summary
    print(f'    Extracted {records:,} records from divergence summary table file.\n    Discarded {discard:,} records.')
    return divsum

def parse_crossmatch_table(crossmatch_f):
    '''Parse the cross_match table output from RepeatMasker'''
    fh = open(crossmatch_f, 'r')
    if crossmatch_f.endswith('.gz'):
        fh = gzip.open(crossmatch_f, 'rt')
    print('Parsing cross match table...')
    cross_match = list()
    records = 0
    kept = 0
    # Parse file
    for i, line in enumerate(fh):
        line = line.strip('\n')
        # Skip comments and empty lines
        if line.startswith('#'):
            continue
        if len(line) == 0:
            continue
        # Split into different fields, remove headers
        fields = line.split()
        if not fields[0].isnumeric():
            continue
        # Process the record line if possible
        try:
            # Now, process the record lines
            records += 1
            sw_score   = int(fields[0])
            chromosome = fields[4]
            start_bp   = int(fields[5])
            end_bp     = int(fields[6])
            re_name    = fields[9]
            re_class   = fields[10]
            strand     = fields[8]
            if strand == 'C':
                strand = '-'
        except IndexError:
            print(f"Error: line {i} of {crossmatch_f} is not following the standard format:\n\n{line}")
        annotation = RepeatAnnot(re_name, re_class, chromosome, start_bp, end_bp, strand, sw_score)
        # Add to the annotation dict
        kept += 1
        cross_match.append(annotation)
    print(f'    Read {records:,} records from the cross_match table file.\n    Retained {kept:,} records.')
    return cross_match

def merge_cross_divsum(cross_match, divsum, outdir='.', basename=f'RepeatMasker_{DATE}'):
    '''Merge the cross_match and divsum datasets into a final output.'''
    print('\nMatching cross_match and divsum records...')
    fh = open(f'{outdir}/{basename}.repeat_masked_merged.tsv', 'w')
    header = ['#Chromosome', 'StartBP', 'EndBP', 'Strand', 'Name', 'Class', 'Family', 'WellCharLen', 'Kimura','SwScore']
    header = '\t'.join(header)
    fh.write(f'{header}\n')
    matches = 0
    # Process each cross_match annotation
    for annotation in sorted(cross_match, key=lambda a: a.start):
        assert isinstance(annotation, RepeatAnnot)
        wchar_len = None
        kimura = None
        # Find the divsum for the annotation, if available
        divergence = divsum.get(annotation.name, None)
        if divergence is not None:
            assert isinstance(divergence, RepeatDivsum)
            # Make sure the class and family match exactly
            assert annotation.r_class == divergence.r_class, f'{annotation.r_class} {divergence.r_class}'
            # assert annotation.r_family == divergence.r_family, f'{annotation.r_family} {divergence.r_family}'
            # Add the given values
            wchar_len = divergence.wchar_len
            kimura = f'{divergence.kimura:0.6g}'
        # Process only fully matching records
        matches += 1
        ['#Chromosome', 'StartBP', 'EndBP', 'Strand', 'Name', 'Class', 'Family', 'WellCharLen', 'Kimura','SwScore']
        row = f'{annotation.chrom}\t{annotation.start}\t{annotation.end}\t{annotation.strand}\t{annotation.name}\t{annotation.r_class}\t{annotation.r_family}\t{wchar_len}\t{kimura}\t{annotation.score:0.6g}\n'
        fh.write(row)
    print(f'    Exported a total of {matches:,} matching records.')
    fh.close()

def clean_divsum_file(divsum, outdir, basename):
    '''Print a cleaner DIVSUM output file.'''
    print('\nPrinting clean divum file...')
    fh = open(f'{outdir}/{basename}.divsum.tsv', 'w')
    header = 'Class\tRepeat\tabsLen\twellCharLen\tKimura\n'
    fh.write(header)
    for element in sorted(divsum):
        div = divsum[element]
        assert isinstance(div, RepeatDivsum)
        row = f'{div.r_family}\t{div.r_id}\t{div.abs_len}\t{div.wchar_len}\t{div.kimura:0.08g}\n'
        fh.write(row)
    fh.close()

def main():
    print(f'{PROG} started on {now()}\n')
    args = parse_args()
    # Parse cross match table
    cross_match = parse_crossmatch_table(args.cross_match)
    # Parse divsum table
    divsum = parse_divsum_table(args.divsum, args.min_length)
    # Join cross_match and divsum
    merge_cross_divsum(cross_match, divsum, args.outdir, args.basename)
    # Print clean divsum
    clean_divsum_file(divsum, args.outdir, args.basename)
    # Finish
    print(f'\n{PROG} finished on {now()}')


# Run Code
if __name__ == '__main__':
    main()
