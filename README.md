# mapcount

The Perl program **mapcount** counts next generation sequencing reads that map against a reference gene set.

**mapcount** extracts a list of expected sequence IDs from a user-provided FASTA file (such as a file of gene models).  It then parses a user-provided mapping file in the standard SAM format and counts the number of reads that map to each gene.  Tabular output is written to a user-defined text file.

Optional flags allow the user to:

+ Trim the sequence ID before a space (*e.g.*, "comp18530_c0_seq1 len=1406 path=[683:0-1405]" would instead be represented as "comp18530_c0_seq1").
+ Define the minimum match length to include the read (*e.g.*, all matching reads smaller than some specified value, such as 25 base pairs, can be excluded from the output table).
+ Count paired end reads as only a single match. (Although there might be two paired reads, these represent only one underlying nucleotide fragment).


Program usage is as follows:

```
map_count -f|fasta fasta_file -s|sam SAM_file -o|output filename [-t|trim] [-l|length 25] [-p|pair]
```
```
Required flags:
-f|fasta    fasta reference file
-s|sam      SAM file
-o|output   output filename
```
```
Optional flags:
-t|trim     trim ID before a space
-l|length   minimum match length for counting
-p|pair     count paired matches as one match
```
