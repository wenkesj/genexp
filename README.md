# genexp

[![GoDoc](https://godoc.org/github.com/wenkesj/genexp?status.svg)](https://godoc.org/github.com/wenkesj/genexp)

Package `genexp` includes an API for gene expression analysis using methods I
learn over the bioinformatics course. It comes with MT transcripts and
chromosome for the Human genome. Everything is a WIP, don't expect to find
anything amazing.

```sh
go get github.com/wenkesj/genexp
```

**NOTE**: Download and install <https://golang.org/doc/install>.

## genexpp

Package `genexpp` is a CLI gene predictor and coding potential calculator.

This algorithm uses a simple strategy:
+ Find all longest possible ORFs of a given length threshold.
+ Assign a coding potential based on the usage of given codon triplets. The
  default generates the codon usage table for human MT.
+ The coding potential is a product of probabilities of finding possible codons
  for the given threshold.

This project also includes an API friendly in which the user can change the
parameters of the codon profiling to a specific list codons to get higher coding
potential scores.

### Get Started with `genexpp`

- Install the CLI for this project, run this from the command-line:
  ```sh
  $ go get github.com/wenkesj/genexp/genexpp
  ```

- Run the program:
  ```sh
  $ genexpp -transcripts=./transcripts \
    -genome=./genomes/Homo_sapiens.GRCh38.dna.chromosome.MT.fa \
    -orflengths=50,2000
  ```

- Expected output:
  ```
  genexpp -transcripts=./transcripts \
  >     -genome=./genomes/Homo_sapiens.GRCh38.dna.chromosome.MT.fa \
  >     -orflengths=50,2000
  Matches: 15598
  Matching: 99.99%
  > Homo_sapiens.GRCh38.dna.chromosome.MT.fa
  +--------------------------------+----------------------+
  | TRANSCRIPT (DISTANCE, LOWER IS BETTER) | CODING POTENTIAL (%) |
  +--------------------------------+----------------------+
  | Q9NPG2.fasta (250)             | 2.17E-80%            |
  | ENSG00000212907.txt (219)      | 1.84E-89%            |
  | ENSG00000198840.txt (265)      | 4.51E-90%            |
  | ENSG00000198899.txt (465)      | 4.02E-91%            |
  | ENSG00000198899.txt (470)      | 1.22E-90%            |
  ...
  | ENSG00000212907.txt (204)      | 6.30E-83%            |
  | ENSG00000198840.txt (248)      | 3.99E-88%            |
  | ENSG00000198840.txt (250)      | 4.02E-88%            |
  | Q9NPG2.fasta (251)             | 2.78E-78%            |
  | Q9NPG2.fasta (247)             | 1.10E-76%            |
  | Q9NPG2.fasta (254)             | 2.81E-76%            |
  +--------------------------------+----------------------+
  ```

## genexpa

Package `genexpa` is a CLI sequence aligner using local alignment method.

### Get Started with `genexpa`
- Install the CLI for this project, run this from the command-line:
  ```sh
  $ go get github.com/wenkesj/genexp/genexpa
  ```

- Run the program:
  ```sh
  $ genexpa ATCTGTTCTT CTTTATGTT
  ```

- Expected output (DP Table, alignments, etc.):
  ```
  $ genexpa ATCTGTTCTT CTTTATGTT
    C T T T A T G T T
  0 0 0 0 0 0 0 0 0 0
  A 0 0 0 0 0 0 0 0 0 0
  T 0 0 1 1 1 0 1 0 1 1
  C 0 0 0 0 0 0 0 0 0 0
  T 0 0 1 1 1 0 1 0 1 1
  G 0 0 0 0 0 0 0 2 1 0
  T 0 0 1 1 1 0 1 1 3 2
  T 0 0 1 2 2 1 1 0 2 4
  C 0 0 0 1 1 1 0 0 1 3
  T 0 0 1 1 2 1 2 1 1 2
  T 0 0 1 2 2 1 2 1 2 2
  Alignments:
  A:  ...__TCTGTT...
  B:  ...CTTATGTT...
  Score:  4
  ```
