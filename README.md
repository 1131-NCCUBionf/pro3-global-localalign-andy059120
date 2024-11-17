[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/pWmxMLzQ)
# pro3.global|local-Aln
* name : 陳昶安 Chen,Chang-An
* student ID : 113753121

## Description

* Write a Python script to perform a global or local alignment.
* Creating your own script, i.e. hw3.py.
* In this program, library Biostrings is only used to parse input fasta file.
* Packages you can use: numpy, pandas, Bio
* You should write a program with a function named alignment, ie.
```
def alignment(input_path, score_path, output_path, aln, gap):
    .
    .
    .
    .
```
* If there is more than one local alignment with the same highest score, you should output local alignments with the maximum length. 
* If there is more than one local alignment with the same highest score, you should output those local alignments in string sequential order according to protein1 and then protein2, i.e., 
  ```
  >protein1
  local alignment1
  >protein2
  local alignment1
  >protein1
  local alignment2
  >protein2
  local alignment2
  ```
## Parameters

* input: .fasta file (ex. test_global.fasta)
* score: score file (ex. pam250.txt)
* aln: global|local
* gap: gap score
* output: .fasta file

## Files

* hw3_ref.py: You can start from this reference code, and try to write your own comment in English.
* pam100.txt
* pam250.txt
* test_global.fasta
* result_global.fasta: You should output your alignment in FASTA format.
* test_local.fasta
* result_local.fasta
## Command

Executing your code with the following command.


```Python
alignment("examples/test_global.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10)
alignment("examples/test_local.fasta", "examples/pam100.txt", "examples/result_local.fasta", "local", -10)
```

## Evaluation

10 testing data(5 public, 5 private)

The correct answer gets 10 points for each testing data.



### Penalty

* High code similarity to others: YOUR SCORE = 0

## References

* no reference





