# privateHLA

## Introduction

In recent years, several HLA imputation methods on local servers have been developed, using single nucleotide polymorphisms (SNPs) as inputs to predict Human Leukocyte Antigen (HLA) genotypes. However, these methods require HLA reference panels which are often unavailable and memory-intensive. Cloud-based outsourced HLA imputation overcomes these limitations by utilizing built-in reference panels on the cloud server, removing local storage needs. However, uploading genotype data online raises privacy concerns. Although secure from third parties, the uploaded data can be misused by cloud server administrators.  Additionally, reference panels and their HLA imputation model on the cloud servers must be safeguarded against malicious clients. To address these privacy issues, we developed privateHLA, a secure HLA imputation method using homomorphic encryption. privateHLA securely performs HLA imputation on an outsourced server, protecting both client’s data and enhances model privacy. privateHLA outperformed SNP2HLA but had slightly lower accuracy than CookHLA, both plaintext-based methods.

## Installation

This repository contains Go code for privateHLA. One can run it as any other Go program in all platforms supported by Go runtime, including Windows, macOS and Linux on x86 and ARM.

First, install Go, following the official [instruction](https://go.dev/doc/install).

For Mac users, using package managers like homebrew might be perferrable:
```
brew install go
```

For Linux users, Go is likely to be available out-of-the-box in your favourite package manager.

Finally, clone this repository and run `cmd/main.go` as follows.
```
go run cmd/main.go -dataset="HapMap_EUR" -prefix="A"
```
The `main.go` file executes the main components of `privateHLA`, including:

- Encryption and decryption
- Evaluation of the following functions:
  - Linear
  - Top2
  - Sigmoid
  - Thresholding

Currently `dataset` flag accepts `HapMap_EUR, Korean`, and `prefix` flag accepts `A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1`.
To automatically run all dataset and prefix, run `out.sh`.
```
sh cmd/out.sh
```
Benchmarks using Go benchmark is also available. Run
```
go test . -run=^$ -bench=.
```

### Dataset Structure
This repository also provides toy example datasets, each consisting of 9 individuals from the HapMap and Korean populations.
Input data must be placed inside the specified dataset folder (e.g., `HapMap_EUR`).
Both the genotype and the corresponding weights folders must be included in the same folder.

Example:
```
├── HapMap_EUR/
│   ├── X/
│   ├── weights/
```

### Output to File
To execute and save the output to a text file, run:

```bash
go run cmd/main.go -dataset="HapMap_EUR" -prefix="A" > output.txt
```

This will process the dataset located in the `HapMap_EUR` folder, using the target allele prefix `"A"`, and save the output to `output.txt`.
