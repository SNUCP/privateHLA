# privateHLA

## Introduction

In recent years, several HLA imputation methods on local servers have been developed, using single nucleotide polymorphisms (SNPs) as inputs to predict Human Leukocyte Antigen (HLA) genotypes. However, these methods require HLA reference panels which are often unavailable and memory-intensive. Cloud-based outsourced HLA imputation overcomes these limitations by utilizing built-in reference panels on the cloud server, removing local storage needs. However, uploading genotype data online raises privacy concerns. Although secure from third parties, the uploaded data can be misused by cloud server administrators.  Additionally, reference panels and their HLA imputation model on the cloud servers must be safeguarded against malicious clients. To address these privacy issues, we developed privateHLA, a secure HLA imputation method using homomorphic encryption. privateHLA securely performs HLA imputation on an outsourced server, protecting both client’s data and enhances model privacy. privateHLA outperformed SNP2HLA but had slightly lower accuracy than CookHLA, both plaintext-based methods.

## Installation

This repository contains Go code for privateHLA. One can run it as any other Go program.

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
Currently `dataset` flag accepts `HapMap_EUR, Korean`, and `prefix` flag accepts `A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1`. To automatically run all dataset and prefix, run `out.sh`.
```
sh cmd/out.sh
```
Benchmarks using Go benchmark is also available. Run
```
go test . -run=^$ -bench=.
```
