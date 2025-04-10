#!/bin/bash

datasets=("HapMap_EUR" "Korean")
prefixes=("A" "B" "C" "DPA1" "DPB1" "DQA1" "DQB1" "DRB1")

for dataset in "${datasets[@]}"
do
    echo "" > ./out/"$dataset".out
    for prefix in "${prefixes[@]}"
    do
        go run cmd/main.go -prefix="$prefix" -dataset="$dataset" >> "./out/$dataset.out"
    done
done
