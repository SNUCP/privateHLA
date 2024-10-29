#!/bin/bash

arr=("A" "B" "C" "DPA1" "DPB1" "DQA1" "DQB1" "DRB1")

for prefix in "${arr[@]}"
do
    taskset -c 40-48 go run cmd/main.go -prefix="$prefix" -dataset=1 >> ./out/1.out
done
