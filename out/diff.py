import re, collections

prefix = ["A" ,"B", "C", "DPA1", "DPB1" ,"DQA1" ,"DQB1", "DRB1"]

plain = dict()
encrypted = dict()

for p in prefix:
    plain[p] = []
    encrypted[p] = []

with open("./out/1.real.out") as f:
    for line in f:
        l = line.strip().split(" ")
        id = l[0]
        p = l[2]
        allele0, allele1 = l[4].split(",")
        plain[p].append((id, allele0))
        plain[p].append((id, allele1))

with open("./out/1.out") as f:
    for line in f:
        l = line.strip().split(" ")
        id = l[0]
        p = l[2]
        allele0, allele1 = l[4].split(",")
        encrypted[p].append((id, allele0))
        encrypted[p].append((id, allele1))

for p in prefix:
    plain_counter = collections.Counter(plain[p])
    encrypted_counter = collections.Counter(encrypted[p])

    intersection = plain_counter & encrypted_counter

    print(p, len(intersection), len(intersection) / len(plain[p]))
