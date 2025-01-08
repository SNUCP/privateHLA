package main

import (
	"flag"
	"fmt"
	"hla"
)

type PredictionResult struct {
	Allele string
	Pred   float64
}

func main() {
	var prefix string
	var dataset string
	flag.StringVar(&prefix, "prefix", "", "HLA prefix")
	flag.StringVar(&dataset, "dataset", "", "Dataset")
	flag.Parse()

	switch dataset {
	case "HapMap_EUR":
		hla.CompareThreshold = 17
	case "Korean":
		hla.CompareThreshold = 15
	}

	data, err := hla.LoadHLA(dataset, prefix)
	if err != nil {
		return
	}

	client := hla.NewClient(data)
	evk := client.Encryptor.GenEvaluationKeyParallel()

	server := hla.NewServer(data, evk)
	server.Encryptor = client.Encryptor

	// cnt = 0.0t
	for _, id := range data.ID {

		ctSnip := client.EncryptSnips(id)

		ctPred := server.Predict(ctSnip)

		allele0, allele1 := client.DecryptResults(prefix, ctPred)

		fmt.Printf("%v %v %v ", id, id, prefix)

		if allele0 == "0" {
			allele0 = "0000"
		}
		if allele1 == "0" {
			allele1 = "0000"
		}

		fmt.Printf("%v,%v %v,%v ", allele0[:2], allele1[:2], allele0, allele1)

		if allele0 == allele1 {
			fmt.Printf("1 1 1\n")
		} else {
			fmt.Printf("0.5 0.5 1\n")
		}
	}
}
