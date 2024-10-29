package hla

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// HLAData is a singleton variable for all the data for HLA prediction.
var HLAData = struct {
	// Alleles is a list of all the alleles.
	//
	// AllelePrefix: []AlleleName
	Alleles map[string][]string
	// ID is a list of all the IDs.
	ID []string
	// Weights is a weight for each alleles (W).
	//
	// AllelePrefix: AlleleName: []Weight
	Weights map[string]map[string][]float64
	// Snips is a list of all the snips for each ID (X).
	//
	// ID: []Snip
	Snips map[string][]float64
	// Predictions is a list of all the predictions for each ID (Raw data of Y).
	//
	// ID: AllelePrefix: Predictions
	Predictions map[string]map[string][2]string
}{}

func LoadHLA(dataset int) {
	var err error

	// Fill the HLAData with the data from the CSV files.
	allelesFile, err := os.Open(fmt.Sprintf("data/%v/weights/alleles.csv", dataset))
	if err != nil {
		panic(err)
	}
	defer allelesFile.Close()

	HLAData.Alleles = make(map[string][]string, 0)
	allelesReader := csv.NewReader(bufio.NewReader(allelesFile))

	for {
		row, err := allelesReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			panic(err)
		}
		alleleCut := strings.Split(row[0], "*")

		HLAData.Alleles[alleleCut[0]] = append(HLAData.Alleles[alleleCut[0]], alleleCut[1])
	}

	HLAData.Weights = make(map[string]map[string][]float64, 0)
	for allelePrefix, alleleNames := range HLAData.Alleles {
		HLAData.Weights[allelePrefix] = make(map[string][]float64, 0)
		for _, allele := range alleleNames {
			weightFileName := fmt.Sprintf("data/%v/weights/w_%s_%s.csv", dataset, allelePrefix, allele)
			weightFile, err := os.Open(weightFileName)
			if err != nil {
				panic(err)
			}

			weight := make([]float64, 0)
			weightReader := csv.NewReader(bufio.NewReader(weightFile))

			_, err = weightReader.Read() // Skip the header
			if err != nil {
				panic(err)
			}

			for {
				row, err := weightReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					panic(err)
				}
				w, err := strconv.ParseFloat(row[0], 64)
				if err != nil {
					panic(err)
				}

				weight = append(weight, w)
			}

			HLAData.Weights[allelePrefix][allele] = weight
		}
	}

	snipsFile, err := os.Open(fmt.Sprintf("data/%v/X/X_test.csv", dataset))
	if err != nil {
		panic(err)
	}

	HLAData.ID = make([]string, 0)
	HLAData.Snips = make(map[string][]float64)
	snipsReader := csv.NewReader(bufio.NewReader(snipsFile))

	_, err = snipsReader.Read() // Skip the header
	if err != nil {
		panic(err)
	}

	for {
		row, err := snipsReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			panic(err)
		}

		rowSplit := strings.Split(row[0], "\t")

		HLAData.ID = append(HLAData.ID, rowSplit[0])
		HLAData.Snips[rowSplit[0]] = make([]float64, 1)
		HLAData.Snips[rowSplit[0]][0] = 1 // Bias term
		for _, snip := range rowSplit[1:] {
			s, _ := strconv.ParseFloat(snip, 64)
			HLAData.Snips[rowSplit[0]] = append(HLAData.Snips[rowSplit[0]], s)
		}
	}

	// HLAData.Predictions = make(map[string]map[string][2]string, 0)
	// for _, id := range HLAData.ID {
	// 	HLAData.Predictions[id] = make(map[string][2]string, 0)
	// }

	// for allelePrefix := range HLAData.Alleles {
	// 	predictionFileName := fmt.Sprintf("data/Y_val_%s.csv", allelePrefix)
	// 	predictionFile, err := os.Open(predictionFileName)
	// 	if err != nil {
	// 		panic(err)
	// 	}

	// 	predictionReader := csv.NewReader(bufio.NewReader(predictionFile))

	// 	alleleRow, err := predictionReader.Read() // Skip the header
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	alleles := strings.Split(alleleRow[0], "\t")[1:]

	// 	for {
	// 		row, err := predictionReader.Read()
	// 		if err != nil {
	// 			if err == io.EOF {
	// 				break
	// 			}
	// 			panic(err)
	// 		}

	// 		rowSplit := strings.Split(row[0], "\t")

	// 		id := rowSplit[0]
	// 		ptr := 0
	// 		prediction := [2]string{}
	// 		for i, p := range rowSplit[1:] {
	// 			if p == "1" {
	// 				prediction[ptr] = alleles[i]
	// 				ptr++
	// 			} else if p == "2" {
	// 				prediction[0] = alleles[i]
	// 				prediction[1] = alleles[i]
	// 				break
	// 			}
	// 		}

	// 		HLAData.Predictions[id][allelePrefix] = prediction
	// 	}
	// }
}
