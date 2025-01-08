package hla

import (
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// HLAData stores all the data for HLA prediction for specific dataset and allele prefix.
type HLAData struct {
	// Alleles is a list of all the alleles.
	Alleles []string

	// ID is a list of all the IDs.
	ID []string

	// Weights is a weight for each alleles (W).
	//
	// AlleleName: []Weight
	Weights map[string][]float64

	// Snips is a list of all the snips for each ID (X).
	//
	// ID: []Snip
	Snips map[string][]float64
}

// LoadHLA loads the HLA data for the given dataset and allele prefix.
func LoadHLA(dataset, prefix string) (HLAData, error) {
	var err error

	data := HLAData{
		Alleles: make([]string, 0),
		ID:      make([]string, 0),
		Weights: make(map[string][]float64, 0),
		Snips:   make(map[string][]float64, 0),
	}

	filePath := fmt.Sprintf("data/%v/weights", dataset)
	fileNames, err := os.ReadDir(filePath)
	if err != nil {
		return HLAData{}, err
	}

	for _, file := range fileNames {
		name := strings.TrimSuffix(file.Name(), filepath.Ext(file.Name()))
		nameSplit := strings.Split(name, "*")
		if nameSplit[0] == prefix {
			data.Alleles = append(data.Alleles, nameSplit[1])
		}
	}

	if len(data.Alleles) == 0 {
		return HLAData{}, errors.New("no alleles found")
	}

	for _, allele := range data.Alleles {
		weightFileName := fmt.Sprintf("data/%v/weights/%v*%v.csv", dataset, prefix, allele)
		weightFile, err := os.Open(weightFileName)
		if err != nil {
			return HLAData{}, err
		}
		defer weightFile.Close()

		weightReader := csv.NewReader(weightFile)

		_, err = weightReader.Read() // Skip the header
		if err != nil {
			return HLAData{}, err
		}

		weight := make([]float64, 0)
		for {
			row, err := weightReader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				return HLAData{}, err
			}
			w, err := strconv.ParseFloat(row[0], 64)
			if err != nil {
				return HLAData{}, err
			}

			weight = append(weight, w)
		}
		data.Weights[allele] = weight
	}

	snipsFile, err := os.Open(fmt.Sprintf("data/%v/X/X_test.csv", dataset))
	if err != nil {
		return HLAData{}, err
	}
	defer snipsFile.Close()

	snipsReader := csv.NewReader(snipsFile)

	_, err = snipsReader.Read() // Skip the header
	if err != nil {
		return HLAData{}, err
	}

	for {
		row, err := snipsReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return HLAData{}, err
		}

		rowSplit := strings.Split(row[0], "\t")

		data.ID = append(data.ID, rowSplit[0])
		data.Snips[rowSplit[0]] = make([]float64, 1)
		data.Snips[rowSplit[0]][0] = 1 // Bias term
		for _, snip := range rowSplit[1:] {
			s, err := strconv.ParseFloat(snip, 64)
			if err != nil {
				return HLAData{}, err
			}
			data.Snips[rowSplit[0]] = append(data.Snips[rowSplit[0]], s)
		}

	}

	return data, nil
}
