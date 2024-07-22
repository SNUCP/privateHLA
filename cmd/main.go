package main

import (
	"cmp"
	"fmt"
	"hla"
	"math"
	"slices"
	"time"

	"github.com/sp301415/tfhe-go/math/vec"
)

type PredictionResult struct {
	Allele string
	Pred   float64
}

func sigmoid(x float64) float64 {
	return 1 / (1 + math.Exp(-x))
}

func main() {
	prefix := "B"

	cnt := 0
	predictions := make(map[string][]PredictionResult)
	predictionsFinal := make(map[string][2]PredictionResult)
	for _, id := range hla.HLAData.ID {
		predictions[id] = make([]PredictionResult, len(hla.HLAData.Alleles[prefix]))
		for j, allele := range hla.HLAData.Alleles[prefix] {
			predictions[id][j] = PredictionResult{
				Allele: hla.HLAData.Alleles[prefix][j],
				Pred:   sigmoid(vec.Dot(hla.HLAData.Weights[prefix][allele], hla.HLAData.Snips[id])),
			}
		}

		slices.SortFunc(predictions[id], func(i, j PredictionResult) int {
			return -cmp.Compare(i.Pred, j.Pred)
		})

		allele0, pred0 := predictions[id][0].Allele, predictions[id][0].Pred
		allele1, pred1 := predictions[id][1].Allele, predictions[id][1].Pred

		if pred0 > hla.CompareThreshold*pred1 {
			allele1 = allele0
			pred1 = pred0
		} else if pred1 > hla.CompareThreshold*pred0 {
			allele0 = allele1
			pred0 = pred1
		}

		if allele0 > allele1 {
			allele0, allele1 = allele1, allele0
			pred0, pred1 = pred1, pred0
		}

		predictionsFinal[id] = [2]PredictionResult{
			{Allele: allele0, Pred: pred0},
			{Allele: allele1, Pred: pred1},
		}

		if allele0 == hla.HLAData.Predictions[id][prefix][0] && allele1 == hla.HLAData.Predictions[id][prefix][1] {
			cnt++
		}
	}

	fmt.Println("Plain Accuracy:", float64(cnt)/float64(len(hla.HLAData.ID)))

	client := hla.NewClient()
	evk := client.Encryptor.GenEvaluationKeyParallel()

	server := hla.NewServer(prefix, evk)

	cnt = 0
	var elapsed time.Duration
	for i, id := range hla.HLAData.ID {
		fmt.Printf("Test %v of %v\n", i+1, len(hla.HLAData.ID))
		fmt.Println("ID:", id)

		ctSnip := client.EncryptSnips(id)

		now := time.Now()
		ctPred := server.Predict(ctSnip)
		infTime := time.Since(now)
		fmt.Println("Inference Time:", infTime)
		elapsed += infTime

		idx00 := server.IndexEvaluator.DecodeLWE(client.Encryptor.DecryptLWEPlaintext(ctPred.Idx00))
		idx01 := server.IndexEvaluator.DecodeLWE(client.Encryptor.DecryptLWEPlaintext(ctPred.Idx01))
		idx10 := server.IndexEvaluator.DecodeLWE(client.Encryptor.DecryptLWEPlaintext(ctPred.Idx10))
		idx11 := server.IndexEvaluator.DecodeLWE(client.Encryptor.DecryptLWEPlaintext(ctPred.Idx11))

		allele0 := hla.HLAData.Alleles[prefix][idx01*16+idx00]
		allele1 := hla.HLAData.Alleles[prefix][idx11*16+idx10]

		predRaw0 := float64(int64(client.Encryptor.DecryptLWEPlaintext(ctPred.Top0).Value)) / (1 << 62)
		predRaw1 := float64(int64(client.Encryptor.DecryptLWEPlaintext(ctPred.Top1).Value)) / (1 << 62)

		// Reverse ReLU
		predRaw0 = predRaw0*(hla.PredictionBound-hla.NormalizeBound) + hla.NormalizeBound
		predRaw1 = predRaw1*(hla.PredictionBound-hla.NormalizeBound) + hla.NormalizeBound

		pred0 := sigmoid(predRaw0)
		pred1 := sigmoid(predRaw1)

		if pred0 > hla.CompareThreshold*pred1 {
			allele1 = allele0
			pred1 = pred0
		} else if pred1 > hla.CompareThreshold*pred0 {
			allele0 = allele1
			pred0 = pred1
		}

		if allele0 > allele1 {
			allele0, allele1 = allele1, allele0
			pred0, pred1 = pred1, pred0
		}

		allele0Real := hla.HLAData.Predictions[id][prefix][0]
		allele1Real := hla.HLAData.Predictions[id][prefix][1]

		fmt.Println("Encrypted :", allele0, allele1)
		fmt.Println("Real      :", allele0Real, allele1Real)

		if allele0 == allele0Real && allele1 == allele1Real {
			fmt.Println("Correct :)")
			cnt++
		} else {
			fmt.Println("Incorrect :(")
			fmt.Println("Raw Data:")
			fmt.Printf("Encrypted                    %v - %v\n", allele0, pred0)
			fmt.Printf("Encrypted                    %v - %v\n", allele1, pred1)
			fmt.Printf("Encrypted (Plain Recomputed) %v - %v\n", allele0, sigmoid(vec.Dot(hla.HLAData.Snips[id], hla.HLAData.Weights[prefix][allele0])))
			fmt.Printf("Encrypted (Plain Recomputed) %v - %v\n", allele1, sigmoid(vec.Dot(hla.HLAData.Snips[id], hla.HLAData.Weights[prefix][allele1])))
			fmt.Printf("Plain                        %v - %v\n", predictionsFinal[id][0].Allele, predictionsFinal[id][0].Pred)
			fmt.Printf("Plain                        %v - %v\n", predictionsFinal[id][1].Allele, predictionsFinal[id][1].Pred)
			fmt.Printf("Real                         %v\n", allele0Real)
			fmt.Printf("Real                         %v\n", allele1Real)
		}
	}

	fmt.Println("Encrypted Accuracy:", float64(cnt)/float64(len(hla.HLAData.ID)))
	fmt.Println("Average Inference Time:", elapsed/time.Duration(len(hla.HLAData.ID)))
}
