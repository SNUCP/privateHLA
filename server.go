package hla

import (
	"math"
	"sync"

	"github.com/sp301415/tfhe-go/math/num"
	"github.com/sp301415/tfhe-go/math/poly"
	"github.com/sp301415/tfhe-go/tfhe"
)

// Server is a struct for predicting each allele prefix.
type Server struct {
	// Only for debugging purposes.
	Encryptor *tfhe.Encryptor[uint64]

	// Parameters is a "default" parameter for floating-point operations.
	Parameters tfhe.Parameters[uint64]
	// IndexParameters is a parameter for indices of alleles.
	IndexParameters tfhe.Parameters[uint64]
	// SignParameters is a parameter for computing the comparison.
	SignParameters tfhe.Parameters[uint64]

	// Evaluator is an evaluator for this component.
	Evaluator *tfhe.Evaluator[uint64]
	// EvaluatorPool is a pool of evaluators for this component.
	EvaluatorPool []*tfhe.Evaluator[uint64]
	// IndexEvaluator is an evaluator for indices of alleles.
	IndexEvaluator *tfhe.Evaluator[uint64]
	// SignEvaluator is an evaluator for computing the comparison.
	SignEvaluator *tfhe.Evaluator[uint64]

	// NormalizeLUT is a LUT for normalizing the predictions.
	NormalizeLUT tfhe.LookUpTable[uint64]
	// UnNormalizeLUT is a LUT for unnormalizing the predictions.
	UnNormalizeLUT tfhe.LookUpTable[uint64]
	// UnNormalizeThresholdLUT is UnNormalizeLUT / CompareThreshold.
	UnNormalizeThresholdLUT tfhe.LookUpTable[uint64]
	// SignLUT is a LUT for computing the sign.
	SignLUT tfhe.LookUpTable[uint64]
	// MulLUT is a LUT for comparing the predictions.
	MulLUT tfhe.LookUpTable[uint64]
	// MulIndexLUT is a LUT for comparing the indices.
	MulIndexLUT tfhe.LookUpTable[uint64]

	// Weights is a slice of encoded weights.
	// Weights[i] is the weight for the i-th allele (ordered as the same order in HLAData).
	// Since the number of weights exceed the polynomial degree,
	// we split the weights into multiple plaintexts.
	Weights [][]poly.FourierPoly

	ctSignSub tfhe.LWECiphertext[uint64]
	ctSign    tfhe.LWECiphertext[uint64]
	ctSignNeg tfhe.LWECiphertext[uint64]

	ctTop0      tfhe.LWECiphertext[uint64]
	ctTop1      tfhe.LWECiphertext[uint64]
	ctTop0Cmp   tfhe.LWECiphertext[uint64]
	ctTop1Cmp   tfhe.LWECiphertext[uint64]
	ctMaxMinA   tfhe.LWECiphertext[uint64]
	ctMaxMinB   tfhe.LWECiphertext[uint64]
	ctMaxMinAdd tfhe.LWECiphertext[uint64]

	ctTopIdx00     tfhe.LWECiphertext[uint64]
	ctTopIdx01     tfhe.LWECiphertext[uint64]
	ctTopIdx10     tfhe.LWECiphertext[uint64]
	ctTopIdx11     tfhe.LWECiphertext[uint64]
	ctMaxMinIdxA   tfhe.LWECiphertext[uint64]
	ctMaxMinIdxB   tfhe.LWECiphertext[uint64]
	ctMaxMinIdxAdd tfhe.LWECiphertext[uint64]

	ctIdx0       tfhe.LWECiphertext[uint64]
	ctIdx1       tfhe.LWECiphertext[uint64]
	ctBuffResult tfhe.LWECiphertext[uint64]
}

// NewServer creates a new ServerComponent.
func NewServer(data HLAData, evk tfhe.EvaluationKey[uint64]) *Server {
	params := FloatParamsLiteral.Compile()

	indexParams := IntParamsLiteral.Compile()
	signParams := IntParamsLiteral.WithMessageModulus(1 << 1).Compile()

	evaluator := tfhe.NewEvaluator(params, evk)
	indexEvaluator := tfhe.NewEvaluator(indexParams, evk)
	signEvaluator := tfhe.NewEvaluator(signParams, evk)
	evPool := make([]*tfhe.Evaluator[uint64], num.Sqrt(len(data.Alleles)))
	for i := range evPool {
		evPool[i] = evaluator.ShallowCopy()
	}

	// Normalize:
	// We use a ReLU-like function.
	normalizeLUT := evaluator.GenLookUpTable(func(x int) int {
		// Currently input is mapped as [-PredictionBound, PredictionBound] -> [0, 2*PredictionBound] -> [0, MessageModulus/2].
		f := float64(x) / float64(params.MessageModulus())
		f = 2*PredictionBound*f - PredictionBound
		if f < NormalizeBound {
			return 0
		}
		f = (f - NormalizeBound) / (PredictionBound - NormalizeBound)
		return int(math.Round(f * float64(params.MessageModulus()/2)))
	})

	// Unnormalize:
	// Inverse of normalizeLUT + sigmoid.
	unNormalizeLUT := evaluator.GenLookUpTable(func(x int) int {
		f := float64(2*x) / float64(params.MessageModulus())
		f = NormalizeBound + (PredictionBound-NormalizeBound)*f
		f = sigmoid(f)
		return int(math.Round(f * float64(params.MessageModulus()/2)))
	})

	unNormalizeThresholdLUT := evaluator.GenLookUpTable(func(x int) int {
		f := float64(2*x) / float64(params.MessageModulus())
		f = NormalizeBound + (PredictionBound-NormalizeBound)*f
		f = sigmoid(f) / CompareThreshold
		return int(math.Round(f * float64(params.MessageModulus()/2)))
	})

	// Sign:
	// We map sign(x).
	signLUT := tfhe.NewLookUpTable(signParams)
	for i := 0; i < signParams.LookUpTableSize()/2; i++ {
		signLUT.Value[i] = 1 << 62
	}

	// Mul:
	// If the first bit of x is 1, map x -> 1 * LSB(x) = LSB(x).
	// Otherwise, we map 0 * x = 0.
	mulLUT := evaluator.GenLookUpTable(func(x int) int {
		bound := int(params.MessageModulus() / 2)
		if x >= bound {
			return x - bound
		}
		return 0
	})
	mulIndexLUT := indexEvaluator.GenLookUpTable(func(x int) int {
		bound := int(indexParams.MessageModulus() / 2)
		if x >= bound {
			return x - bound
		}
		return 0
	})

	// Encode Weights
	weights := make([][]poly.FourierPoly, len(data.Alleles))
	for i, allele := range data.Alleles {
		weight := data.Weights[allele]
		chunkCount := int(math.Round(float64(len(weight)) / float64(params.PolyDegree())))
		weights[i] = make([]poly.FourierPoly, chunkCount)
		for j := range chunkCount {
			start := j * params.PolyDegree()
			end := min((j+1)*params.PolyDegree(), len(weight))

			pt := poly.NewPoly[uint64](params.PolyDegree())
			for k, kk := start, 0; k < end; k, kk = k+1, kk+1 {
				pt.Coeffs[kk] = uint64(math.Round(weight[k] * ScaleWeight))
			}

			weights[i][j] = evaluator.PolyEvaluator.ToFourierPoly(pt)
		}
	}

	return &Server{
		Parameters:      params,
		IndexParameters: indexParams,
		SignParameters:  signParams,

		Evaluator:      evaluator,
		EvaluatorPool:  evPool,
		IndexEvaluator: indexEvaluator,
		SignEvaluator:  signEvaluator,

		NormalizeLUT:            normalizeLUT,
		UnNormalizeLUT:          unNormalizeLUT,
		UnNormalizeThresholdLUT: unNormalizeThresholdLUT,
		SignLUT:                 signLUT,
		MulLUT:                  mulLUT,
		MulIndexLUT:             mulIndexLUT,

		Weights: weights,

		ctSignSub: tfhe.NewLWECiphertext(params),
		ctSign:    tfhe.NewLWECiphertext(signParams),
		ctSignNeg: tfhe.NewLWECiphertext(signParams),

		ctTop0:      tfhe.NewLWECiphertext(params),
		ctTop1:      tfhe.NewLWECiphertext(params),
		ctTop0Cmp:   tfhe.NewLWECiphertext(params),
		ctTop1Cmp:   tfhe.NewLWECiphertext(params),
		ctMaxMinA:   tfhe.NewLWECiphertext(params),
		ctMaxMinB:   tfhe.NewLWECiphertext(params),
		ctMaxMinAdd: tfhe.NewLWECiphertext(params),

		ctTopIdx00:     tfhe.NewLWECiphertext(indexParams),
		ctTopIdx01:     tfhe.NewLWECiphertext(indexParams),
		ctTopIdx10:     tfhe.NewLWECiphertext(indexParams),
		ctTopIdx11:     tfhe.NewLWECiphertext(indexParams),
		ctMaxMinIdxA:   tfhe.NewLWECiphertext(indexParams),
		ctMaxMinIdxB:   tfhe.NewLWECiphertext(indexParams),
		ctMaxMinIdxAdd: tfhe.NewLWECiphertext(indexParams),

		ctIdx0:       tfhe.NewLWECiphertext(indexParams),
		ctIdx1:       tfhe.NewLWECiphertext(indexParams),
		ctBuffResult: tfhe.NewLWECiphertext(params),
	}
}

// ServerResult holds the result of the prediction.
type ServerResult struct {
	Idx00 tfhe.LWECiphertext[uint64]
	Idx01 tfhe.LWECiphertext[uint64]

	Idx10 tfhe.LWECiphertext[uint64]
	Idx11 tfhe.LWECiphertext[uint64]
}

type LinearResult struct {
	Preds []tfhe.LWECiphertext[uint64]
}

type Top2Result struct {
	Top0 tfhe.LWECiphertext[uint64]
	Top1 tfhe.LWECiphertext[uint64]

	TopIdx00 tfhe.LWECiphertext[uint64]
	TopIdx01 tfhe.LWECiphertext[uint64]

	TopIdx10 tfhe.LWECiphertext[uint64]
	TopIdx11 tfhe.LWECiphertext[uint64]

	IdxExceedsBound bool
}

// Linear executes the linear part of the prediction.
func (s *Server) Linear(snips []tfhe.FourierGLWECiphertext[uint64]) LinearResult {
	preds := make([]tfhe.LWECiphertext[uint64], len(s.Weights))
	for i := range s.Weights {
		if len(s.Weights[i]) != len(snips) {
			panic("Invalid length of snips or weights")
		}

		ctProd := s.Evaluator.FourierPolyMulFourierGLWE(snips[0], s.Weights[i][0])
		for j := 1; j < len(snips); j++ {
			s.Evaluator.FourierPolyMulAddFourierGLWEAssign(snips[j], s.Weights[i][j], ctProd)
		}

		preds[i] = (s.Evaluator.ToGLWECiphertext(ctProd).ToLWECiphertext(s.Parameters.PolyDegree() - 1))
	}

	return LinearResult{
		Preds: preds,
	}
}

// Top2 executes the top2 part of the prediction.
func (s *Server) Top2(linearResult LinearResult) Top2Result {
	preds := make([]tfhe.LWECiphertext[uint64], len(linearResult.Preds))
	for i := range linearResult.Preds {
		preds[i] = linearResult.Preds[i].Copy()
	}

	normalizeIdxChan := make(chan int)
	go func() {
		for i := range preds {
			normalizeIdxChan <- i
		}
		close(normalizeIdxChan)
	}()

	var wg sync.WaitGroup
	for i := range s.EvaluatorPool {
		wg.Add(1)
		go func(idx int) {
			defer wg.Done()
			ev := s.EvaluatorPool[idx]
			for i := range normalizeIdxChan {
				// Add Q/4 to the prediction to make it positive.
				preds[i].Value[0] += 1 << 62
				// Normalize the prediction.
				ev.BootstrapLUTAssign(preds[i], s.NormalizeLUT, preds[i])
			}
		}(i)
	}
	wg.Wait()

	// Trivial encryption of index.
	s.ctTopIdx00.Clear()
	s.ctTopIdx00.Value[0] = s.IndexEvaluator.EncodeLWE(0).Value
	s.ctTopIdx01.Clear()

	s.ctTopIdx10.Clear()
	s.ctTopIdx10.Value[0] = s.IndexEvaluator.EncodeLWE(1).Value
	s.ctTopIdx11.Clear()

	idxExceedsBound := len(preds) > int(s.IndexParameters.MessageModulus()/2)

	s.SignBitAssign(preds[0], preds[1], s.ctSign, s.ctSignNeg)
	s.MaxMinAssign(s.ctSign, s.ctSignNeg, preds[0], preds[1], s.ctTop0, s.ctTop1)
	s.MaxMinIndexAssign(s.ctSign, s.ctSignNeg, s.ctTopIdx00, s.ctTopIdx10, s.ctTopIdx00, s.ctTopIdx10)
	if idxExceedsBound {
		s.MaxMinIndexAssign(s.ctSign, s.ctSignNeg, s.ctTopIdx01, s.ctTopIdx11, s.ctTopIdx01, s.ctTopIdx11)
	}

	for i := 2; i < len(preds); i++ {
		s.ctIdx0.Clear()
		s.ctIdx0.Value[0] = s.IndexEvaluator.EncodeLWE(i % (int(s.IndexParameters.MessageModulus()) / 2)).Value
		s.ctIdx1.Clear()
		s.ctIdx1.Value[0] = s.IndexEvaluator.EncodeLWE(i / (int(s.IndexParameters.MessageModulus()) / 2)).Value

		s.SignBitAssign(preds[i], s.ctTop0, s.ctSign, s.ctSignNeg)

		s.MaxMinAssign(s.ctSign, s.ctSignNeg, preds[i], s.ctTop0, s.ctTop0, s.ctBuffResult)
		s.MaxMinIndexAssign(s.ctSign, s.ctSignNeg, s.ctIdx0, s.ctTopIdx00, s.ctTopIdx00, s.ctIdx0)
		if idxExceedsBound {
			s.MaxMinIndexAssign(s.ctSign, s.ctSignNeg, s.ctIdx1, s.ctTopIdx01, s.ctTopIdx01, s.ctIdx1)
		}

		s.SignBitAssign(s.ctBuffResult, s.ctTop1, s.ctSign, s.ctSignNeg)

		s.MaxMinAssign(s.ctSign, s.ctSignNeg, s.ctBuffResult, s.ctTop1, s.ctTop1, s.ctBuffResult)
		s.MaxMinIndexAssign(s.ctSign, s.ctSignNeg, s.ctIdx0, s.ctTopIdx10, s.ctTopIdx10, s.ctIdx0)
		if idxExceedsBound {
			s.MaxMinIndexAssign(s.ctSign, s.ctSignNeg, s.ctIdx1, s.ctTopIdx11, s.ctTopIdx11, s.ctIdx1)
		}
	}

	return Top2Result{
		Top0: s.ctTop0.Copy(),
		Top1: s.ctTop1.Copy(),

		TopIdx00: s.ctTopIdx00.Copy(),
		TopIdx01: s.ctTopIdx01.Copy(),

		TopIdx10: s.ctTopIdx10.Copy(),
		TopIdx11: s.ctTopIdx11.Copy(),

		IdxExceedsBound: idxExceedsBound,
	}
}

// Threshold executes the threshold part of the prediction.
func (s *Server) Threshold(top2Result Top2Result) ServerResult {
	ctTop0 := top2Result.Top0
	ctTop1 := top2Result.Top1

	ctTopIdx00 := top2Result.TopIdx00
	ctTopIdx01 := top2Result.TopIdx01

	ctTopIdx10 := top2Result.TopIdx10
	ctTopIdx11 := top2Result.TopIdx11

	idxExceedsBound := top2Result.IdxExceedsBound

	// Since CompareThreshold * top0 > top0 > top1,
	// we only compare CompareThreshold * top1 and top0.
	// or, equivalently, top1 and top0 / CompareThreshold.
	s.Evaluator.BootstrapLUTAssign(ctTop0, s.UnNormalizeThresholdLUT, s.ctTop0Cmp)
	s.Evaluator.BootstrapLUTAssign(ctTop1, s.UnNormalizeLUT, s.ctTop1Cmp)

	s.SignBitAssign(s.ctTop1Cmp, s.ctTop0Cmp, s.ctSign, s.ctSignNeg)
	s.MaxMinIndexAssign(s.ctSign, s.ctSignNeg, ctTopIdx10, ctTopIdx00, ctTopIdx10, s.ctBuffResult)

	if idxExceedsBound {
		s.MaxMinIndexAssign(s.ctSign, s.ctSignNeg, ctTopIdx11, ctTopIdx01, ctTopIdx11, s.ctBuffResult)
	}

	return ServerResult{
		Idx00: ctTopIdx00.Copy(),
		Idx01: ctTopIdx01.Copy(),

		Idx10: ctTopIdx10.Copy(),
		Idx11: ctTopIdx11.Copy(),
	}
}

// Predict predicts the value of the given snips.
func (s *Server) Predict(snips []tfhe.FourierGLWECiphertext[uint64]) ServerResult {
	linearResult := s.Linear(snips)
	top2Result := s.Top2(linearResult)
	return s.Threshold(top2Result)
}

// SignBitAssign returns 2^(p-1) if ct0 > ct1, and 0 otherwise.
func (s *Server) SignBitAssign(ct0, ct1, ctSign, ctSignNeg tfhe.LWECiphertext[uint64]) {
	s.Evaluator.SubLWEAssign(ct0, ct1, s.ctSignSub)
	s.SignEvaluator.BootstrapLUTAssign(s.ctSignSub, s.SignLUT, s.ctSign)
	s.SignEvaluator.NegLWEAssign(s.ctSign, ctSignNeg)
	ctSignNeg.Value[0] += 1 << 62
}

// MaxMinAssign reorders the ciphertexts in the order of the values, so that
// ctOut0 := max(ct0, ct1) and ctOut1 := min(ct0, ct1).
func (s *Server) MaxMinAssign(ctSign, ctSignNeg, ct0, ct1, ctOut0, ctOut1 tfhe.LWECiphertext[uint64]) {
	s.Evaluator.AddLWEAssign(ct0, ctSign, s.ctMaxMinA)
	s.Evaluator.AddLWEAssign(ct1, ctSignNeg, s.ctMaxMinB)

	s.Evaluator.BootstrapLUTAssign(s.ctMaxMinA, s.MulLUT, s.ctMaxMinA)
	s.Evaluator.BootstrapLUTAssign(s.ctMaxMinB, s.MulLUT, s.ctMaxMinB)

	s.Evaluator.AddLWEAssign(ct0, ct1, s.ctMaxMinAdd)
	s.Evaluator.AddLWEAssign(s.ctMaxMinA, s.ctMaxMinB, ctOut0)
	s.Evaluator.SubLWEAssign(s.ctMaxMinAdd, ctOut0, ctOut1)
}

// MaxMinIndexAssign is the same as MaxMinAssign, but for indices.
func (s *Server) MaxMinIndexAssign(ctSign, ctSignNeg, ct0, ct1, ctOut0, ctOut1 tfhe.LWECiphertext[uint64]) {
	s.IndexEvaluator.AddLWEAssign(ct0, ctSign, s.ctMaxMinIdxA)
	s.IndexEvaluator.AddLWEAssign(ct1, ctSignNeg, s.ctMaxMinIdxB)

	s.IndexEvaluator.BootstrapLUTAssign(s.ctMaxMinIdxA, s.MulIndexLUT, s.ctMaxMinIdxA)
	s.IndexEvaluator.BootstrapLUTAssign(s.ctMaxMinIdxB, s.MulIndexLUT, s.ctMaxMinIdxB)

	s.IndexEvaluator.AddLWEAssign(ct0, ct1, s.ctMaxMinIdxAdd)
	s.IndexEvaluator.AddLWEAssign(s.ctMaxMinIdxA, s.ctMaxMinIdxB, ctOut0)
	s.IndexEvaluator.SubLWEAssign(s.ctMaxMinIdxAdd, ctOut0, ctOut1)
}
