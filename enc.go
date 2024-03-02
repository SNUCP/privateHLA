package main

import (
	"math"
	"sync"

	"github.com/sp301415/tfhe-go/math/poly"
	"github.com/sp301415/tfhe-go/tfhe"
)

// Evaluator is a buffer for model evaluation.
type Evaluator struct {
	Evaluator       *tfhe.Evaluator[uint64]
	EvaluatorPool   []*tfhe.Evaluator[uint64]
	EvaluatorSmall  *tfhe.Evaluator[uint64]
	EvaluatorMedium *tfhe.Evaluator[uint64]

	Length int
	Bound  float64
	Delta  float64

	normLUT     tfhe.LookUpTable[uint64]
	signLUT     tfhe.LookUpTable[uint64]
	mulLUT      tfhe.LookUpTable[uint64]
	mulExactLUT tfhe.LookUpTable[uint64]

	ctProbs       []tfhe.LWECiphertext[uint64]
	ctProd        tfhe.GLWECiphertext[uint64]
	ctFourierProd tfhe.FourierGLWECiphertext[uint64]

	ctCmp0 tfhe.LWECiphertext[uint64]
	ctCmp1 tfhe.LWECiphertext[uint64]

	ctCmpIdx0 tfhe.LWECiphertext[uint64]
	ctCmpIdx1 tfhe.LWECiphertext[uint64]

	ctSub    tfhe.LWECiphertext[uint64]
	ctSubNeg tfhe.LWECiphertext[uint64]
	ctAdd    tfhe.LWECiphertext[uint64]
}

// EvaluateResult is a result of model evaluation.
type EvaluateResult struct {
	ctA0   tfhe.LWECiphertext[uint64]
	ctIdx0 tfhe.LWECiphertext[uint64]
	ctA1   tfhe.LWECiphertext[uint64]
	ctIdx1 tfhe.LWECiphertext[uint64]
}

// NewEvaluator creates a new Evaluator.
// weights are the model parameters.
// bound is a expected maximum value for X * data.
// All values should be in the range [-bound, bound] after linear transformation.
func NewEvaluator(e, eSmall, eMedium *tfhe.Evaluator[uint64], N int, bound, delta float64) *Evaluator {
	evaluatorPool := make([]*tfhe.Evaluator[uint64], 8)
	for i := range evaluatorPool {
		evaluatorPool[i] = e.ShallowCopy()
	}

	// TODO: We have to optimize this.
	// 1. Consider selecting optimal threshold value.
	// 2. Should we "erase" the error here? Currently we include the rounding error.
	eps := 1 / float64(e.Parameters.PolyLargeDegree()) // Arbitrary small value which is larger than the bootstrapping error
	threshold := 0.1
	normLUT := tfhe.NewLookUpTable(e.Parameters)
	for i := 0; i < e.Parameters.PolyLargeDegree(); i++ {
		in := float64(i) / float64(e.Parameters.PolyLargeDegree())
		in = 2*bound*in - bound
		out := eps
		if in > threshold {
			out = (in-threshold)*(1-eps)/(bound-threshold) + eps
		}

		normLUT.Coeffs[i] = uint64(math.Round(out * (1 << 62))) // Leave 2 bits of padding
	}

	signLUT := tfhe.NewLookUpTable(eSmall.Parameters)
	for i := 0; i < eSmall.Parameters.PolyLargeDegree()/2; i++ {
		signLUT.Coeffs[i] = 1 << 62
	}

	mulLUT := e.GenLookUpTable(func(x int) int {
		msb := (x >> (e.Parameters.MessageModulusLog() - 1)) & 1
		if msb == 0 {
			return 0
		}
		return x & (1<<(e.Parameters.MessageModulusLog()-1) - 1)
	})

	mulExactLUT := eMedium.GenLookUpTable(func(x int) int {
		msb := (x >> (eMedium.Parameters.MessageModulusLog() - 1)) & 1
		if msb == 0 {
			return 0
		}
		return x & (1<<(eMedium.Parameters.MessageModulusLog()-1) - 1)
	})

	ctProbs := make([]tfhe.LWECiphertext[uint64], N)
	for i := range ctProbs {
		ctProbs[i] = tfhe.NewLWECiphertext(e.Parameters)
	}

	return &Evaluator{
		Evaluator:       e,
		EvaluatorPool:   evaluatorPool,
		EvaluatorSmall:  eSmall,
		EvaluatorMedium: eMedium,

		Length: N,
		Bound:  bound,
		Delta:  delta,

		normLUT:     normLUT,
		signLUT:     signLUT,
		mulLUT:      mulLUT,
		mulExactLUT: mulExactLUT,

		ctProbs:       ctProbs,
		ctProd:        tfhe.NewGLWECiphertext(e.Parameters),
		ctFourierProd: tfhe.NewFourierGLWECiphertext(e.Parameters),

		ctCmp0: tfhe.NewLWECiphertext(e.Parameters),
		ctCmp1: tfhe.NewLWECiphertext(e.Parameters),

		ctCmpIdx0: tfhe.NewLWECiphertext(e.Parameters),
		ctCmpIdx1: tfhe.NewLWECiphertext(e.Parameters),

		ctSub:    tfhe.NewLWECiphertext(e.Parameters),
		ctSubNeg: tfhe.NewLWECiphertext(e.Parameters),
		ctAdd:    tfhe.NewLWECiphertext(e.Parameters),
	}
}

// Evaluate evalutes the model with the encoded data.
func (e *Evaluator) Evaluate(enc *tfhe.Encryptor[uint64], weights []poly.FourierPoly, data tfhe.FourierGLWECiphertext[uint64]) EvaluateResult {
	// y = X * data
	for i := 0; i < e.Length; i++ {
		e.InnerProduct(weights[i], data, e.ctProbs[i])
		// Move the encrypted values to the range [0, 2 * Bound).
		e.ctProbs[i].Value[0] += uint64(math.Round(e.Bound * e.Delta))
	}

	normJobs := make(chan int)
	go func() {
		defer close(normJobs)
		for i := 0; i < e.Length; i++ {
			normJobs <- i
		}
	}()

	var wg sync.WaitGroup
	wg.Add(len(e.EvaluatorPool))
	for i := range e.EvaluatorPool {
		go func(idx int) {
			defer wg.Done()
			for i := range normJobs {
				e.EvaluatorPool[idx].BootstrapLUTAssign(e.ctProbs[i], e.normLUT, e.ctProbs[i])
			}
		}(i)
	}
	wg.Wait()

	// Top 2 ciphertexts
	ctA0 := tfhe.NewLWECiphertext(e.Evaluator.Parameters)
	ctA1 := tfhe.NewLWECiphertext(e.Evaluator.Parameters)

	// Top 2 indices
	ctIdx0 := tfhe.NewLWECiphertext(e.Evaluator.Parameters)
	ctIdx1 := tfhe.NewLWECiphertext(e.Evaluator.Parameters)

	// Start with y[0] and y[1].
	e.SignBit(e.ctProbs[0], e.ctProbs[1], e.ctSub)
	e.Evaluator.NegLWEAssign(e.ctSub, e.ctSubNeg)
	e.ctSubNeg.Value[0] += 1 << 62 // 1 - ctSub
	e.MaxMin(e.ctSub, e.ctSubNeg, e.ctProbs[0], e.ctProbs[1], ctA0, ctA1)

	e.ctCmpIdx0.Clear()
	e.ctCmpIdx1.Clear()
	e.ctCmpIdx1.Value[0] = 1 << e.EvaluatorMedium.Parameters.DeltaLog()
	e.MaxMinExact(e.ctSub, e.ctSubNeg, e.ctCmpIdx0, e.ctCmpIdx1, ctIdx0, ctIdx1)

	// Loop for the rest of the ciphertexts,
	// Comparing y[i] and the current two maximums.
	for i := 2; i < e.Length; i++ {
		e.ctCmpIdx0.Clear()
		e.ctCmpIdx0.Value[0] = uint64(i) << e.EvaluatorMedium.Parameters.DeltaLog()

		e.SignBit(e.ctProbs[i], ctA0, e.ctSub)
		e.Evaluator.NegLWEAssign(e.ctSub, e.ctSubNeg)
		e.ctSubNeg.Value[0] += 1 << 62 // 1 - ctSub
		e.MaxMin(e.ctSub, e.ctSubNeg, e.ctProbs[i], ctA0, ctA0, e.ctProbs[i])
		e.MaxMinExact(e.ctSub, e.ctSubNeg, e.ctCmpIdx0, ctIdx0, ctIdx0, e.ctCmpIdx0)

		e.SignBit(e.ctProbs[i], ctA1, e.ctSub)
		e.Evaluator.NegLWEAssign(e.ctSub, e.ctSubNeg)
		e.ctSubNeg.Value[0] += 1 << 62 // 1 - ctSub
		e.MaxMin(e.ctSub, e.ctSubNeg, e.ctProbs[i], ctA1, ctA1, e.ctProbs[i])
		e.MaxMinExact(e.ctSub, e.ctSubNeg, e.ctCmpIdx0, ctIdx1, ctIdx1, e.ctCmpIdx0)
	}

	return EvaluateResult{
		ctA0:   ctA0,
		ctIdx0: ctIdx0,
		ctA1:   ctA1,
		ctIdx1: ctIdx1,
	}
}

// InnerProduct calculates the inner product of x and y.
func (e *Evaluator) InnerProduct(x poly.FourierPoly, y tfhe.FourierGLWECiphertext[uint64], ctOut tfhe.LWECiphertext[uint64]) {
	e.Evaluator.FourierPolyMulFourierGLWEAssign(y, x, e.ctFourierProd)
	e.Evaluator.ToStandardGLWECiphertextAssign(e.ctFourierProd, e.ctProd)
	e.Evaluator.SampleExtractAssign(e.ctProd, e.Evaluator.Parameters.PolyDegree()-1, ctOut)
}

// SignBit returns 2^(p-1) if ct0 > ct1, and 0 otherwise.
func (e *Evaluator) SignBit(ct0, ct1 tfhe.LWECiphertext[uint64], ctOut tfhe.LWECiphertext[uint64]) {
	e.EvaluatorSmall.SubLWEAssign(ct0, ct1, e.ctSub)
	e.EvaluatorSmall.BootstrapLUTAssign(e.ctSub, e.signLUT, ctOut)
}

// MaxMin reorders the ciphertexts in the order of the values, so that
// ctOut0 := max(ct0, ct1) and ctOut1 := min(ct0, ct1).
func (e *Evaluator) MaxMin(ctSign, ctSignNeg, ct0, ct1, ctOut0, ctOut1 tfhe.LWECiphertext[uint64]) {
	e.Evaluator.AddLWEAssign(ct0, ct1, e.ctAdd)

	e.Evaluator.AddLWEAssign(ct0, ctSign, e.ctCmp0)
	e.Evaluator.AddLWEAssign(ct1, ctSignNeg, e.ctCmp1)

	e.Evaluator.BootstrapLUTAssign(e.ctCmp0, e.mulLUT, ctOut0)
	e.Evaluator.BootstrapLUTAssign(e.ctCmp1, e.mulLUT, ctOut1)

	e.Evaluator.AddLWEAssign(ctOut0, ctOut1, ctOut0)
	e.Evaluator.SubLWEAssign(e.ctAdd, ctOut0, ctOut1)
}

// MaxMin reorders the ciphertexts in the order of the values, so that
// ctOut0 := max(ct0, ct1) and ctOut1 := min(ct0, ct1).
func (e *Evaluator) MaxMinExact(ctSign, ctSignNeg, ct0, ct1, ctOut0, ctOut1 tfhe.LWECiphertext[uint64]) {
	e.EvaluatorMedium.AddLWEAssign(ct0, ct1, e.ctAdd)

	e.EvaluatorMedium.AddLWEAssign(ct0, ctSign, e.ctCmp0)
	e.EvaluatorMedium.AddLWEAssign(ct1, ctSignNeg, e.ctCmp1)

	e.EvaluatorMedium.BootstrapLUTAssign(e.ctCmp0, e.mulExactLUT, ctOut0)
	e.EvaluatorMedium.BootstrapLUTAssign(e.ctCmp1, e.mulExactLUT, ctOut1)

	e.EvaluatorMedium.AddLWEAssign(ctOut0, ctOut1, ctOut0)
	e.EvaluatorMedium.SubLWEAssign(e.ctAdd, ctOut0, ctOut1)
}
