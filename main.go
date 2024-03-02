package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"

	"github.com/sp301415/tfhe-go/math/poly"
	"github.com/sp301415/tfhe-go/math/vec"
	"github.com/sp301415/tfhe-go/tfhe"
)

var params = tfhe.ParametersLiteral[uint64]{
	LWEDimension:    2033,
	GLWEDimension:   1,
	PolyDegree:      2048,
	PolyLargeDegree: 16384,

	LWEStdDev:  math.Exp(-51),
	GLWEStdDev: math.Exp(-52),

	BlockSize: 19,

	MessageModulus: 1 << 10,

	BootstrapParameters: tfhe.GadgetParametersLiteral[uint64]{
		Base:  1 << 17,
		Level: 2,
	},
	KeySwitchParameters: tfhe.GadgetParametersLiteral[uint64]{
		Base:  1 << 4,
		Level: 5,
	},

	BootstrapOrder: tfhe.OrderKeySwitchBlindRotate,
}

func randFloat() float64 {
	return 2*rand.Float64() - 1
}

func encodeFloat(x float64, delta float64) uint64 {
	return uint64(math.Round(x * delta))
}

func decodeFloat(x uint64, delta float64) float64 {
	return float64(int64(x)) / delta
}

func main() {
	// We use three parameters: Small for sign bits, Medium for indices.
	paramsSmallLiteral := params
	paramsSmallLiteral.PolyLargeDegree = params.PolyDegree

	paramsMediumLiteral := params
	paramsMediumLiteral.MessageModulus = 1 << 7
	paramsMediumLiteral.PolyLargeDegree = int(paramsMediumLiteral.MessageModulus) << 5 // Rounding Error: 5bits

	params := params.Compile()
	paramsSmall := paramsSmallLiteral.Compile()
	paramsMedium := paramsMediumLiteral.Compile()

	enc := tfhe.NewEncryptor(params)

	evk := enc.GenEvaluationKeyParallel()
	tfheEval := tfhe.NewEvaluator(params, evk)
	tfheEvalSmall := tfhe.NewEvaluator(paramsSmall, evk)
	tfheEvalMedium := tfhe.NewEvaluator(paramsMedium, evk)

	L := 32
	// Bound is an expected bound for X * delta.
	// All values should be in the range (-bound, bound).
	bound := 8.0
	deltaCt := math.Exp2(40)
	deltaPt := math.Exp2(62) / (deltaCt * bound)
	deltaFinal := deltaPt * deltaCt
	// 2 * bound * deltaFinal should have at least 1 bit of padding.

	eval := NewEvaluator(tfheEval, tfheEvalSmall, tfheEvalMedium, L, bound, deltaFinal)

	tests := 30
	correct := 0
	for t := 0; t < tests; t++ {
		fmt.Println("Test", t, "/", tests)

		data := make([]float64, L)
		for i := 0; i < L; i++ {
			data[i] = randFloat()
		}

		X := make([][]float64, L)
		for i := 0; i < L; i++ {
			X[i] = make([]float64, L)
			for j := 0; j < L; j++ {
				X[i][j] = randFloat()
			}
		}

		y, _, i0, _, i1 := plain(X, data)
		fmt.Printf("(%v, %v)\n", i0, i1)

		XEcd := make([]poly.FourierPoly, L)
		pBuff := enc.PolyEvaluator.NewPoly()
		for i := 0; i < L; i++ {
			for j := 0; j < L; j++ {
				pBuff.Coeffs[j] = encodeFloat(X[i][j], deltaPt)
			}
			XEcd[i] = enc.FourierEvaluator.ToFourierPoly(pBuff)
		}

		ptData := tfhe.NewGLWEPlaintext(params)
		for i := 0; i < L; i++ {
			ptData.Value.Coeffs[i] = encodeFloat(data[i], deltaCt)
		}
		vec.ReverseInPlace(ptData.Value.Coeffs)
		ctData := enc.EncryptFourierGLWEPlaintext(ptData)

		now := time.Now()
		ctResult := eval.Evaluate(enc, XEcd, ctData)
		fmt.Println("Encrypted", time.Since(now))

		// Output Bounds are always 1<<62.
		A0 := decodeFloat(enc.DecryptLWEPlaintext(ctResult.ctA0).Value, 1<<62)
		A1 := decodeFloat(enc.DecryptLWEPlaintext(ctResult.ctA1).Value, 1<<62)

		I0 := enc.DecodeLWECustom(enc.DecryptLWEPlaintext(ctResult.ctIdx0), paramsMedium.MessageModulus(), paramsMedium.Delta())
		I1 := enc.DecodeLWECustom(enc.DecryptLWEPlaintext(ctResult.ctIdx1), paramsMedium.MessageModulus(), paramsMedium.Delta())

		fmt.Printf("(%v, %v)\n", I0, I1)

		if i0 == I0 && i1 == I1 {
			fmt.Println("Correct!")
			correct++
		} else {
			fmt.Println("Incorrect!")
			fmt.Printf("%v, %v\n", i0, y[i0])
			fmt.Println("vs")
			fmt.Printf("%v, %v - %v\n", I0, y[I0], A0)
			fmt.Printf("%v, %v\n", i1, y[i1])
			fmt.Println("vs")
			fmt.Printf("%v, %v - %v\n", I1, y[I1], A1)
		}
	}

	fmt.Println("Correct:", correct, "/", tests)
}
