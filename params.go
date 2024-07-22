package hla

import (
	"math"

	"github.com/sp301415/tfhe-go/tfhe"
)

var (
	// FloatParamsLiteral is the parameters for floating-point operations.
	// It has much lower failure rate. (Around 2^-20).
	FloatParamsLiteral = tfhe.ParametersLiteral[uint64]{
		LWEDimension:    935,
		GLWEDimension:   1,
		PolyDegree:      2048,
		PolyLargeDegree: 4096,

		LWEStdDev:  0.000000277202526552345,
		GLWEStdDev: 0.00000000000000029403601535432533,

		BlockSize: 5,

		MessageModulus: 1 << 7,

		BootstrapParameters: tfhe.GadgetParametersLiteral[uint64]{
			Base:  1 << 23,
			Level: 1,
		},
		KeySwitchParameters: tfhe.GadgetParametersLiteral[uint64]{
			Base:  1 << 4,
			Level: 5,
		},

		BootstrapOrder: tfhe.OrderKeySwitchBlindRotate,
	}

	// IntParamsLiteral is the parameters for integer operations.
	// This has standard failure rate, around 2^-40.
	IntParamsLiteral = tfhe.ParametersLiteral[uint64]{
		LWEDimension:    935,
		GLWEDimension:   1,
		PolyDegree:      2048,
		PolyLargeDegree: 2048,

		LWEStdDev:  0.000000277202526552345,
		GLWEStdDev: 0.00000000000000029403601535432533,

		BlockSize: 5,

		MessageModulus: 1 << 5,

		BootstrapParameters: tfhe.GadgetParametersLiteral[uint64]{
			Base:  1 << 23,
			Level: 1,
		},
		KeySwitchParameters: tfhe.GadgetParametersLiteral[uint64]{
			Base:  1 << 4,
			Level: 5,
		},

		BootstrapOrder: tfhe.OrderKeySwitchBlindRotate,
	}

	// PredictionBound is a bound for the predicted values.
	// All predictions must be in the range (-PredictionBound, PredictionBound).
	PredictionBound = 40.0

	// ScaleSnip is the scale for the ciphertexts ("Snips").
	ScaleSnip = math.Exp2(40)

	// ScaleWeight is the scale for the plaintexts ("Weights").
	// Due to correctness, we require that 2*PredictionBound*ScaleSnip*ScaleWeight < Q/4,
	// i.e. has at least 1 bit of padding (on top of the usual negacyclic padding).
	ScaleWeight = math.Exp2(62) / (PredictionBound * ScaleSnip)

	// NormalizeBound is a bound for "Normalizing" the predictions.
	// We use a ReLU-like function.
	NormalizeBound = -10.0

	// CompareThreshold is a threshold for comparing the predictions.
	CompareThreshold = 7.0
)