package main

import (
	"math"
)

func plain(X [][]float64, data []float64) ([]float64, float64, int, float64, int) {
	// y = X * data
	y := matVecMul(X, data)
	yy := make([]float64, len(y))

	// y = sigmoid(y)
	for i := range y {
		yy[i] = sigmoid(y[i])
	}

	// Top 2 elements
	t0, t1 := math.Inf(-1), math.Inf(-1)
	i0, i1 := -1, -1
	for i, v := range yy {
		if v > t0 {
			t1, i1 = t0, i0
			t0, i0 = v, i
		} else if v > t1 {
			t1, i1 = v, i
		}
	}

	return y, t0, i0, t1, i1
}

func matVecMul(A [][]float64, x []float64) []float64 {
	n := len(A)

	y := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			y[i] += A[i][j] * x[j]
		}
	}
	return y
}

func sigmoid(x float64) float64 {
	return 1 / (1 + math.Exp(-x))
}
