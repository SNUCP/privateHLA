package hla

import (
	"math"

	"github.com/sp301415/tfhe-go/math/vec"
	"github.com/sp301415/tfhe-go/tfhe"
)

// Client is the client that requires HLA imputation.
type Client struct {
	// Parameters is the parameters for client.
	Parameters tfhe.Parameters[uint64]
	// Encryptor is the encryptor that encrypts the data.
	Encryptor *tfhe.Encryptor[uint64]
}

// NewClient creates a new client.
func NewClient() *Client {
	params := FloatParamsLiteral.Compile()
	return &Client{
		Parameters: params,
		Encryptor:  tfhe.NewEncryptor(params),
	}
}

// EncryptSnips encrypts the snips of the given ID.
// Since the number of snips exceed the polynomial degree,
// we split the snips into multiple ciphertexts.
func (c *Client) EncryptSnips(id string) []tfhe.FourierGLWECiphertext[uint64] {
	snip := HLAData.Snips[id]

	chunkCount := int(math.Round(float64(len(snip)) / float64(c.Parameters.PolyDegree())))
	ctSnip := make([]tfhe.FourierGLWECiphertext[uint64], chunkCount)
	for i := range chunkCount {
		start := i * c.Parameters.PolyDegree()
		end := min((i+1)*c.Parameters.PolyDegree(), len(snip))

		pt := tfhe.NewGLWEPlaintext(c.Parameters)
		for j, jj := start, 0; j < end; j, jj = j+1, jj+1 {
			pt.Value.Coeffs[jj] = uint64(math.Round(snip[j] * ScaleSnip))
		}
		// Reverse packed values for later inner product.
		vec.ReverseInPlace(pt.Value.Coeffs)
		ctSnip[i] = c.Encryptor.EncryptFourierGLWEPlaintext(pt)
	}

	return ctSnip
}
