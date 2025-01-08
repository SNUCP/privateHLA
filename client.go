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
	// IdxParameters is the parameters for indices.
	IdxParameters tfhe.Parameters[uint64]
	// Encryptor is the encryptor that encrypts the data.
	Encryptor *tfhe.Encryptor[uint64]

	// Data is the HLA data.
	Data HLAData
}

// NewClient creates a new client.
func NewClient(data HLAData) *Client {
	params := FloatParamsLiteral.Compile()
	return &Client{
		Parameters:    params,
		IdxParameters: IntParamsLiteral.Compile(),
		Encryptor:     tfhe.NewEncryptor(params),
		Data:          data,
	}
}

// EncryptSnips encrypts the snips of the given ID.
// Since the number of snips exceed the polynomial degree,
// we split the snips into multiple ciphertexts.
func (c *Client) EncryptSnips(id string) []tfhe.FourierGLWECiphertext[uint64] {
	snip := c.Data.Snips[id]

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

func sigmoid(x float64) float64 {
	return 1 / (1 + math.Exp(-x))
}

// DecryptResults decrypts the prediction results from the server.
func (c *Client) DecryptResults(prefix string, res ServerResult) (string, string) {
	idx00 := c.Encryptor.DecodeLWECustom(c.Encryptor.DecryptLWEPlaintext(res.Idx00), c.IdxParameters.MessageModulus(), c.IdxParameters.Scale())
	idx01 := c.Encryptor.DecodeLWECustom(c.Encryptor.DecryptLWEPlaintext(res.Idx01), c.IdxParameters.MessageModulus(), c.IdxParameters.Scale())
	idx10 := c.Encryptor.DecodeLWECustom(c.Encryptor.DecryptLWEPlaintext(res.Idx10), c.IdxParameters.MessageModulus(), c.IdxParameters.Scale())
	idx11 := c.Encryptor.DecodeLWECustom(c.Encryptor.DecryptLWEPlaintext(res.Idx11), c.IdxParameters.MessageModulus(), c.IdxParameters.Scale())

	idx0 := idx01*int(c.IdxParameters.MessageModulus()/2) + idx00
	idx1 := idx11*int(c.IdxParameters.MessageModulus()/2) + idx10

	if idx0 > len(c.Data.Alleles) || idx1 > len(c.Data.Alleles) {
		// This should not happen, but if it did the result is wrong.
		// Just output the first allele.
		return c.Data.Alleles[0], c.Data.Alleles[0]
	}

	allele0 := c.Data.Alleles[idx0]
	allele1 := c.Data.Alleles[idx1]

	return allele0, allele1
}
