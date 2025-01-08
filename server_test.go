package hla_test

import (
	"fmt"
	"hla"
	"testing"
)

var (
	datasets     = [...]string{"HapMap_EUR", "Korean"}
	prefixs      = [...]string{"A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"}
	linearResult = hla.LinearResult{}
	top2Result   = hla.Top2Result{}
	result       = hla.ServerResult{}
)

func BenchmarkServer(b *testing.B) {
	for _, dataset := range datasets {
		for _, prefix := range prefixs {
			data, err := hla.LoadHLA(dataset, prefix)
			if err != nil {
				continue
			}

			client := hla.NewClient(data)
			evk := client.Encryptor.GenEvaluationKeyParallel()

			server := hla.NewServer(data, evk)
			// server.Encryptor = client.Encryptor

			benchName := fmt.Sprintf("dataset=%v/prefix=%v", dataset, prefix)
			b.Run(benchName+"/Full", func(b *testing.B) {
				id := data.ID[0]
				ctSnip := client.EncryptSnips(id)
				b.ResetTimer()

				for i := 0; i < b.N; i++ {
					result = server.Predict(ctSnip)
				}
			})

			b.Run(benchName+"/Linear", func(b *testing.B) {
				id := data.ID[0]
				ctSnip := client.EncryptSnips(id)
				b.ResetTimer()

				for i := 0; i < b.N; i++ {
					linearResult = server.Linear(ctSnip)
				}
			})

			b.Run(benchName+"/Top2", func(b *testing.B) {
				id := data.ID[0]
				ctSnip := client.EncryptSnips(id)
				linearResult := server.Linear(ctSnip)
				b.ResetTimer()

				for i := 0; i < b.N; i++ {
					top2Result = server.Top2(linearResult)
				}
			})

			b.Run(benchName+"/Threshold", func(b *testing.B) {
				id := data.ID[0]
				ctSnip := client.EncryptSnips(id)
				linearResult := server.Linear(ctSnip)
				top2Result := server.Top2(linearResult)
				b.ResetTimer()

				for i := 0; i < b.N; i++ {
					result = server.Threshold(top2Result)
				}
			})
		}
	}
}
