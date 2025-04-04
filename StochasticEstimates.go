package main

import (
	"math/rand"
)

func GenerateRandomString(length int, alphabet []string) string {
	if length == 0 {
		return ""
	}
	return GenerateRandomString(length-1, alphabet) + alphabet[rand.Intn(len(alphabet))]

}

func AverageDistanceEstimator(num_draws int, size int, alphabet []string, gap_cost int, BLOSOM_matrix map[string]map[string]int) float32 {
	sum := float32(0)
	for range num_draws {
		result, _, _, _ := GappedEditDistanceBLOSOM(GenerateRandomString(size, alphabet),
			GenerateRandomString(size, alphabet),
			make(map[edit_key]edit_info),
			make(map[edit_key]edit_info),
			make(map[edit_key]edit_info),
			make(map[edit_key]edit_info),
			gap_cost,
			BLOSOM_matrix)
		sum += float32(result[edit_key{size, size}].distance)
	}
	return sum / float32(num_draws)
}
