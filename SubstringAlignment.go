package main

import (
	"maps"
)

func SuffixAlignment(str1 string, str2 string, del_cost int, BLOSOM_matrix map[string]map[string]int, edit_map map[edit_key]int) map[edit_key]int {

	m := len(str1)
	n := len(str2)
	key := edit_key{m, n}

	if str1 == "" {
		edit_map[key] = 0
		return edit_map
	}
	if str2 == "" {
		edit_map[key] = 0
		return edit_map
	}

	prev_ins_key := edit_key{m, n - 1}
	_, t := edit_map[prev_ins_key]
	if !t {
		ins_map := SuffixAlignment(str1, str2[0:n-1], del_cost, BLOSOM_matrix, edit_map)
		maps.Copy(edit_map, ins_map)
	}

	prev_del_key := edit_key{m - 1, n}
	_, t = edit_map[prev_del_key]
	if !t {
		del_map := SuffixAlignment(str1[0:m-1], str2, del_cost, BLOSOM_matrix, edit_map)
		maps.Copy(edit_map, del_map)
	}

	prev_subs_key := edit_key{m - 1, n - 1}
	_, t = edit_map[prev_subs_key]
	if !t {
		subs_map := SuffixAlignment(str1[0:m-1], str2[0:n-1], del_cost, BLOSOM_matrix, edit_map)
		maps.Copy(edit_map, subs_map)
	}

	edit_map[key] = max(0,
		edit_map[prev_ins_key]+del_cost,
		edit_map[prev_del_key]+del_cost,
		edit_map[prev_subs_key]+BLOSOM_matrix[str1[m-1:]][str2[n-1:]])
	return edit_map
}

func SubstringAlignment(str1 string, str2 string, del_cost int, BLOSOM_matrix map[string]map[string]int) int {
	score := 0
	for _, value := range SuffixAlignment(str1, str2, del_cost, BLOSOM_matrix, make(map[edit_key]int)) {
		if value > score {
			score = value
		}
	}
	return score
}
