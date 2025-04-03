package main

import (
	"maps"
	"strings"
)

type edit_info struct {
	transcript string
	distance   int
}

type edit_key struct {
	key1 int
	key2 int
}

func EditDistance(str1 string, str2 string, edit_map map[edit_key]edit_info) map[edit_key]edit_info {
	m := len(str1)
	n := len(str2)
	if str1 == "" {
		edit_map[edit_key{0, n}] = edit_info{strings.Repeat("I", n), n}
		return edit_map
	}
	if str2 == "" {
		edit_map[edit_key{m, 0}] = edit_info{strings.Repeat("D", m), m}
		return edit_map
	}

	if edit_map[edit_key{m - 1, n}].transcript == "" {
		delete_map := EditDistance(str1[0:m-1], str2, edit_map)
		maps.Copy(edit_map, delete_map)
	}
	delete_transcript := edit_map[edit_key{m - 1, n}].transcript + "D"
	delete_distance := edit_map[edit_key{m - 1, n}].distance + 1

	if edit_map[edit_key{m, n - 1}].transcript == "" {
		insert_map := EditDistance(str1, str2[0:n-1], edit_map)
		maps.Copy(edit_map, insert_map)
	}
	insert_transcript := edit_map[edit_key{m, n - 1}].transcript + "I"
	insert_distance := edit_map[edit_key{m, n - 1}].distance + 1

	if edit_map[edit_key{m - 1, n - 1}].transcript == "" {
		substitute_map := EditDistance(str1[0:m-1], str2[0:n-1], edit_map)
		maps.Copy(edit_map, substitute_map)
	}
	substitute_transcript := edit_map[edit_key{m - 1, n - 1}].transcript
	substitute_distance := edit_map[edit_key{m - 1, n - 1}].distance
	if str1[m-1:] == str2[n-1:] {
		substitute_transcript += "M"
	} else {
		substitute_transcript += "R"
		substitute_distance += 1
	}

	distance := min(delete_distance, insert_distance, substitute_distance)

	switch distance {
	case delete_distance:
		edit_map[edit_key{m, n}] = edit_info{delete_transcript, distance}
		return edit_map
	case insert_distance:
		edit_map[edit_key{m, n}] = edit_info{insert_transcript, distance}
		return edit_map
	case substitute_distance:
		edit_map[edit_key{m, n}] = edit_info{substitute_transcript, distance}
		return edit_map
	default:
		return edit_map
	}
}

func EditDistanceBLOSUM(str1 string, str2 string, edit_map map[edit_key]edit_info, BLOSOM_matrix map[string]map[string]int) map[edit_key]edit_info {

	m := len(str1)
	n := len(str2)

	if str1 == "" {
		edit_map[edit_key{0, n}] = edit_info{strings.Repeat("I", n), (n * -8)}
		return edit_map
	}
	if str2 == "" {
		edit_map[edit_key{m, 0}] = edit_info{strings.Repeat("D", m), (m * -8)}
		return edit_map
	}

	if edit_map[edit_key{m - 1, n}].transcript == "" {
		delete_map := EditDistanceBLOSUM(str1[0:m-1], str2, edit_map, BLOSOM_matrix)
		maps.Copy(edit_map, delete_map)
	}
	delete_transcript := edit_map[edit_key{m - 1, n}].transcript + "D"
	delete_distance := edit_map[edit_key{m - 1, n}].distance - 8

	if edit_map[edit_key{m, n - 1}].transcript == "" {
		insert_map := EditDistanceBLOSUM(str1, str2[0:n-1], edit_map, BLOSOM_matrix)
		maps.Copy(edit_map, insert_map)
	}
	insert_transcript := edit_map[edit_key{m, n - 1}].transcript + "I"
	insert_distance := edit_map[edit_key{m, n - 1}].distance - 8

	if edit_map[edit_key{m - 1, n - 1}].transcript == "" {
		substitute_map := EditDistanceBLOSUM(str1[0:m-1], str2[0:n-1], edit_map, BLOSOM_matrix)
		maps.Copy(edit_map, substitute_map)
	}
	a := str1[m-1:]
	b := str2[n-1:]
	substitute_transcript := edit_map[edit_key{m - 1, n - 1}].transcript
	substitute_distance := edit_map[edit_key{m - 1, n - 1}].distance + BLOSOM_matrix[a][b]
	if a == b {
		substitute_transcript += "M"
	} else {
		substitute_transcript += "R"
	}

	distance := max(delete_distance, insert_distance, substitute_distance)

	switch distance {
	case delete_distance:
		edit_map[edit_key{m, n}] = edit_info{delete_transcript, distance}
		return edit_map
	case insert_distance:
		edit_map[edit_key{m, n}] = edit_info{insert_transcript, distance}
		return edit_map
	case substitute_distance:
		edit_map[edit_key{m, n}] = edit_info{substitute_transcript, distance}
		return edit_map
	default:
		return edit_map
	}
}

func GappedEditDistanceBLOSOM(str1 string, str2 string,
	edit_map map[edit_key]edit_info,
	aux_ins map[edit_key]edit_info,
	aux_del map[edit_key]edit_info,
	aux_subs map[edit_key]edit_info,
	gap_cost int,
	BLOSOM_matrix map[string]map[string]int) (map[edit_key]edit_info, map[edit_key]edit_info, map[edit_key]edit_info, map[edit_key]edit_info) {

	m := len(str1)
	n := len(str2)
	key := edit_key{m, n}

	if str1 == "" && str2 == "" {
		edit_map[key] = edit_info{"", 0}
		aux_ins[key] = edit_info{"", 0}
		aux_del[key] = edit_info{"", 0}
		aux_subs[key] = edit_info{"", 0}
		return edit_map, aux_ins, aux_del, aux_subs
	}
	if str1 == "" {
		edit_map[key] = edit_info{strings.Repeat("I", n), gap_cost}
		aux_ins[key] = edit_map[key]
		aux_del[key] = edit_info{"", 3 * gap_cost}
		aux_subs[key] = edit_info{"", 3 * gap_cost}
		return edit_map, aux_ins, aux_del, aux_subs
	}
	if str2 == "" {
		edit_map[key] = edit_info{strings.Repeat("D", m), gap_cost}
		aux_ins[key] = edit_info{"", 3 * gap_cost}
		aux_del[key] = edit_map[key]
		aux_subs[key] = edit_info{"", 3 * gap_cost}
		return edit_map, aux_ins, aux_del, aux_subs
	}

	prev_del_key := edit_key{m - 1, n}
	if edit_map[prev_del_key].transcript == "" {
		del_edit_map, del_aux_ins, del_aux_del, del_aux_subs := GappedEditDistanceBLOSOM(str1[0:m-1], str2, edit_map, aux_ins, aux_del, aux_subs, gap_cost, BLOSOM_matrix)
		maps.Copy(edit_map, del_edit_map)
		maps.Copy(aux_ins, del_aux_ins)
		maps.Copy(aux_del, del_aux_del)
		maps.Copy(aux_subs, del_aux_subs)
	}
	if aux_del[prev_del_key].distance < edit_map[prev_del_key].distance+gap_cost {
		aux_del[key] = edit_info{edit_map[prev_del_key].transcript + "D", edit_map[prev_del_key].distance + gap_cost}
	} else {
		aux_del[key] = edit_info{aux_del[prev_del_key].transcript + "D", aux_del[prev_del_key].distance}
	}

	prev_ins_key := edit_key{m, n - 1}
	if edit_map[prev_ins_key].transcript == "" {
		ins_edit_map, ins_aux_ins, ins_aux_del, ins_aux_subs := GappedEditDistanceBLOSOM(str1, str2[0:n-1], edit_map, aux_ins, aux_del, aux_subs, gap_cost, BLOSOM_matrix)
		maps.Copy(edit_map, ins_edit_map)
		maps.Copy(aux_ins, ins_aux_ins)
		maps.Copy(aux_del, ins_aux_del)
		maps.Copy(aux_subs, ins_aux_subs)
	}
	if aux_ins[prev_ins_key].distance < edit_map[prev_ins_key].distance+gap_cost {
		aux_ins[key] = edit_info{edit_map[prev_ins_key].transcript + "I", edit_map[prev_ins_key].distance + gap_cost}
	} else {
		aux_ins[key] = edit_info{aux_ins[prev_ins_key].transcript + "I", aux_ins[prev_ins_key].distance}
	}

	prev_subs_key := edit_key{m - 1, n - 1}
	if edit_map[prev_subs_key].transcript == "" {
		subs_edit_map, subs_aux_ins, subs_aux_del, subs_aux_subs := GappedEditDistanceBLOSOM(str1[0:m-1], str2[0:n-1], edit_map, aux_ins, aux_del, aux_subs, gap_cost, BLOSOM_matrix)
		maps.Copy(edit_map, subs_edit_map)
		maps.Copy(aux_ins, subs_aux_ins)
		maps.Copy(aux_del, subs_aux_del)
		maps.Copy(aux_subs, subs_aux_subs)
	}
	a := str1[m-1:]
	b := str2[n-1:]
	prev_subs_transcript := edit_map[prev_subs_key].transcript
	subs_distance := edit_map[prev_subs_key].distance + BLOSOM_matrix[a][b]
	if a == b {
		aux_subs[key] = edit_info{prev_subs_transcript + "M", subs_distance}
	} else {
		aux_subs[key] = edit_info{prev_subs_transcript + "R", subs_distance}
	}

	distance := max(aux_del[key].distance, aux_ins[key].distance, aux_subs[key].distance)

	switch distance {
	case aux_del[key].distance:
		edit_map[key] = edit_info{aux_del[key].transcript, distance}
		return edit_map, aux_ins, aux_del, aux_subs
	case aux_ins[key].distance:
		edit_map[key] = edit_info{aux_ins[key].transcript, distance}
		return edit_map, aux_ins, aux_del, aux_subs
	case aux_subs[key].distance:
		edit_map[key] = edit_info{aux_subs[key].transcript, distance}
		return edit_map, aux_ins, aux_del, aux_subs
	default:
		return edit_map, aux_ins, aux_del, aux_subs
	}
}
