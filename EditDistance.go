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
