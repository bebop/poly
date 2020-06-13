package main

import (
	"lukechampine.com/blake3"
	"strings"
	"encoding/hex"
)

func least_rotation(s string) int {
	s += s
	k := 0
	f := make([]int, len(s))
	for i := range f {
		f[i] = -1
	}
	for j:=1; j<len(s); j++ {
		sj := s[j]
		i := f[j - k - 1]
		for (i != -1 && sj != s[k +i +1]) {
			if sj < s[k + i + 1] {
				k = j - i - 1
			}
			i = f[i]
		}
		if sj != s[k + i + 1] {
			if sj < s[k] {
				k = j
			}
			f[j - k] = -1
		} else {
			f[j - k] = i + 1
		}
	}
	return k
}

func Seqhash(sequence string, circular bool) string {
	s := strings.ToUpper(sequence)
	if circular {
		r := least_rotation(s)
		concat_seq := s+s
		s = concat_seq[r:r+len(s)]
	}
	b := blake3.Sum256([]byte(s))
	return hex.EncodeToString(b[:])
}

