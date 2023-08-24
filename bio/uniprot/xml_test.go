package uniprot

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestIntListType_UnmarshalText(t *testing.T) {
	list := IntListType{}
	err := list.UnmarshalText([]byte("a"))
	assert.EqualError(t, err, `strconv.Atoi: parsing "a": invalid syntax`)
}
