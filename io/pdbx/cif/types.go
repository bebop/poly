package cif

// A CIF represents a complete CIF file, which
// is itself a collection of named and unordered
// DataBlocks.
type CIF struct {
	DataBlocks map[string]DataBlock
}

func NewCIF() CIF {
	return CIF{
		DataBlocks: make(map[string]DataBlock),
	}
}

// A DataBock is the highest-level component of a CIF
// and contains DataItems and SaveFrames.
type DataBlock struct {
	Name       string
	DataItems  map[string]any
	SaveFrames map[string]SaveFrame
}

func NewDataBlock(name string) DataBlock {
	return DataBlock{
		Name:       name,
		DataItems:  make(map[string]any),
		SaveFrames: make(map[string]SaveFrame),
	}
}

// A SaveFrame is a collection of data items that resides
// within a parent DataBlock.
type SaveFrame struct {
	Name      string
	DataItems map[string]any
}

func NewSaveFrame(name string) SaveFrame {
	return SaveFrame{
		Name:      name,
		DataItems: make(map[string]any),
	}
}

// A SpecialValue is a non-numeric, non-string value.
type SpecialValue string

const (
	// Inapplicable indicates the value is not applicable.
	Inapplicable SpecialValue = "."
	// Unknown indicates the value is unknown.
	Unknown SpecialValue = "?"
)
