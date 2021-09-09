PRAGMA foreign_keys = ON;
CREATE TABLE codon (
  codon TEXT PRIMARY KEY,
  aa TEXT
);

CREATE TABLE seq (
  pos INT PRIMARY KEY
);

-- Weights are set on a per position basis for codon harmonization at a later point
CREATE TABLE weights (
  pos INTEGER REFERENCES seq(pos),
  codon TEXT NOT NULL,
  weight INTEGER,
  FOREIGN KEY(codon) REFERENCES codon(codon)
);

CREATE TABLE codonbias (
  gcbias TEXT CHECK(gcbias IN ('NA', 'GC', 'AT')),
  fromcodon TEXT NOT NULL,
  tocodon TEXT NOT NULL,
  FOREIGN KEY(fromcodon) REFERENCES codon(codon),
  FOREIGN KEY(tocodon) REFERENCES codon(codon)
);

CREATE TABLE suggestedfix (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  step INTEGER,
  start INT NOT NULL,
  end INT NOT NULL,
  gcbias TEXT,
  quantityfixes INTEGER,
  suggestiontype TEXT,
  FOREIGN KEY(start) REFERENCES seq(pos),
  FOREIGN KEY(end) REFERENCES seq(pos)
);

CREATE TABLE history (
  pos INTEGER,
  codon TEXT NOT NULL,
  step INT,
  suggestedfix INT,
  FOREIGN KEY(codon) REFERENCES codon(codon),
  FOREIGN KEY(suggestedfix) REFERENCES suggestedfix(id),
  FOREIGN KEY(pos) REFERENCES seq(pos)
      );