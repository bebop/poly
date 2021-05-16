package rhea

import (
	"database/sql"
	"errors"
)

/******************************************************************************

Database functions

These functions take in the rhea.rdf.gz dump file and import them into an sqlite
database. This database can be used locally for mapping functions, but more importantly,
it is used because it can enforce the relationships between each part of Rhea and alert
the user if the parser is failing to pick up on any of the relationships in the Rhea
database.

******************************************************************************/

// RheaSqlite inserts the Rhea database into an SQLite database with proper normalization for advanced queries.
func RheaSqlite(db *sql.DB, rhea Rhea) error {
	// Start transaction with database for insertion. This ensures if there are any problems, they are seamlessly rolled back
	tx, err := db.Begin()
	if err != nil {
		return err
	}

	schema := `
	-- Rhea tables --
	CREATE TABLE IF NOT EXISTS chebi (
		accession TEXT PRIMARY KEY,
		subclassof TEXT REFERENCES chebi(accession)
	);
	
	CREATE TABLE IF NOT EXISTS compound (
		id INT,
		accession TEXT PRIMARY KEY,
		position TEXT,
		name TEXT,
		htmlname TEXT,
		formula TEXT,
		charge TEXT,
		chebi TEXT REFERENCES chebi(accession),
		polymerizationindex TEXT,
		compoundtype TEXT NOT NULL CHECK(compoundtype IN ('SmallMolecule', 'Polymer', 'GenericPolypeptide', 'GenericPolynucleotide', 'GenericHeteropolysaccharide'))
	);
	
	CREATE TABLE IF NOT EXISTS reactivepart (
		id INT,
		accession TEXT PRIMARY KEY,
		name TEXT,
		htmlname TEXT,
		compound TEXT NOT NULL REFERENCES compound(accession)
	);
	
	CREATE TABLE IF NOT EXISTS reaction (
		id INT,
		directional BOOL NOT NULL DEFAULT false,
		accession TEXT PRIMARY KEY,
		status TEXT,
		comment TEXT,
		equation TEXT,
		htmlequation TEXT,
		ischemicallybalanced BOOL NOT NULL DEFAULT true,
		istransport BOOL NOT NULL DEFAULT false,
		ec TEXT,
		location TEXT
	);
	
	CREATE TABLE IF NOT EXISTS reactionside (
		accession TEXT PRIMARY KEY
	);
	
	CREATE TABLE IF NOT EXISTS reactionsidereaction (
		reaction TEXT NOT NULL REFERENCES reaction(accession),
		reactionside TEXT NOT NULL REFERENCES reactionside(accession),
		reactionsidereactiontype TEXT NOT NULL CHECK(reactionsidereactiontype IN ('substrateorproduct', 'substrate', 'product'))
	);
	
	CREATE TABLE IF NOT EXISTS reactionparticipant (
		compound TEXT REFERENCES compound(accession),
	        reactionside TEXT NOT NULL REFERENCES reactionside(accession),
	        contains INTEGER,
	        containsn BOOL NOT NULL DEFAULT false,
	        minus BOOL NOT NULL DEFAULT false,
	        plus BOOL NOT NULL DEFAULT false
	);`
	_, err = tx.Exec(schema)
	if err != nil {
		return errors.New("Failed to exec schema. Failed on: " + err.Error())
	}

	// First, insert ChEBIs and Compounds
	compoundKeys := make(map[string]bool)
	for _, compound := range rhea.Compounds {
		// Insert root chebi. Ie, what this current compound's subclass is
		_, err = tx.Exec("INSERT OR IGNORE INTO chebi(accession) VALUES (?)", compound.SubclassOfChEBI)
		if err != nil {
			tx.Rollback()
			return err
		}
		// Insert the chebi of the current compound. If it is already inserted, update the subclassification
		_, err = tx.Exec("INSERT INTO chebi(accession, subclassof) VALUES (?, ?) ON CONFLICT (accession) DO UPDATE SET subclassof = ?", compound.ChEBI, compound.SubclassOfChEBI, compound.SubclassOfChEBI)
		if err != nil {
			tx.Rollback()
			return err
		}
		// Insert the compound itself
		_, err = tx.Exec("INSERT INTO compound(id, accession, position, name, htmlname, formula, charge, chebi, polymerizationindex, compoundtype) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING", compound.ID, compound.Accession, compound.Position, compound.Name, compound.HTMLName, compound.Formula, compound.Charge, compound.ChEBI, compound.PolymerizationIndex, compound.CompoundType)
		if err != nil {
			tx.Rollback()
			return err
		}
		// If the compound isn't a small molecule or polymer, that means it would be a reactive part of a larger compound. So we add it in
		if (compound.CompoundType != "SmallMolecule") && (compound.CompoundType != "Polymer") {
			_, err = tx.Exec("INSERT INTO reactivepart(id, accession, name, htmlname, compound) VALUES (?, ?, ?, ?, ?)", compound.CompoundID, compound.CompoundAccession, compound.CompoundName, compound.CompoundHTMLName, compound.Accession)
			if err != nil {
				tx.Rollback()
				return err
			}
		}
		// Add compound.Access to the compoundKeys map
		compoundKeys[compound.Accession] = true
	}

	// Next, insert the ReactionSides and ReactionParticipants
	for _, reactionParticipant := range rhea.ReactionParticipants {
		// Insert ReactionSide, which is needed to insert the ReactionParticipant
		_, err = tx.Exec("INSERT INTO reactionside(accession) VALUES (?) ON CONFLICT DO NOTHING", reactionParticipant.ReactionSide)
		if err != nil {
			tx.Rollback()
			return err
		}
		// Insert the ReactionParticipants
		_, err = tx.Exec("INSERT INTO reactionparticipant(reactionside, contains, containsn, minus, plus, compound) VALUES (?, ?, ?, ?, ?, ?)", reactionParticipant.ReactionSide, reactionParticipant.Contains, reactionParticipant.ContainsN, reactionParticipant.Minus, reactionParticipant.Plus, reactionParticipant.Compound)
		if err != nil {
			tx.Rollback()
			return err
		}
	}

	// Insert the Reactions themselves
	for _, reaction := range rhea.Reactions {
		// Insert Reaction
		_, err = tx.Exec("INSERT INTO reaction(id, directional, accession, status, comment, equation, htmlequation, ischemicallybalanced, istransport, ec, location) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", reaction.ID, reaction.Directional, reaction.Accession, reaction.Status, reaction.Comment, reaction.Equation, reaction.HTMLEquation, reaction.IsChemicallyBalanced, reaction.IsTransport, reaction.Ec, reaction.Location)
		if err != nil {
			tx.Rollback()
			return err
		}

		// Insert ReactionsideReaction. Basically, these represent the substrates, products, or substratesOrProducts of a given reaction
		for _, substrate := range reaction.Substrates {
			_, err = tx.Exec("INSERT INTO reactionsidereaction(reaction, reactionside, reactionsidereactiontype) VALUES (?, ?, 'substrate')", reaction.Accession, substrate)
			if err != nil {
				tx.Rollback()
				return err
			}
		}
		for _, product := range reaction.Products {
			_, err = tx.Exec("INSERT INTO reactionsidereaction(reaction, reactionside, reactionsidereactiontype) VALUES (?, ?, 'product')", reaction.Accession, product)
			if err != nil {
				tx.Rollback()
				return err
			}
		}
		for _, substrateorproduct := range reaction.SubstrateOrProducts {
			_, err = tx.Exec("INSERT INTO reactionsidereaction(reaction, reactionside, reactionsidereactiontype) VALUES (?, ?, 'substrateorproduct')", reaction.Accession, substrateorproduct)
			if err != nil {
				tx.Rollback()
				return err
			}
		}

	}

	// Commit
	err = tx.Commit()
	if err != nil {
		return err
	}
	return nil
}
