CREATE CONSTRAINT IF NOT EXISTS ON (ipro:InterProFamily   ) ASSERT ipro.family_id  IS UNIQUE;
CREATE CONSTRAINT IF NOT EXISTS ON (go  :GOClass          ) ASSERT go.class_id   IS UNIQUE;
CREATE CONSTRAINT IF NOT EXISTS ON (q   :RibosomeStructure) Assert q.rcsb_id    IS UNIQUE;
CREATE CONSTRAINT IF NOT EXISTS ON (pf  :PFAMFamily       ) assert pf.family_id  is unique;
CREATE CONSTRAINT IF NOT EXISTS ON (lig :Ligand           ) assert lig.chemicalId is unique;
CREATE CONSTRAINT IF NOT EXISTS ON (nc  :NomenclatureClass) assert nc.class_id   is unique;
