

Concerning the spread of ligand-prediction with respect to evolutionary divergence:

- Structural divergence between two proteins of the same class
- Sequence divergence between two proteins of the same class

What about ensemble-like learning?

- Structural divergence between two clusters of proteins(larger biomolecule) _given_ their  cumulative evolutionary divergence.






# Paper structure:




1  Introduction1.  
Resolution Revolution
	1. the spiel about universality (a lot of refs in proteovision) + cryoem, axes of variation
	2.  Current repositoriespdb as a current repository, its limitations/segue into nomenclature, mention other ribosomaldatabases

	3.  RiboXYZenables queries of the form [x]:
		Facilitate downstream comparative analysis,
		Provide easy andquick visualization tools,
		Facilitate search across structures and their components at differentresolution (whole structure, individual chain, ligand)1

2  Data infrastructure and processing

	A sentence or two about each of:
	structures, chains(proteins,rnas), ligand neighborhoods.
	Stress the utility of having a dedicated centralized repository as opposed to filtering the PDB:

	 -complete categories custom to the ribosome (protein, rna nomenclature, incoming: mitochondrial, ) 
	 - ease of access - enricheddata - connected graph —¿ more comprehensive querying - modularitydata  ETL  figure:   rough  flow  :   PDB.GQL  -¿  PFAM+Custom  Scripts  (Nomenclature-gen,ligands)  -¿  Neo4j  -¿  organization  of  the  frontend(think  about  how  to  homogenize  classes  in  thefrontend)3  Web server description3.1  Data Modules3.2  Tools4  Examples and ApplicationsRefer to the User Manual Tutorials•alignment and export of rna or relevant chains + figure•visualization•binding sites figureData Processing  Software used5  DiscussionDiscuss modularity and further applications future:•mitochondria•translation cycle stages•in-browser strucutre manip



Considering the scale of the ribosome as a biomolecule, it warrants

- of 177 591 structures, about 1500 of those are ribosomes  
- 144 338 of those are between 10 and  1000 residues in size:  all sorts of DNA, RNA fragments, oxygen transports, hydrolases, ligases, lyases, kinases, transferases, oxidoriductases, nucleases, synthases, isomerases, haloperoxidases, endonucleases, antibiotics and more exotic molecules.
- most of the ribosomes are in between 5000 and 25000 polymer residues-large
- the only things in the pdb that are effectively larger than ribosomes (>30k residues) are viruses.
- a panoply of things smaller: a lot of modularity to build on