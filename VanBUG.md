# VANBUG December 9th 5:30, Vancouver BC

3 minute sprint

## ribosome.xyz

Ribosome xyz is a database application that provides a central point of access to organized ribosome data.
We get our data from the protein data bank and refine the subcomponents 

- the Protein Data Bank has a vast number of deposited models (~200k)
- of those only about 1500 are ribosomes

A given rbosome contains 50-80 odd proteins and some RNA chains, which are the most of the mass.

- ribosomes and viruses are about the largest (by a large margin) molecules in the PDB


# lack of specialization  : only polypeptides or polynucleotides

The problem we are addressing is that of specialization. At this moment PDB is only able to tell you whether a subchain is polypeptide or polynucleotide or neither and outside of that it's just arbitrary alphabetical ids. whereas different classes of protein and rna chains have very specific functions and should ideally be grouped accordingly. 



To this end we implement a classification system that does just that for RNA and proteins internal to the ribosome as well as  a some exogenous molecules: translation/elongation factors, mrna, multiple types of trna as well as various ligands and antibioitcs. 
For proteins -- we group them into classes of proteins proposed by some prominent people in the ribosomal community based on which PFAM family a given protein belongs to: in a nutshell. Arguably this can be more robust by pushing this upstream and checking every protein against templates on deposition.

Having these closed well-delineated categories allows us to look at things at a scale that is new. Some of these include:

- quantifying the conformational heterogeneity of the ribosome. Its parts are always in motion and MD is not available at this scale for a long time yet. 
- exit tunnel. 
- binding sites
- seq & struct homology 


## Limitations 

There are lots of limitations -- for example we don't have the mitoribosomes or the chloroplasts dialed in as a separate category, an interesting thing to look is stages of the translation cycles etc, but we are working to incorporate that soon. 

## Conclude

So yeah, you can download a set of homologs, like, let's say, ribosomal protein uL22 which modulates transcription in all species and run whatever analyses or learning on that dataset that you want to, no parsing and data-wrangling required. Same for every other category or multiple. Check it out at `ribosome.xyz`.
______________________________________________________

Arguably a refinement system like this should exist for every type of sizeable biological molecule. It's not clear how to do this given that a lot of practitioners specialize and it's hard to get infrastructure grants. 

I think it'd be awesome and in some sense inevitable that this will take on a new shape as primary data is being amassed at a blistering rate. Come chat! Suggestions and critique especially welcome!

# Babylon

I love the ribosome. Both as central, sophisticated mechanism and a segue to the more philosophical stuff like the origins of self-replication, but what i'm mostly thinking about these days is this problem of proliferation of types and formats that people always quip about in bioinformatics. 

There is an ELIXIR-funded initiative to integrate existing resources and PDB has developed a GraphQL endpoint recetly, but these are like attempts to mend something that's bursting at the seams and i'm wondering if there is a better way in 2021 to coordinate around this.


- the obvious fact being that there is 

- make illegal states unrepresentable

- prevalent languages being R, Python, Perl, JavaScript. Some purpose-built utilities are in C/C++( genome-assemblers, MD tools), but these are mostly data-terminuses: it either has a user interface or is provided as a binary so there isn't much to talk about in the way of exported data topology.

- maybe pick up 
- formats as opposed to types.
