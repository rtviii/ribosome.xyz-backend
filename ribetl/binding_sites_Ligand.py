from neo4j import GraphDatabase, Result
from dotenv import load_dotenv
import os
from typing import List, Tuple, TypedDict, Union, Callable
import sys
from neo4j.work.simple import Query
import pandas as pd
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
import json
from asyncio import run
from dotenv import load_dotenv


load_dotenv(dotenv_path='/home/rxz/dev/ribxz/.env')

STATIC_ROOT = os.getenv('STATIC_ROOT')
print("Got static root:", STATIC_ROOT)
print("Got static root:", STATIC_ROOT)

def _neoget(CYPHER_STRING:str)->Result:
    driver = GraphDatabase.driver(
    os.getenv( 'NEO4J_URI' ),
    auth=(os.getenv( 'NEO4J_USER' ),os.getenv( 'NEO4J_PASSWORD' )))

    def parametrized_query(tx, **kwargs):
        result:Result = tx.run(CYPHER_STRING, **kwargs)
        return result.values()

    with driver.session() as session:
        session.close()
        return session.read_transaction(parametrized_query)


class ResidueAsDict(TypedDict, total=False):
    resn       :  str
    strand_id  :  str
    resid      :  int
    struct     :  str
    banClass   :  str

class ResidueDict():
    """This is needed to hash and compare two residues inside a list"""

    def __init__(self, res:Residue):
        fid  =  list(res.get_full_id())
        self.model       :int  =  fid[1]
        self.resname     :str  =  fid[3][0]
        self.struct      :str  =  fid[0]
        self.chemicalName:str  =  res.get_resname()
        self.strand_id   :str  =  fid[2]
        self.residue_id  :int  =  [*fid[3]][1]
        self.banClass    :str  = ""

    def __eq__(self, other):
        return self.residue_id == other.residue_id and self.strand_id == other.strand_id

    def __hash__(self):
        return hash(( 'strand_id',self.strand_id,'residue_id', self.residue_id ))
    
    def toJSON(self):
            return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def asdict(self)->ResidueAsDict:
        return {
            "resn"       :  self.resname,
            "strand_id"  :  self.strand_id,
            "resid"      :  self.residue_id,
            "struct"     :  self.struct,
            "banClass"   :  self.banClass}

def getLigandResIds(ligchemid:str, struct: Structure)->List[Residue]:
    """Returns a list of dictionaries specifying each _ligand_ of type @ligchemid as a biopython-residue inside a given @struct."""
    """*ligchemids are of type https://www.rcsb.org/ligand/IDS"""
    ligandResidues: List[Residue] = list(filter(lambda x: x.get_resname() == ligchemid, list( struct.get_residues() )))
    return ligandResidues

async def matchStrandToClass(pdbid:str, strand_id:str)->Union[str, None]:
    """Request Ban nomenclature classes from the db given a protein's entity_poly_strand_id."""
    CYPHER="""match (r:RibosomeStructure{{_rcsb_id: "{}"}})-[]-(rp:RibosomalProtein{{entity_poly_strand_id:"{}"}})-[]-(n:NomenclatureClass)
    return n.class_id""".format(pdbid.upper(), strand_id)
    resp = _neoget(CYPHER)
    if len(resp) > 0:
        return resp[0]
    else:
        return None
def addBanClass(x:ResidueDict)->ResidueDict:
    """Tag a residue dictionary with the nomenclature of the strand it belongs to, if any is found"""
    banClass:str            =  run(matchStrandToClass(x.struct,x.strand_id))
    x.banClass = banClass
    return x

def getLigandNbrs(resids: List[Residue], struct:Structure)->List[ResidueDict]:
    """KDTree search the neighbors of a given list of residues(which constitue a ligand) 
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues. """

    ns   = NeighborSearch(list( struct.get_atoms() ))
    nbrs = []

    for r in resids:
        # a ligand consists of residues
        resatoms = r.child_list[0]
        #  each residue has an atom plucked at random
        for nbrresidues in ns.search(resatoms.get_coord(), 5,level='R'):
            # we grab all residues in radius around that atom and extend the list of neighbors with those
            nbrs.extend([nbrresidues])

    # Filter out the residues that constitute the ligand itself
    filtered = [] 
    for neighbor in nbrs:
        present = 0
        for constit in resids:
            if ResidueDict(constit)==ResidueDict( neighbor ):
                present = 1
        if present == 0:
            filtered.append(ResidueDict(neighbor))

    return [ * map(lambda x: addBanClass(x) ,  set(filtered) ) ]

def openStructutre(pdbid:str, cifpath: str)->Structure:
    return FastMMCIFParser(QUIET=True).get_structure(pdbid,cifpath)

def parseLigandNeighborhoods(pdbid:str,pathtostruct:str)->None:
    pdbid  =  pdbid.upper()

    if ("." in pdbid):
        print("Provide a PDB *ID*, not the file. The filepaths are defined in the .env.")
        return

    db_response     = _neoget("""match (l:Ligand)-[]-(r:RibosomeStructure{{rcsb_id:"{pdbid}"}}) 
    return {{struct: r.rcsb_id, ligs: collect({{ id:l.chemicalId, name: l.chemicalName }})}}""".format_map({ "pdbid":pdbid }))[0]

    class ResponseLigand(TypedDict):
        id:str;name:str
    presentLigands:List[ResponseLigand]

    if len(db_response)  == 0:
        print(f"No ligands for {pdbid} the DB. Exiting..")
        return
    else:
        print("Received ligands for {}: ".format(pdbid), db_response)
        presentLigands:List[ResponseLigand] = db_response[0]['ligs']

    # filtering the ions out
    dropIon       :Callable[[ResponseLigand], bool ]  =  lambda x: True if "ion" not in x['name'] else False
    presentLigands:List[ResponseLigand]                   =  [*filter(dropIon,presentLigands)]

    for ligand in presentLigands:

        # savepath = os.path.join(STATIC_ROOT, pdbid, 'LIGAND_{}.json'.format(x))
        savepath = os.path.join('LIGAND_{}.json'.format(ligand['id']))
        # ! These are a few of the things that turned out to be "problematic". Either whole structs not rendering or failing silently.
        # if pdbid in ['4U3N'] or x in ['A', 'OHX', ]:
        #     print("Skipping problematic {}".format(pdbid))
        #     continue
        # if pdbid in ['5TGM'] or x in ['A', 'OHX', 'LEU']:
        #     continue

        if os.path.exists(savepath):
            print(savepath, " already exists. Skipping rendering.")
            continue
        else:
            struct:Structure = openStructutre(pdbid, pathtostruct)
            # Parsing residues of ligand x
            ligand_as_residues    :List[Residue]      =  getLigandResIds(ligand['id'], struct )
            ligand_as_residuedicts:List[ResidueDict]  =  [ ResidueDict(res) for res in ligand_as_residues ]

            internal_residues   =  [addBanClass(ResidueDict(x)) for x in ligand_as_residues ]

            nbrs        =  getLigandNbrs(ligand_as_residues, struct)
            for nbr in nbrs:
                run(matchStrandToClass(nbr.struct,nbr.strand_id))

            ligandProfile  =  {
            'constituents': [ *map(lambda x: x.asdict(),internal_residues) ],
            'nbrs':         [ *map(lambda x: x.asdict(),nbrs) ]}
            with open(savepath, 'w') as json_file:
                json.dump(ligandProfile,json_file)
                print(f'Wrote report for {struct}/{ligand} to {savepath}')
    print('Done.')


if __name__ == "__main__":
    pdbid         =  sys.argv[1]
    pathtostruct  =  sys.argv[2]
    parseLigandNeighborhoods(pdbid,pathtostruct)
   