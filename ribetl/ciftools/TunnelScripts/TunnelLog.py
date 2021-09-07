#!/usr/bin/env python3

import os,sys,numpy,json,math
from neo4j import Neo4jDriver
from dotenv import load_dotenv
import numpy as np

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')

import pandas as pd
import matplotlib.pyplot as plt
from typing import  Iterator, List
from dotenv import load_dotenv
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from ciftools.Structure import fetchStructure
from ciftools.Neoget import _neoget

Nucleotides  = ['A', 'U', 'G', 'C', 'T']

# AAs by assumed charge
AMINO_ACIDS={"ALA":0,'ARG':1,'ASN':0,'ASP':-1,'CYS':0,'GLN':0,'GLU':-1,'GLY':0,'HIS':0,'ILE':0,'LEU':0,'LYS':1,'MET':0,'PHE':0,'PRO':0,'SER':0,'THR':0,'TRP':0,'TYR':0,'VAL':0,'SEC':0,'PYL':0}

def get_CA_or(res:Residue)->Atom:

    """Returns the alpha carbon for the resiude if there is one, else the first atom in the list"""
    atoms:List[Atom] =list( res.get_atoms() )
    alphacarbons = list(filter(lambda atom: True if atom.get_name() == 'CA' else False ,atoms))
    return atoms[0] if len(alphacarbons) == 0 else alphacarbons[0]


class Log:

    """
    Log file itself is the interface to the csv files produced by MOLE.
    """
    def __init__(self, path:str)->None  :

        """Logging utility for keeping track of the ribosomal tunnels,
        cosuming and concatenating MOLE csv outputs,
        miscellaneous comments. \
        Initiate from a ${PROJECT_ROOT}/ciftools/TUNNELS/TUNNEL_LOG.csv
        _______________________________________________________________
        An assumption is made that the TUNNELS_LOG.csv resides in the 
        top-level directory which contains all the other tunnel files.
        """
        self.log :pd.DataFrame   = pd.read_csv(path)
        self.path = path
        self.__tunnels_path     = os.path.dirname(path)

    def _write(self)->None:
        self.log.to_csv(self.path,index=False)
        print("Has written successfully to {}".format(self.path))

    def all_structs(self):
        return self.log['pdbid'].tolist()
        
    def drop_column(self,colname)->None:
        self.log = self.log.drop([colname], axis=1)

    def add_column(self,colname:str)->None:

        if colname in self.log.columns.values:
            print("Column {} exists already".format(colname))
            return

        dflen = len(self.log['pdbid'])
        x     = np.zeros(dflen)
        self.log[colname]=x
        self._write()
        
    def get_struct(self, pdbid:str)->pd.DataFrame:
        pdbid = pdbid.upper()
        row   = self.log.loc[self.log['pdbid'] ==pdbid]

        if row.empty: 
            print("{} is empty.".format(pdbid))
            return pd.DataFrame()
        return row

    def update_struct(self, pdbid:str, dowrite:bool=False, **kwargs)->None:

        pdbid = pdbid.upper()
        row   = self.log.loc[self.log['pdbid'] ==pdbid]
        index = row.index

        if row.empty: 
            newrecord = pd.DataFrame({"pdbid": [ pdbid ],**kwargs})
            self.log  = self.log.append(newrecord, ignore_index=True)

        else:
            for item in kwargs.items():
                self.log.at[index, item[0]]= item[1]
        if dowrite:
            self._write()



class TunnelWalls:

    def __init__(self, pdbid:str, structure: Structure, mole_dataframe:pd.DataFrame) -> None:  
        self.mole_dataframe = mole_dataframe
        self.structure      = structure
        self.pdbid          = pdbid.upper()
        self.strands        = {}
        self.radius         = []
        self.ligands        = []
        self.rescount       = 0

        # self.rna                = {}
        # self.rps                = {}
        # self.nomenclatureMap    = {}
        # self.other              = []
        # self.adjacentRnaStrands = []
        # self.adjacentRPStrands  = []

    def _getDf(self):
        return self.mole_dataframe

    def __getProteinResidues(self):
        return self.rps

    def __getRnaResidues(self):
        return self.rna

    def __deprecated_addResidue(self, res:Residue)->None:
        parentStrand = res.get_parent().get_id()
        
        if res.get_resname() not in [AMINO_ACIDS.keys(), *Nucleotides] and res not in self.ligands:
            self.ligands.append(res)

        if parentStrand not in self.adjacentRnaStrands and parentStrand not in self.adjacentRPStrands:

            response = _neoget("""
            match (n {{entity_poly_strand_id:"{parentStrand}"}})-[]-(r:RibosomeStructure{{rcsb_id:"{pdbid}"}}) \
            return {{type: n.entity_poly_polymer_type, nomenclature:n.nomenclature}};""".format_map({
                "parentStrand": parentStrand,
                "pdbid"       : self.pdbid
            }))

            try:
                profile = response[0]
                if profile['type'] == 'RNA':
                    self.adjacentRnaStrands.append(parentStrand)
                    self.rna[parentStrand] = []

                if profile['type'] == 'Protein':
                    self.adjacentRPStrands.append(parentStrand)
                    self.rps            [parentStrand] = []
                    self.nomenclatureMap[parentStrand] = profile['nomenclature']
            except:
                pass
        if parentStrand in self.adjacentRnaStrands:
            if res in self.rna[parentStrand] :
                None
            else:
                self.rna[parentStrand].append(res)
        elif parentStrand in self.adjacentRPStrands:
            if res in self.rps[parentStrand]:
                None
            else:
                self.rps[parentStrand].append(res)
        print("added residue")
        self.rescount +=1
    def addResidue(self, res:Residue)->None:
        parentStrand = res.get_parent().get_id()
        
        if res.get_resname() not in [AMINO_ACIDS.keys(), *Nucleotides] and res not in self.ligands:
            self.ligands.append(res)
            return
            
        if parentStrand not in self.strands.keys():
            self.strands[parentStrand]=[]
        if res not in self.strands[parentStrand]:
            self.strands[parentStrand].append(res)

            # response = _neoget("""
            # match (n {{entity_poly_strand_id:"{parentStrand}"}})-[]-(r:RibosomeStructure{{rcsb_id:"{pdbid}"}}) \
            # return {{type: n.entity_poly_polymer_type, nomenclature:n.nomenclature}};""".format_map({
            #     "parentStrand": parentStrand,
            #     "pdbid"       : self.pdbid
            # }))

            # try:
            #     profile = response[0]

            #     if profile['type'] == 'RNA':
            #         self.adjacentRnaStrands.append(parentStrand)
            #         self.rna[parentStrand] = []

            #     if profile['type'] == 'Protein':
            #         self.adjacentRPStrands.append(parentStrand)
            #         self.rps            [parentStrand] = []
            #         self.nomenclatureMap[parentStrand] = profile['nomenclature']
            # except:
            #     pass
        # if parentStrand in self.adjacentRnaStrands:
        #     if res in self.rna[parentStrand] :
        #         None
        #     else:
        #         self.rna[parentStrand].append(res)
        # elif parentStrand in self.adjacentRPStrands:
        #     if res in self.rps[parentStrand]:
        #         None
        #     else:
        #         self.rps[parentStrand].append(res)
        # print("added residue")
        self.rescount +=1

    def consumeMoleDataframe(self, radius:float):

        """
        Takes a dataframe which must contain X,Y,Z columns assuming the centerline of the tunnel,
        although other columns can be present. Aimed at MOLE's mergd csv results.
        Iterates over each row and applies neighbor search on each focus, appends non-redundant residues
        to appropriate registry on the object. 
        """

        atoms        =  list(self.structure.get_atoms())
        ns           =  NeighborSearch(atoms,bucket_size=3)
        self.radius  =  radius


        def getVicinity(row):
            
            self.rescount+=1
            x = row['X'];
            y = row['Y'];
            z = row['Z']

            res:List[Residue.Residue] = ns.search(numpy.array([x,y,z]), radius,level='R')

            for nbr in res:
                self.addResidue(nbr)

        self.mole_dataframe.apply(getVicinity, axis=1)

    def get_ptc_residues(self, ptc_reslist:List[int])->List[Residue]:

        #ecoli=[ 2055 , 2056 , 2451 , 2452 , 2507 , 2506 ]
        def belongs_to_ptc(x:Residue):
            return int(x.get_id()[1]) in ptc_reslist

        PTC_residues = filter(belongs_to_ptc, [*self.structure.get_residues()]) 
        return [* PTC_residues ]
    
    def getResInfo(self,res): 

        res:Residue
        parent       : Chain.Chain = res.get_parent()
        parentStrand: str          = parent.get_id()
        resid        : int         = res.get_id()[1]
        resname      : str         = res.get_resname()
        rescoord                   = get_CA_or(res).get_coord().tolist()


        # nom   = nomMap[parentStrand] if polytype == "Protein" else None 
        islig = False if resname.upper() in [*Nucleotides, AMINO_ACIDS.keys()] else True

        return {
            "strand"  : parentStrand,
            "resid"   : resid,
            "resname" : resname,
            # "polytype": polytype,
            # "nom"     : nom,
            "isligand": islig,
            "rescoord": rescoord
            }

    def generateReport(self, write_to_path:str=""):
        """Consume Mole Dataframe first. Things are empty otherwise. Should be in the appropriate record."""
        # def PTC_coordinates():
        #     self.structure.

        wall_lining  =  {}
        # rnawall   =  {}


        for tpl in self.strands.items():
            wall_lining[tpl[0]] = [self.getResInfo(x) for x in tpl[1]]

        # for tpl in self.__getProteinResidues().items():
        #     protwall[ tpl[0] ] = [self.getResInfo(x,nomMap=self.nomenclatureMap, polytype="Protein") for x in tpl[1]]

        # for tpl in self.__getRnaResidues().items():
        #     rnawall[ tpl[0] ] = [self.getResInfo(x,nomMap=self.nomenclatureMap, polytype="RNA") for x in tpl[1]]

        report = {
              "pdbid"           : self.pdbid,
              "probeRadius"     : self.radius,
              "adjacent_strands": wall_lining,
            #   "ligands"         : presentLigands,
            # "rna"             : rnawall,
            # "proteins"        : protwall,
            # "nomMap"          : self.nomenclatureMap
            }

        if write_to_path != "":
            OUTPATH   =  write_to_path

            with open(OUTPATH, "w") as outfile:
                json.dump(report, outfile)
                print("Has written to path {}".format(OUTPATH))

        return report

