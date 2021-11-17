import os,sys
from pprint import pprint
import json
from typing import List
from Bio import Align
from Bio.SeqRecord import SeqRecord
import dotenv
from pymol import cmd
import argparse
import glob
from ribetl.ciftools.bsite_mixed import BindingSite
from Bio.Align import MultipleSeqAlignment
import numpy as np
from Bio import pairwise2
import matplotlib.pyplot as plt
# import pylab as plt

def get_noms(a:dict):
	_ = []
	for i in a:
		for name in a[i]['nomenclature']:
			if name in _:
				print("DUPLICATE")
			else:
				_ = [*_, [name,len(a[i]['residues'])]  ]
	return _
sys.path.append('/home/rxz/dev/riboxyzbackend/')
dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get("STATIC_ROOT")
nomid       = sys.argv[1]
bsites      = []

profiles = []
for  _ in os.listdir(STATIC_ROOT):
	if len(_) > 4:
		continue
	else:
		profiles = [*profiles, os.path.join(STATIC_ROOT, _, "{}.json".format(_))]


def get_liglike_and_paths(profile_path:str):
	profile={};
	_ = []
	with open(profile_path, 'rb') as infile:
		profile  = json.load(infile)
		polymers=[]
		nonpolymers=[]

		if profile['proteins'] != None :
			polymers = [*polymers, *profile['proteins']]
		if profile['rnas'] != None :
			polymers = [*polymers, *profile['rnas']]

		try:nonpolymers = [*profile['ligands']]
		except:...


	for poly in polymers:
		if bool(poly['ligand_like']) == True:
			_.append({ 
				'description': poly['rcsb_pdbx_description'],
				'parent'     : profile['rcsb_id'],
				'chain'      : poly['auth_asym_id'],
				'path'       : os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(),"POLYMER_{}.json".format(poly['auth_asym_id'])) })

	for np in nonpolymers:
		if not "ion" in np['chemicalName'].lower():

			_.append({ 
				'description': np['chemicalName'],
				'parent'     : profile['rcsb_id'],
				'path'       : os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(),"LIGAND_{}.json".format(np['chemicalId'])) })
	return _
	

def get_matches_lig(ligpath:str,description:str, nomclass:str):

	try:
		with open(ligpath, 'rb') as lig_infile:
			data = json.load(lig_infile)
	except:
		return [ ]

	for chain in data:
		if nomclass in data[chain]['nomenclature']:
			# return data[chain]['sequence']
			return  [description, data[chain]['residues'] ]

	return []


residue_hits_all = []

for i in profiles:
	associated_liglike:List[dict] = get_liglike_and_paths(i)
	for ll in  associated_liglike:
		# matched_seqs = get_matches_lig(ll['path'], nomid)
		# sought_seqs  = [*sought_seqs, matched_seqs]

		matched_resids = get_matches_lig(ll['path'],ll['description'], nomid)
		if len(matched_resids) > 0:
			matched_resids = [matched_resids[0], *map(lambda _ : _['residue_id'], matched_resids[1])]
			residue_hits_all  = [*residue_hits_all, matched_resids]

# sought_seqs = [*filter(lambda x : x != None,sought_seqs)]
# pprint(sought_seqs)
# sought_seqs = [*map(lambda _ : SeqRecord(_), sought_seqs)]
# print(residue_hits)

candidates_len   = len(residue_hits_all)
candidates_names = [* map(lambda y : y[0], residue_hits_all)]
# print(candidates_len)
# pprint([ *enumerate(candidates_names) ])


# !-------------------※⋈------------- ∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷※⋈∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷--------⋮Ͽ⋈--------※⋈----Ͼ-----⋮∷-------------------------
conss ={
	'uS7'  : 'MARRRRAEVRQLQPDLVYGDVLVTAFINKIMRDGKKNLAARIFYDACKIIQEKTGQEPLKVFKQAVENVKPRMEVRSRRVGGANYQVPMEVSPRRQQSLALRWLVQAANQRPERRAAVRIAHELMDAAEGKGGAVKKKEDVERMAEANRAYAHYRW',
	'uL4'  : 'MELVLKDAQSALTVSETTFGRDFNEALVHQVVVAYAAGARQGTRAQKTRAEVTGSGKKPWRQKGTGRARSGSIKSPIWRSGGVTFAARPQDHSQKVNKKMYRGALKSILSELVRQDRLIVVEKFSVEAPKTKLLAQKLKDMALEDVLIITGELDENLFLAARNLHKVDVRDATGIDPVSLIAFDKVVMTADAVKQVEEMLA',
	'uL22' : 'METIAKHRHARSSAQKVRLVADLIRGKKVSQALDILTYTNKKAAVLVKKVLESAIANAEHNDGADIDDLKVTKIFVDEGPSMKRIMPRAKGRADRILKRTSHITVVVSDR',
	'uS4'  : 'MGRYIGPVCRLCRREGVKLYLKGERCYSPKCAMERRPYPPGQHGQKRARRPSDYAVRLREKQKLRRIYGISERQFRNLFEEASKKKGVTGSVFLGLLESRLDNVVYRLGFAVSRRQARQLVRHGHITVNGRRVDLPSYRVRPGDEIAVAEKSRNLELIRQNLEAMKGRKVGPWLSLDVEGMKGKFLRLPDREDLALPVNEQLVIEFYSR',
	'uL10' : "MALNLQDKQAIVAEVSEVAKGALSAVVADSRGVTVDKMTELRKAGREAGVYMRVVRNTLLRRAVEGTPFECLKDAFVGPTLIAYVTEHPGAAARLFKEFAKANAKFEVKAAAFEGELIPASQIDRLATLPTYEEAIARLMATMKEASAGKLVRTLAAVRDAKEAA",
	"bS20" : "MAQKKPKRNLSALKRHRQSLKRRLRNKAKKSAIKTLSKKAIQLAQEGKAEEALKIMRKAESLIDKAAKGSTLHKNAAARRKSRLMRKVRQLLEAAGAPLIGGGLSA",
	"bL9"  : "MKVILLEPLENLGDVGQVVDVKPGYARNYLLPRGLAVLATESNLKALEARIRAQAKRLAERKAEAERLKEILENLTLTIPVRAGETKIYGSVTAKDIAEALSRQHGVTIDPKRLALEKPIKELGEYVLTYKPHPEVPIQLKVSVVAQE",
	"bL12" : "MLPAAARPLWGPCLGLRAAAFRLARRQVPCVCAVRHMRSSGHQRCEALAGAPLDNAPKEYPPKIQQLVQDIASLTLLEISDLNELLKKTLKIQDVGLVPMGGVMSGAVPAAAAQEAVEEDIPIAKERTHFTVRLTEAKPVDKVKLIKEIKNYIQGINLVQAKKLVESLPQEIKANVAKAEAEKIKAALEAVGGTVVLE",
	"eS24" : "MNDTVTIRTRKFMTNRLLQRKQMVIDVLHPGKATVPKTEIRELAKMYKTTPDVIFVFGFRTHFGGGKTTGFGMIYDSLDYAKKNEPKHRLARHGLYEKKKTSRKQRKERKNRMKKVRGTAKANVGAGKKEPRG",
	"eS1"  : "MAVGKNKRLTKGGKKGAKKKVVDPFSKKDWYDVKAPAMFNIRNIGKTLVTRTQGTKIASDGLKGRVFEVSLADLQNDEVAFRKFKLITEDVQGKNCLTNFHGMDLTRDKMCSMVKKWQTMIEAHVDVKTTDGYLLRLFCVGFTKKRNNQIRKTSYAQHQQVRQIRKKMMEIMTREVQTNDLKEVVNKLIPDSIGKDIEKACQSIYPLHDVFVRKVKMLKKPKFELGKLMELHGEGSSSGKATGDETGAKVERADGYEPPVQESV",
	"eS4"  : "MRVKMHVKKGDTVLVASGKYKGRVGKVKEVLPKKYAVIVEGVNIVKKAVRVSPKYPQGGFIEKEAPLHASKVRPICPACGKPTRVRKKFLENGKKIRVCAKCGGALDTEE",
	"eS10" : "MLMPKEDRNKIHQYLFQEGVVVAKKDFNQAKHEEIDTKNLYVIKALQSLTSKGYVKTQFSWQYYYYTLTEEGVEYLREYLNLPEHIVPATYIQERNPTQRPQRRY",
	"eS26" : "PKKRASNGRNKKGRGHVKPVRCVNCSKSIPKDKAIKRMAIRNIVEAAAVRDLSEASVYPEYALPKTYNKLHYCVSCAIHARIVRVRSREDRKNRAPP",
	"eS30" : "AKVHGSLARAGKVKSQTPKVEKTEKPKKPKGRAYKRLLYTRRFVNVTLVNGKRRMNPGPSVQ",
	"RACK1": "ASNEVLVLRGTLEGHNGWVTSLATSAGQPNLLLSASRDKTLISWKLTGDDQKFGVPVRSFKGHSHIVQDCTLTADGAYALSASWDKTLRLWDVATGETYQRFVGHKSDVMSVDIDKKASMIISGSRDKTIKVWTIKGQCLATLLGHNDWVSQVRVVPNEKADDDSVTIISAGNDKMVKAWNLNQFQIEADFIGHNSNINTLTASPDGTLIASAGKDGEIMLWNLAAKKAMYTLSAQDEVFSLAFSPNRYWLAAATATGIKVFSLDPQYLVDDLRPEFAGYSKAAEPHAVSLAWSADGQTLFAGYTDNVIRVWQVMTAN"
}
# !-------※⋈------------------------- ∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷∷--------⋮Ͽ⋈--------※⋈----Ͼ-----⋮∷-------------------------


consensus_sequence = conss[nomid]
substrate = np.zeros(shape=(candidates_len,len(consensus_sequence)))
i         = 0

residue_hits_all = [*sorted(residue_hits_all, key=lambda x: x[0])]

fig, ax = plt.subplots(figsize=(24,10))
yheight = 4
# ax.text(-15, 4, 'unicode: Institut für Festkörperphysik', fontsize='4', color='blue')
for hit in residue_hits_all:

	# ax.text(-15, yheight, hit[0], fontsize='4', color='blue')
	# yheight+=2
	for resid in hit[1:]:

		if resid >= len(consensus_sequence):
			continue
		substrate[i,resid] = 1
	i+=1
	print(hit)

plt.yticks(range(0,candidates_len),[*map(lambda _: _[0], residue_hits_all)], fontsize="5")
ax.imshow(substrate,  aspect='auto')
plt.show()






# Align.MultipleSeqAlignment(sought_seqs)



# paromomycin : 50
# ery : 13
# kirromycin: 7
# viomycin : 12
# Apidaecin : 5
# listerin : 5
# neomycin :5 

# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0X 6ND6      --- ERY ---> 7aqc
# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0z 6S0X 6ND6 --- PAR ---> 4lfz
# 6XZB 6XZA 6XZ7 6OF1                               --- DI0 ---> 3j9w 