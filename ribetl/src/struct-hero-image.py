import argparse
import os,sys
import dotenv
sys.path.append("/usr/lib/python3/dist-packages")
from pymol import cmd, util


parser = argparse. ArgumentParser(description='generate a ribosome image')

parser. add_argument ("-s"    , "--structure"          , type= str ,required=True, help="RCSB ID of structure to process"                               )
parser. add_argument ("-env"  , "--dotenv_path"        , type= str ,required=True, help="Fallback dotenv path. Needed to locate the static files folder")
parser. add_argument ("-dest" , "--react_public_folder", type= str ,required=True, help="Static folder for react app, from where the images are ultimately served.")
# destination = f"/home/rxz/dev/ribosome.xyz-frontend.ts/public/ray_templates/_ray_{pdbid}.png"

args  = parser.parse_args()
pdbid = args.structure.upper()

dotenv.load_dotenv(dotenv_path=args.dotenv_path)
STATIC_ROOT = os.environ.get('STATIC_ROOT')

cmd.load(os.path.join(STATIC_ROOT, pdbid, f"{pdbid}.cif"))
cmd.reset()
cmd.spectrum('chain')
cmd.ray(500,500)
cmd.png(os.path.join(args.react_public_folder,f"_ray_{pdbid}.png"))
print('Saved {}'.format(os.path.join(args.react_public_folder,f"_ray_{pdbid}.png")))



