from asyncio import gather
from dataclasses import dataclass
import json
import os
from pprint import pprint
import sys
from typing import Dict, List, Tuple
import dotenv
import numpy as np





# Void --> List[Path]
def gather_taxa():
    """"""
    structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
    profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))

    org_id_arrays = []

    for profile in profiles:
        org_id_arrays.append(profile_taxa(profile))
    return org_id_arrays



