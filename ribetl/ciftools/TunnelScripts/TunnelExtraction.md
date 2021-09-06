Some of the parameters like scoop radius, structfile and executable paths figure in the `.env` which is in root.

### Whole procedure

For each structure:

1. Locate the [ptc](#PTC). Write barycenter(just averaging out the coordinates now hoping for something better than just a single nucleotide's position) to log.
2. Locate the [ constriction site ](#constriction). Write to Log.
3. Scoop a SCOOP_RADIUS around the ptc and save as `.pdb` file. (Mole doesn't work with .cifs): `scoop.py`
4. Iniate mole on the PTC.
5. Curate manually([ `tview.py` ](tview.py) to help).  
6. [ Generate ](#Report) "PDBID_TUNNEL_REPORT.json" based on the chosen csv tunnel. -> [ Tunnel Generation ](WallsReportGeneration.py)
7. Various [ plotting options ](PlotMultiple.py) in the works.




### Constriction

Contsriction site is obtained via interpolation between the closest points between uL22 and uL4(kdtrees) in [ `constriction_site` ](constriction_site.py). *Written to the Log.*

      - An assumption is made that constriction site is the unique closure of uL22 and uL4


### PTC

PTC is obtained from the following [ conserved residues ](#https://www.pnas.org/content/pnas/suppl/2010/09/23/1007988107.DCSupplemental/pnas.1007988107_SI.pdf?targetid=ST2) and initiates mole search.

The following positions are 100% and near-100% conserved in bacteria and archea respectively.

      2055
      2451
      2452
      2504
      2505
      2506
      2507





### Implementing now with the following mole setups for ecoli, Dec 22 after a go-ahead from Khanh:

```pythonh's feedback on tunnel cutoffs / the edge of the scoop.

      args['Points']             = origins
      args['ProbeRadius']        = "12"
      args['InteriorThreshold']  = "1.25"
      args['SurfaceCoverRadius'] = "10"
      args['OriginRadius']       = "5"
      args['BottleneckRadius']   = "1"
      args['exports']            = "t"
```

 
### The "biocuration step" (ha) involves just a comment to the pd dataframe for now with the following meanings:
Decided to stick with a single tunnel choice for now.
Assigining a negative value to the moletunnel of the dataframe to indicate that:

    - (0) mole was unable to find the correct tunnel
    - (-1) something is obviously blocking the tunnel 
    - (-2) the ribosome is not canonically assembled (PTC is fragmented/constriction site is far away from ptc)


There is a pretty urgent question of how to structure the interactions of MOLE's initiations with the common log, how much dependency there should be.
I'd like to keep PTC, constriction sites logged and accessible.

-----------------------

### Report Generation(Tunnel Walls)
