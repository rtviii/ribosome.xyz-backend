--------------------------------------
General XML input shape:


<Tunnels>
  <Input Attributes="...">filename</Input>
  <ChargeSources>
    <ScalarGrid Attributes="..." />
    <AtomValues Attributes="..." />
    <RandomGrid Attributes="..." />
    <!--More optional Field elements-->
  </ChargeSources>
  <WorkingDirectory>path</WorkingDirectory>
  <Params>
    <Cavity Attributes="..." />
    <Tunnel Attributes="..." />
    <Filter>PatternQuery expression (see http://webchem.ncbr.muni.cz/Wiki for help)</Filter>
  </Params>
  <Export>
    <Formats Attributes="..." />
    <Types Attributes="..." />
    <Mesh Attributes="..." />
    <PyMol Attributes="..." />
    <Chimera Attributes="..." />
    <VMD Attributes="..." />
  </Export>
  <CustomVdw>
    <Radius Element="C" Value="1.75" />
    <!--More optional Radius elements-->
  </CustomVdw>
  <NonActiveParts>
    <Residue SequenceNumber="123" Chain="A" />
    <Query>PatternQuery expression (see http://webchem.ncbr.muni.cz/Wiki for help)</Query>
    <!--More optional Residue or Query elements-->
  </NonActiveParts>
  <CustomExits>
    <Exit>Point specification</Exit>
    <!--More optional Exit elements-->
  </CustomExits>
  <Origins Auto="0 or 1">
    <Origin>Point specification</Origin>
    <!--More optional Origin elements-->
  </Origins>
  <Paths>
    <!--Optional-->
    <Path>
      <Start>Point specification</Start>
      <End>Point specification</End>
    </Path>
    <!--More optional Path elements.-->
  </Paths>
</Tunnels>




Parameter types:
  [string] - a sequence of characters
  [real]   - a floating point number (1.23)
  [int]    - an integer (123)
  [bool]   - 0/1 or True/False

<Input> attributes: 
  SpecificChains
    Type [string] 
      One or more characters that specify which chains should participate in
      the computation. If empty or non present, all chains are loaded.
  ReadAllModels
    Type [bool], Default value = False 
      Determines whether to read all models from the PDB file. All models
      are read automatically for PDB assemblies (.pdb0 extension).


<ChargeSources> elements: 
  <AtomValues> attributes:
    Name
      Type [string] 
        Name of the field.
    Source
      Type [string] 
        Filename of the source. Supported formats: OpenDX.
    Method
      Type [string], Default value = NearestValue 
        Method used for charge computation [NearestValue, RadiusSum,
        RadiusSumDividedByDistance, RadiusMultiplicativeScale,
        RadiusAdditiveScale, KNearestSum,
        KNearestSumDividedByDistance, Lining, WholeStructure,
        RadiusResidueSum, RadiusResidueSumDividedByDistance,
        KNearestResidueSum, KNearestResidueSumDividedByDistance,
        AllResidues].
    Radius
      Type [real], Default value = 5 
        Radius used for the Radius* methods.
    K
      Type [int], Default value = 5 
        Number of neighbors used used for the KNearest* methods.
    IgnoreHydrogens
      Type [bool], Default value = False 
        Determines if to consider hydrogens for value assignment.
  <RandomGrid> attributes:
    Name
      Type [string] 
        Name of the field.
    MinValue
      Type [real], Default value = -1 
        Minimum value to generate.
    MaxValue
      Type [real], Default value = 1 
        Maximum value to generate.
  <ScalarGrid> attributes:
    Name
      Type [string] 
        Name of the field.
    Source
      Type [string] 
        Filename of the source. Supported formats: OpenDX.


<Params> elements: 
  <Cavity> attributes:
    IgnoreHETAtoms
      Type [bool], Default value = False 
        Allows to exclude the HET atoms from the calculation.
    IgnoreHydrogens
      Type [bool], Default value = False 
        Allows to exclude the hydrogen atoms from the calculation.
    InteriorThreshold
      Type [real], Default value = 1.25 
        Minimum radius of void inside the protein structure, so that
        the void would be considered a cavity.
    MinDepth
      Type [int], Default value = 8 
        Minimum cavity depth in the number of tetrahedrons.
    MinDepthLength
      Type [real], Default value = 5 
        Minimum depth of a cavity in angstroms.
    ProbeRadius
      Type [real], Default value = 3 
        Regulates level of detail of the molecular surface. Higher
        Probe Radius produces less detail.

  <Tunnel> attributes:
    AutoOriginCoverRadius
      Type [real], Default value = 10 
        The minimal distance between two auto origins in a given
        cavity.
    BottleneckRadius
      Type [real], Default value = 1.25 
        Minimal radius of a valid tunnel if BottleneckLength is 0.
    BottleneckTolerance
      Type [real], Default value = 0 
        Maximum length of a valid tunnel for which the radius is
        less than Bottleneck Radius.
    FilterBoundaryLayers
      Type [bool], Default value = False 
        Determines whether to remove layers with boundary residues
        from the tunnel.
    MaxAutoOriginsPerCavity
      Type [int], Default value = 5 
        The maximum number of automatically computed origins per
        cavity.
    MaxTunnelSimilarity
      Type [real], Default value = 0.9 
        Maximum degree of similarity between two tunnels before one
        tunnel is discarded.
    MinPoreLength
      Type [real], Default value = 0 
        Determines the minimal length (in ang) of a pore.
    MinTunnelLength
      Type [real], Default value = 0 
        Determines the minimal length (in ang) of a tunnel.
    OriginRadius
      Type [real], Default value = 5 
        If the user defined a tunnel start point, expand the search
        for tunnel start points to a sphere of radius.
    SurfaceCoverRadius
      Type [real], Default value = 10 
        Regulates the density of exit points tested at each outer
        boundary. Higher Surface Cover radius produces a lower
        density of exit points.
    UseCustomExitsOnly
      Type [bool], Default value = False 
        Only user defined exits are used for tunnel computation.
    WeightFunction
      Type [string], Default value = VoronoiScale 
        Determines the weight function used to compute channels
        [VoronoiScale, LengthAndRadius, Length, Constant].
  <Filter>
    A PatternQuery lambda expression that returns True to keep
    the given channel or False to discard it.

<Export> elements: 
  <Chimera> attributes:
    PDBId
      Type [string] 
        If this value is present, the option for downloading the
        structure is incorporated in the Chimera visualization
        script.
    SurfaceType
      Type [string], Default value = Spheres 
        Controls, if channels are displayed in Chimera as isosurface
        or as a set of spheres [Surface, Spheres].
  <Formats> attributes:
    ChargeSurface
      Type [bool], Default value = True 
        Controls if XML representation of surface with charges is
        created.
    Chimera
      Type [bool], Default value = False 
        Controls whether a Chimera script will be generated, for
        subsequent visualization in Chimera
    CSV
      Type [bool], Default value = False 
        Controls if CSV export of channel profiles created.
    JSON
      Type [bool], Default value = False 
        Controls if JSON output is created.
    Mesh
      Type [bool], Default value = False 
        Controls storing information about the mesh of detected
        tunnels, for subsequent visualization in PyMol or Jmol.
    PDBProfile
      Type [bool], Default value = False 
        Controls if channel profiles are exported in PDB format.
    PDBStructure
      Type [bool], Default value = False 
        Controls if channel residues (surrounding atoms) are
        exported in PDB format.
    PyMol
      Type [bool], Default value = False 
        Controls whether a PyMol script will be generated, for
        subsequent visualization in PyMol.
    VMD
      Type [bool], Default value = False 
        Controls whether a VMD script will be generated, for
        subsequent visualization in VMD.
  <Mesh> attributes:
    Compress
      Type [bool], Default value = False 
        Controls storing of mesh information in a GZip file archive.
    Density
      Type [real], Default value = 1.33 
        Level of detail of mesh surface. The higher the Mesh
        Density, the lower the level of detail in the visualization.
  <PyMol> attributes:
    ChargePalette
      Type [string], Default value = RedWhiteBlue 
        Determines the palette used for charge coloring
        [RedWhiteBlue, BlueWhiteRed].
    PDBId
      Type [string] 
        If this value is present, the option for downloading the
        structure is incorporated in the PyMOL visualization script.
    SurfaceType
      Type [string], Default value = Surface 
        Controls, if channels are displayed in PyMOL as isosurface
        or as a set of spheres [Surface, Spheres].
  <Types> attributes:
    Cavities
      Type [bool], Default value = True 
        Controls whether cavities will be exported.
    PoresAuto
      Type [bool], Default value = False 
        Controls whether 'auto' pores will be computed.
    PoresMerged
      Type [bool], Default value = False 
        Controls whether merged pores will be computed.
    PoresUser
      Type [bool], Default value = False 
        Controls whether user exit pores will be computed.
    Tunnels
      Type [bool], Default value = True 
        Controls whether tunnels will be computed.
  <VMD> attributes:
    PDBId
      Type [string] 
        If this value is present, the option for downloading the
        structure is incorporated in the VMD visualization script.
    SurfaceType
      Type [string], Default value = Spheres 
        Controls, if channels are displayed in VMD as isosurface or
        as a set of spheres [Surface, Spheres].


<Origin> attributes: 
  Auto
    Type [bool], Default value = False 
      Determines whether to (also) use automatically computed tunnel start
      points.

Points can be specified in 4 ways: 

- One or more residue elements
    <Residue SequenceNumber="123" Chain="A" InsertionCode=" " />
- One or more 3D points that 'snaps' to the closest residue
    <ResidueFromPoint X="1.0" Y="2.0" Z="3.0" />
- One or more 3D points
    <Point X="1.0" Y="2.0" Z="3.0" />
- PatternQuery expression (see http://webchem.ncbr.muni.cz/Wiki for help)
    <Query>expression</Query>
    Creates a point from each motif returned by the query.
The final start point is defined as the centroid of atomic centers
or defined points.