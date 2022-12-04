#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions      # Tutorial clean functions
#------------------------------------------------------------------------------

# Remove previous meshes
rm -r constant/polyMesh
rm -r constant/fluid/polyMesh
rm -r constant/solid/polyMesh

# Mesh conversion
fluentMeshToFoam  ../../Meshes/Mesh_80k.msh  -writeSets -writeZones
setsToZones -noFlipMap

# Renumbering
renumberMesh -overwrite

# Scaling (if needed)
# transformPoints -scale "(1e-3 1e-3 1e-3)" -overwrite

# Split regions
splitMeshRegions -cellZonesOnly -overwrite 

# Remove unsplitted mesh
rm -r constant/polyMesh

#Rename boundaries
sed -i 's/sideeast-fluid/sideEast/g' constant/fluid/polyMesh/boundary
sed -i 's/sidewest-fluid/sideWest/g' constant/fluid/polyMesh/boundary
sed -i 's/facenorth-fluid/faceNorth/g' constant/fluid/polyMesh/boundary
sed -i 's/facesouth-fluid/faceSouth/g' constant/fluid/polyMesh/boundary
sed -i 's/sideeast-solid/sideEast/g' constant/solid/polyMesh/boundary
sed -i 's/sidewest-solid/sideWest/g' constant/solid/polyMesh/boundary
sed -i 's/facenorth-solid/faceNorth/g' constant/solid/polyMesh/boundary
sed -i 's/facesouth-solid/faceSouth/g' constant/solid/polyMesh/boundary
sed -i 's/symmetry/symmetryPlane/g' constant/fluid/polyMesh/boundary
sed -i 's/symmetry/symmetryPlane/g' constant/solid/polyMesh/boundary


# Check mesh
checkMesh -region fluid
checkMesh -region solid

# Manual operations
# Remember to change the boundary names and types

