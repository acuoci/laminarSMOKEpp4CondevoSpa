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
blockMesh

# Renumbering
renumberMesh -overwrite

# Split regions
splitMeshRegions -cellZonesOnly -overwrite 

# Remove unsplitted mesh
rm -r constant/polyMesh

# Check mesh
checkMesh -region fluid
checkMesh -region solid

# Manual operations
# Remember to change the boundary names and types

