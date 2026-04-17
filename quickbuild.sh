#/usr/bin/env bash

main() {

export DART=/Users/hkershaw/DART/Bugs/quads/DART
source "$DART"/build_templates/buildfunctions.sh

MODEL="none"
EXTRA="$DART"/models/template/threed_model_mod.f90
LOCATION="threed_sphere"

programs=(
)

serial_programs=(
test_quad_search
)

arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build DART
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
