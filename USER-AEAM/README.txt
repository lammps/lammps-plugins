This folder contains a loader and support files to build the AEAM pair style
as plugin with CMake. This requires that LAMMPS was compiled with the
PLUGIN package enabled. To configure:

cmake -S . -B build -DLAMMPS_SOURCE_DIR=/path/to/lammps/src
cmake --build .

For more information on plugins please see: https://docs.lammps.org/Developer_plugins.html
