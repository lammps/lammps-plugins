This repository contains source code written for LAMMPS from
various sources that is not part of the LAMMPS distribution
for a number of reasons, but ported to be compatible with the
23 June 2022 LAMMPS release and converted to create plugins.

Folders and origins
USER-AEAM: https://github.com/psaidi/AEAM
USER-VCSGC: https://gitlab.com/materials-modeling/vcsgc-lammps

To compile and configure all plugins

cmake -S . -B build -DLAMMPS_SOURCE_DIR=/path/to/lammps/src
cmake --build build

You can also perform the same commands inside the indivdual
plugin folders to configure and build only those plugin.

Please see README files in the individual folders for additional information.

For more information on plugins please see: https://docs.lammps.org/Developer_plugins.html
