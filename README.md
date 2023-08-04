This repository contains source code written for LAMMPS from
various sources that is not part of the LAMMPS distribution
for a number of reasons, but ported to be compatible with the
2 August 2023 LAMMPS release and converted to create plugins.

| Folders      | Origin of source code                                                       |
|--------------|-----------------------------------------------------------------------------|
| USER-AEAM    | https://github.com/psaidi/AEAM                                              |
| USER-BFIELD  | https://matsci.org/t/discrepency-beween-old-version-and-new-version/50014/5 |
| USER-REBOMOS | https://matsci.org/t/pair-rebomos/30503                                     |

As of the 22 December 2022 LAMMPS release, the USER-VCSGC package
code has been merged with upstream and the corresponding fix is now
available in the MC package. The plugin package has thus been removed.

To compile and configure all plugins

```
cmake -S . -B build -DLAMMPS_SOURCE_DIR=/path/to/lammps/src
cmake --build build
```

You can also perform the same commands inside the indivdual
plugin folders to configure and build only those plugin.

Please see README files in the individual folders for additional information.

For more information on plugins please see: https://docs.lammps.org/Developer_plugins.html

**WARNING:**
No attempt has been made by the LAMMPS developers to validate or
document the source code made available here.  All we have checked
is that it can run the included example inputs without crashing.
