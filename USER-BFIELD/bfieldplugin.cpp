
#include "lammpsplugin.h"

#include "version.h"

#include "fix_bfield.h"

using namespace LAMMPS_NS;

static Fix *bfieldcreator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixBfield(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register bfield fix style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "bfield";
  plugin.info = "fix bfield plugin v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &bfieldcreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
