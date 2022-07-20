
#include "lammpsplugin.h"
#include "version.h"

#include "fix_semigrandcanonical_mc.h"

using namespace LAMMPS_NS;

static Fix *fix_vcscg_creator(LAMMPS *lmp, int argc, char **argv)
{
    return new FixSemiGrandCanonicalMC(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register pace pair style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "sgcmc";
  plugin.info = "VCSGC plugin fix style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &fix_vcscg_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
