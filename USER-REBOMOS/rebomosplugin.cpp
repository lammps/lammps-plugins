
#include "lammpsplugin.h"
#include "version.h"

#include "pair_rebomos.h"

using namespace LAMMPS_NS;

static Pair *pair_rebomos_creator(LAMMPS *lmp)
{
  return new PairREBOMoS(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // register rebomos pair style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "pair";
  plugin.name = "rebomos";
  plugin.info = "REBOMOS plugin pair style v1.0";
  plugin.author = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &pair_rebomos_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
