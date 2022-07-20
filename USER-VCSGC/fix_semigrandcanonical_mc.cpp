/* ----------------------------------------------------------------------
 * Parallel Monte-Carlo code for the semi-grandcanonical ensemble (SGC)
 * and the variance-constrained semi-grandcanonical ensemble (VC-SGC).
 *
 * See Sadigh et al., Phys. Rev. B 85, 184203 (2012) for a
 * description of the algorithm.
 *
 * Code author: Alexander Stukowski (stukowski@mm.tu-darmstadt.de)
 *
 * History:
 *
 * 27-Feb-09 - AS - Original version
 *
 * 11-Mar-09 - AS - The MC fix is now invoked on the POST_FORCE stage of an
 *                  MD step to avoid conflicts with other fixes that change
 *                  the box size.
 *
 * 12-Mar-09 - AS - Added support for an arbitrary number of atom types.
 *
 * 13-Mar-09 - AS - Particle velocity is rescaled after accepted swap
 *                  to conserve kinetic energy.
 *
 * 19-Jul-09 - AS - Made the MC routine compatible with Finnis-Sinclair type
 *                  potentials. Takes now into account that the electron density
 *                  at the swapped atom may change.
 *
 * 19-Jul-09 - AS - Added printing of parameter values to log file.
 *
 * 19-Jul-09 - AS - Added a serial VCSGC mode that doesn't use the second
 *                  rejection step for simulations on a single processor.
 *
 * 17-Nov-09 - AS - The MC routine can now be used with arbitrary potential styles
 *                  in serial mode. In this mode the generic total energy routine
 *                  is used to calculate the energy difference.
 *
 * 15-Apr-10 - AS - Added the Install.sh script and made some small changes
 *                  to make the module compatible with current version of LAMMPS.
 *
 * 17-Nov-10 - AS - Added debug code to check whether the random selection
 *                  of trial particles is evenly distributed in parallel simulations.
 *                  The individual trial count of each particle can be written to a dump file.
 *
 * 17-Nov-10 - AS - Fixed a bug in the printLog() function that led to
 *                  corrupted log files.
 *
 * 15-Feb-11 - AS - Fixed compilation error in the compute_vector() method by replacing
 *                  the max() function with inline code.
 *
 * 24-Sep-11 - AS - Adapted code to new interface of Error::one() function.
 *
 * 18-Sep-12 - AS - Renamed fix style from "sgc_mc" to "sgcmc".
 *
 * 07-Feb-13 - AS - Changed computeEnergyChangeGeneric() such that it no longer screws
 *                  up the stored atomic forces. Now it can actual be used for hybrid MD/MC simulations.
 *
 * 20-Jan-15 - AS - Updated call to Neighbor::request() to include the 'instance_me' identifier.
 *                  This change breaks backward-compatibility with LAMMPS versions prior to "20 Jan 2015".
 *
 * 10-Feb-15 - AS - Changed function pack_comm() to pack_forward_comm(), as required by recent LAMMPS version.
 *
 * 14-Apr-16 - AS - Restricted MC swaps to atoms in the fix group.
 *
------------------------------------------------------------------------- */

#include "fix_semigrandcanonical_mc.h"

#include "atom.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "comm.h"
#include "domain.h"
#include "universe.h"
#include "force.h"
#include "compute.h"
#include "memory.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "integrate.h"

#include "pair_eam.h"
#include "random_park.h"

#if CDEAM_MC_SUPPORT
 #include "pair_eam_cd.h"
#endif
#if TERSOFF_MC_SUPPORT
 #include "pair_tersoff.h"
#endif

#include <cstring>
#include <cstdlib>

using namespace LAMMPS_NS;
using namespace FixConst;

// This is for debugging purposes. The ASSERT() macro is used in the code to check
// if everything runs as expected.
#if SGCMC_DEBUG
    inline void my_noop() {}
    #define ASSERT(cond) ((!(cond)) ? error->one(FLERR, "Assertion failure.") : my_noop())
#else
    #define ASSERT(cond)
#endif

/*********************************************************************
 * Constructs the fix object and parses the input parameters
 * that control the Monte Carlo routine.
 *********************************************************************/
FixSemiGrandCanonicalMC::FixSemiGrandCanonicalMC(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), random(nullptr), localRandom(nullptr), neighborList(nullptr),
    samplingWindowUserSize(0), samplingWindowPosition(5), nAcceptedSwaps(0), nRejectedSwaps(0),
    kappa(0), serialMode(false), compute_pe(nullptr), pairEAM(nullptr)
{
    this->scalar_flag = 0;
    this->vector_flag = 1;
    this->extvector = 0;
    this->global_freq = 1;

    // Specifies the number of output fields this fix produces for thermo output.
    // It calculates the
    //    - Number of accepted trial moves
    //    - Number of rejected trial moves
    //    - Atom counts for each species.
    this->size_vector = 2 + atom->ntypes;

    // Let LAMMPS know the number of data values per atom to transfer in MPI communication.
    this->comm_forward = 4;
    this->comm_reverse = 3;

#if SGCMC_DEBUG
    this->peratom_flag = 1;
    this->size_peratom_cols = 0;
    this->peratom_freq = 1;
    this->create_attribute = 1;
    trialCounters = nullptr;
    grow_arrays(atom->nmax);
    ASSERT(trialCounters != nullptr || atom->nlocal == 0);
    memset(trialCounters, 0, sizeof(trialCounters[0]) * atom->nlocal);
#endif

    if(domain->triclinic)
        error->all(FLERR, "Fix sgcmc does not support non-orthogonal simulation boxes.");

    // Parse fix parameters from input file.
    if(narg < 6)
        error->all(FLERR, "Illegal fix sgcmc command. Not enough parameters.");

    // Counter for reading parameters.
    int iarg = 2;

    // Parse the number of MD timesteps to do between MC.
    iarg++;
    nevery_mdsteps = atoi(arg[iarg]);
    if (comm->me == 0)
      utils::logmesg(lmp, "  SGC - Number of MD timesteps: {}\n", nevery_mdsteps);
    if(nevery_mdsteps <= 0)
        error->all(FLERR, "Illegal fix sgcmc command. Invalid number of MD timesteps.");

    // Parse the fraction of atoms swaps attempted during each cycle.
    iarg++;
    swap_fraction = atof(arg[iarg]);
    if (comm->me == 0)
      utils::logmesg(lmp, "  SGC - Fraction of swap atoms: {}\n", swap_fraction);
    if(swap_fraction < 0 || swap_fraction > 1.0)
        error->all(FLERR, "Illegal fix sgcmc command. Invalid fraction of swap atoms.");

    // Parse temperature for MC.
    iarg++;
    double temperature = atof(arg[iarg]);
    if (comm->me == 0)
      utils::logmesg(lmp, "  SGC - Temperature: %f\n", temperature);
    if(temperature <= 0)
        error->all(FLERR, "Illegal fix sgcmc command. Temperature invalid.");
    double kb = 8.617343e-5;
    beta = 1.0 / ( kb * temperature );

    // Parse chemical potentials.
    iarg++;
    deltamu.resize(atom->ntypes + 1);
    deltamu[0] = 0.0;
    deltamu[1] = 0.0;
    if(atom->ntypes < 2)
        error->all(FLERR, "Illegal fix sgcmc command. Fix can only be used in simulations with at least two atom types.");
    for(int i=2; i<=atom->ntypes; i++, iarg++) {
        if(iarg >= narg)
            error->all(FLERR, "Illegal fix sgcmc command. Too few chemical potentials specified.");
        deltamu[i] = atof(arg[iarg]);
        if (comm->me == 0)
          utils::logmesg(lmp, "  SGC - Chemical potential of species {}: {}\n", i, deltamu[i]);
    }

    // Default values for optional parameters (where applicable).
    numSamplingWindowMoves = 8;
    seed = 324234;

    // Parse extra/optional parameters
    while(iarg < narg) {

        if(strcmp(arg[iarg], "randseed") == 0) {
            // Random number seed.
            if(iarg + 2 > narg) error->all(FLERR, "Illegal fix sgcmc command. Missing parameter after keyword 'randseed'.");
            seed = atoi(arg[iarg+1]);
            if (comm->me == 0)
              utils::logmesg(lmp, "  SGC - Random number seed: {}\n", seed);
            if(seed <= 0)
              error->all(FLERR, "Illegal fix sgcmc command. Random number seed must be positive.");
            iarg += 2;

        }
        else if(strcmp(arg[iarg], "window_moves") == 0) {
            // Parse number of window moves.
            if(iarg + 2 > narg) error->all(FLERR, "Illegal fix sgcmc command. Missing parameter after keyword 'window_moves'.");
            numSamplingWindowMoves = atoi(arg[iarg+1]);
            if (comm->me == 0)
              utils::logmesg(lmp, "  SGC - Number of sampling window moves: {}\n", numSamplingWindowMoves);
            if(numSamplingWindowMoves <= 0)
              error->all(FLERR, "Illegal fix sgcmc command. Invalid number of sampling window moves.");
            iarg += 2;

        }
        else if(strcmp(arg[iarg], "window_size") == 0) {
            // Parse sampling window size parameter.
            if(iarg + 2 > narg) error->all(FLERR, "Missing parameter after keyword 'window_size'.");
            samplingWindowUserSize = atof(arg[iarg+1]);
            if (comm->me == 0)
              utils::logmesg(lmp, "  SGC - Sampling window size: {}\n", samplingWindowUserSize);
            if(samplingWindowUserSize < 0.5 || samplingWindowUserSize > 1.0)
                error->all(FLERR, "Illegal fix sgcmc command. Sampling window size is out of range.");
            iarg += 2;
        }
        else if(strcmp(arg[iarg], "variance") == 0) {
            // Parse parameters for variance constraint ensemble.
            if(iarg + 1 + atom->ntypes > narg)
              error->all(FLERR, "Illegal fix sgcmc command. Too few parameters after keyword 'variance'.");
            iarg++;

            kappa = atof(arg[iarg]);
            if (comm->me == 0) utils::logmesg(lmp, "  SGC - Kappa: {}\n", kappa);
            if(kappa < 0)
              error->all(FLERR, "Illegal fix sgcmc command. Variance constraint parameter must not be negative.");
            iarg++;

            targetConcentration.resize(atom->ntypes + 1);
            targetConcentration[0] = 1.0;
            targetConcentration[1] = 1.0;
            for(int i=2; i<=atom->ntypes; i++, iarg++) {
                targetConcentration[i] = atof(arg[iarg]);
                targetConcentration[1] -= targetConcentration[i];
                if (comm->me == 0)
                        utils::logmesg(lmp, "  SGC - Target concentration of species {}: {}\n", i, targetConcentration[i]);
                if(targetConcentration[i] < 0 || targetConcentration[i] > 1.0)
                    error->all(FLERR, "Illegal fix sgcmc command. Target concentration is out of range.");
            }
            if (comm->me == 0)
              utils::logmesg(lmp, "  SGC - Target concentration of species 1: {}\n", targetConcentration[1]);
            if(targetConcentration[1] < 0)
                error->all(FLERR, "Illegal fix sgcmc command. Target concentration is out of range.");
        }
        else if(strcmp(arg[iarg], "serial") == 0) {
            // Switch off second rejection.
            serialMode = true;
            if (comm->me == 0)
              utils::logmesg(lmp, "  SGC - Using serial MC version without second rejection.\n");
            iarg++;

            if(comm->nprocs != 1)
                error->all(FLERR, "Illegal fix sgcmc command. Cannot use serial mode Monte Carlo in a parallel simulation.");
        }
        else error->all(FLERR, "Illegal fix sgcmc command. Unknown optional parameter.");
    }

    // Initialize random number generators.
    random = new RanPark(lmp, seed);
    localRandom = new RanPark(lmp, seed + universe->me);
}

/*********************************************************************
 * Destructor. Cleans up the random number generators.
 *********************************************************************/
FixSemiGrandCanonicalMC::~FixSemiGrandCanonicalMC()
{
    delete random;
    delete localRandom;

#if SGCMC_DEBUG
    memory->sfree(trialCounters);
#endif
}

#if SGCMC_DEBUG

/*********************************************************************
 * Allocate atom-based array.
 *********************************************************************/
void FixSemiGrandCanonicalMC::grow_arrays(int nmax)
{
    trialCounters = (double*)memory->srealloc(trialCounters, atom->nmax * sizeof(trialCounters[0]), "sgcmc:trialcounters");
    vector_atom = trialCounters;
}

/*********************************************************************
 * Copy values within local atom-based array.
 *********************************************************************/
void FixSemiGrandCanonicalMC::copy_arrays(int i, int j)
{
    trialCounters[j] = trialCounters[i];
}

/*********************************************************************
 * Initialize one atom's array values, called when atom is created.
 *********************************************************************/
void FixSemiGrandCanonicalMC::set_arrays(int i)
{
    trialCounters[i] = 0;
}

/*********************************************************************
 * Pack values in local atom-based array for exchange with another proc.
 *********************************************************************/
int FixSemiGrandCanonicalMC::pack_exchange(int i, double *buf)
{
    buf[0] = trialCounters[i];
    return 1;
}

/*********************************************************************
 * Unpack values in local atom-based array from exchange with another proc.
 *********************************************************************/
int FixSemiGrandCanonicalMC::unpack_exchange(int nlocal, double *buf)
{
    trialCounters[nlocal] = buf[0];
    return 1;
}

#endif

/*********************************************************************
 * The return value of this method specifies at which points the
 * fix is invoked during the simulation.
 *********************************************************************/
int FixSemiGrandCanonicalMC::setmask()
{
    // We want the MC routine to be called in between the MD steps.
    // We need the electron densities for each atom, so after the
    // EAM potential has computed them in the force routine is a good
    // time to invoke the MC routine.
    int mask = 0;
    mask |= POST_FORCE;
    mask |= POST_FORCE_RESPA;
    return mask;
}

/*********************************************************************
 * This gets called by the system before the simulation starts.
 *********************************************************************/
void FixSemiGrandCanonicalMC::init()
{
    // Make sure the user has defined only one Monte-Carlo fix.
    int count = 0;
    for(int i = 0; i < modify->nfix; i++)
        if(strcmp(modify->fix[i]->style,"sgcmc") == 0) count++;
    if(count > 1) error->all(FLERR, "More than one fix sgcmc defined.");

    // Save a pointer to the EAM potential.
    if((pairEAM = dynamic_cast<PairEAM*>(force->pair))) {
#if CDEAM_MC_SUPPORT
        if((pairCDEAM = dynamic_cast<PairEAMCD*>(pairEAM))) {
            if(pairCDEAM->cdeamVersion != 1)
                error->all(FLERR, "The sgcmc fix works only with the one-site concentration version of the CD-EAM potential.");
        }
#endif
    }
#if TERSOFF_MC_SUPPORT
    else if((pairTersoff = dynamic_cast<PairTersoff*>(force->pair))) {}
#endif
    else {
        if (comm->me == 0)
          utils::logmesg(lmp, "  SGC - Using naive total energy calculation for MC -> SLOW!\n");

        // Create a compute that will provide the total energy of the system.
        // This is needed by computeTotalEnergy().
        char* id_pe = (char*)"thermo_pe";
        int ipe = modify->find_compute(id_pe);
        compute_pe = modify->compute[ipe];
    }
    interactionRadius = force->pair->cutforce;
    if (comm->me == 0) utils::logmesg(lmp, "  SGC - Interaction radius: {}\n", interactionRadius);

    // This fix needs a full neighbor list.
    neighbor->add_request(this, NeighConst::REQ_FULL);

    // Count local number of atoms from each species.
    const int *type = atom->type;
    const int *mask = atom->mask;
    std::vector<int> localSpeciesCounts(atom->ntypes+1, 0);
    for(int i = 0; i < atom->nlocal; i++, ++type) {
        ASSERT(*type >= 1 && *type <= atom->ntypes);
        if(mask[i] & groupbit)
            localSpeciesCounts[*type]++;
    }

    // MPI sum to get global concentrations.
    speciesCounts.resize(atom->ntypes+1);
    MPI_Allreduce(&localSpeciesCounts.front(), &speciesCounts.front(), localSpeciesCounts.size(), MPI_INT, MPI_SUM, world);
}

/*********************************************************************
 * Assigns the requested neighbor list to the fix.
 *********************************************************************/
void FixSemiGrandCanonicalMC::init_list(int /*id*/, NeighList *ptr)
{
    this->neighborList = ptr;
}

/*********************************************************************
 * Called after the EAM force calculation during each timestep.
 * This method triggers the MC routine from time to time.
 *********************************************************************/
void FixSemiGrandCanonicalMC::post_force(int /*vflag*/)
{
    if((update->ntimestep % nevery_mdsteps) == 0)
        doMC();
}

/*********************************************************************
 * This routine does one full MC step.
 *********************************************************************/
void FixSemiGrandCanonicalMC::doMC()
{
    /// Reset energy variable to signal the energy calculation routine that
    /// it need to recompute the current total energy.
    totalPotentialEnergy = 0;

    // Allocate array memory.
    changedAtoms.resize(atom->nmax);

    // During the last MD timestep the EAM potential routine has computed the
    // electron densities for all atoms that belong to this processor.
    // They are stored in the rho array of the PairEAM class.
    // But computing the energy change caused by flipping one atom of this processor
    // might require the electron densities of atoms that belong to other processors.
    // So we first need to fetch those electron densities for our ghost atoms now.
    fetchGhostAtomElectronDensities();

    const int *mask = atom->mask;
#if SGCMC_DEBUG
    // This check is for debugging only! Can be safely removed.
    // Check the global concentration counter if it is still in sync on all nodes.
    const int *type = atom->type;
    std::vector<int> localSpeciesCounts(atom->ntypes+1, 0);
    std::vector<int> globalSpeciesCounts(atom->ntypes+1, 0);
    for(int i = 0; i < atom->nlocal; i++, ++type) {
        if(mask[i] & groupbit)
            localSpeciesCounts[*type]++;
    }
    MPI_Allreduce(&localSpeciesCounts.front(), &globalSpeciesCounts.front(), localSpeciesCounts.size(), MPI_INT, MPI_SUM, world);
    for(int i = 1; i <= atom->ntypes; i++) {
        // Note: If this test fails it might be an error in the code or
        //       because LAMMPS has lost atoms. This can happen for bad atom systems.
        ASSERT(globalSpeciesCounts[i] == speciesCounts[i]);
    }
#endif  // End of debugging code

    // Reset counters.
    int nAcceptedSwapsLocal = 0;
    int nRejectedSwapsLocal = 0;

    int oldSpecies, newSpecies;
    std::vector<int> deltaN(atom->ntypes+1, 0);         //< Local change in number of atoms of each species.
    std::vector<int> deltaNGlobal(atom->ntypes+1, 0);   //< Global change in number of atoms of each species.

    for(int i = 0; i < numSamplingWindowMoves; i++) {

        // Reset flag array that keeps track of changed per-atom quantities.
        std::fill(changedAtoms.begin(), changedAtoms.end(), false);

        // Position the sampling window within the node's boundaries.
        // By default the size of the sampling window is the size of the processor bounds minus two cutoff radii.
        // This ensures that changing atoms in the sampling windows of two adjacent processors cannot affect
        // the same atoms in the region between the two sampling windows.
        // For debugging purposes the sampling window can be chosen larger than the default size. Then it is
        // considered an 'oversize' window and we have to exchange atom information after each and
        // and every swap step, which is very slow.
        bool oversizeWindow = placeSamplingWindow();

        /// The number of times we want to swap an atom.
        int nDice = (int)(swap_fraction * numFixAtomsLocal / numSamplingWindowMoves);

        // This number must be synchronized with the other nodes. We take the largest
        // of all nodes and skip trial moves later.
        int largestnDice;
        MPI_Allreduce(&nDice, &largestnDice, 1, MPI_INT, MPI_MAX, world);

        // The probability to do one swap step.
        double diceProbability = (double)nDice / (double)largestnDice;

        // Inner MC loop that swaps atom types.
        for(int j = 0; j < largestnDice; j++) {

            double deltaE = 0;
            std::fill(deltaN.begin(), deltaN.end(), 0);
            int selectedAtom = -1, selectedAtomNL = -1;

            // This is only needed for debugging purposes:
            int selectedAtomDebug = -1;
            double deltaEDebug = 0;

            // As already said above, we have to do swap steps only with a certain probability
            // to keep nodes in sync.
            if(localRandom->uniform() <= diceProbability) {

                // Choose a random atom from the pool of atoms that are inside the sampling window.
                int index = (int)(localRandom->uniform() * (double)samplingWindowAtoms.size());
                ASSERT(index < samplingWindowAtoms.size());
                selectedAtomNL = samplingWindowAtoms[index];

                // Get the real atom index.
                ASSERT(selectedAtomNL < neighborList->inum);
                selectedAtom = neighborList->ilist[selectedAtomNL];
                ASSERT(selectedAtom < atom->nlocal);
                ASSERT(selectedAtom == selectedAtomNL);  // This assumption may be wrong.
                oldSpecies = atom->type[selectedAtom];

#if SGCMC_DEBUG
                trialCounters[selectedAtom]++;
#endif

                // Choose the new type for the swapping atom by random.
                if(atom->ntypes > 2) {
                    // Use a random number to choose the new species if there are three or more atom types.
                    newSpecies = (int)(localRandom->uniform() * (atom->ntypes-1)) + 1;
                    if(newSpecies >= oldSpecies) newSpecies++;
                }
                else {
                    // If there are only two atom types, then the decision is clear.
                    newSpecies = (oldSpecies == 1) ? 2 : 1;
                }
                ASSERT(newSpecies >= 1 && newSpecies <= atom->ntypes && newSpecies != oldSpecies);
                deltaN[oldSpecies] = -1;
                deltaN[newSpecies] = +1;

                // Compute the energy difference that swapping this atom would cost or gain.
                if(pairEAM) {
#if CDEAM_MC_SUPPORT
                    if(pairCDEAM)  // Concentration dependent EAM case:
                        deltaE = computeEnergyChangeCDEAM(selectedAtom, selectedAtomNL, oldSpecies, newSpecies);
                    else    // Standard EAM case:
#endif
                        deltaE = computeEnergyChangeEAM(selectedAtom, selectedAtomNL, oldSpecies, newSpecies);
                }
#if TERSOFF_MC_SUPPORT
                else if(pairTersoff) {
                    deltaE = computeEnergyChangeTersoff(selectedAtom, selectedAtomNL, oldSpecies, newSpecies);
                }
#endif
                else {
                    // Generic case:
                    deltaE = computeEnergyChangeGeneric(selectedAtom, oldSpecies, newSpecies);
                }

                // This is only needed for debugging purposes.
                // Save the values so they can later be checked.
                selectedAtomDebug = selectedAtom;
                deltaEDebug = deltaE;

                // Perform inner MC acceptance test.
                double dm = 0.0;
                if(serialMode && kappa != 0.0) {
                    for(int i = 2; i <= atom->ntypes; i++)
                        dm += (deltamu[i] + kappa / atom->natoms * (2.0 * speciesCounts[i] + deltaN[i])) * deltaN[i];
                }
                else {
                    for(int i = 2; i <= atom->ntypes; i++)
                        dm += deltamu[i] * deltaN[i];
                }
                double deltaB = -(deltaE + dm) * beta;
#if SGCMC_DEBUG
                std::cout << "Old species: " << oldSpecies << " new species: " << newSpecies << " deltaE: " << deltaE << " dm: " << dm << " deltaB: " << deltaB;
#endif
                if(deltaB < 0.0) {
                    if(deltaB < log(localRandom->uniform())) {
                        std::fill(deltaN.begin(), deltaN.end(), 0);
                        selectedAtom = -1;
                        deltaE = 0;
#if SGCMC_DEBUG
                        std::cout << " REJECTED";
#endif
                    }
                }

#if SGCMC_DEBUG
                cout << endl;
#endif
            }

#if 0
            // This is for debugging purposes only to check the energy change calculation routine.
            // Check the return value by calculating the energy change the slow way: difference between total energy before and after the swap.
            // The following code should be deactivated for full performance.
            // WARNING: Calling this function screws up the stored forces and will therefore change the time evolution of the system in
            // the MD simulation.
            double deltaECheck = computeEnergyChangeGeneric(selectedAtomDebug, oldSpecies, newSpecies);
            double deltaETotal, deltaECheckTotal;
            MPI_Allreduce(&deltaEDebug, &deltaETotal, 1, MPI_DOUBLE, MPI_SUM, world);
            MPI_Allreduce(&deltaECheck, &deltaECheckTotal, 1, MPI_DOUBLE, MPI_SUM, world);
            double posx = selectedAtomDebug >= 0 ? atom->x[selectedAtomDebug][0] : 0.0;
            double posy = selectedAtomDebug >= 0 ? atom->x[selectedAtomDebug][1] : 0.0;
            double posz = selectedAtomDebug >= 0 ? atom->x[selectedAtomDebug][2] : 0.0;
            printLog("Checking MC routine. DeltaE=%f DeltaE(check)=%f Atom: %i [%f %f %f] Old species: %i New species: %i\n",
                    deltaETotal, deltaECheckTotal, selectedAtomDebug, posx, posy, posz, oldSpecies, newSpecies);
            if(fabs(deltaECheckTotal - deltaETotal) > 1e-6) {   // Delta E must be equal to the total energy difference.
                error->one(FLERR, "Error in MC energy routine detected. Computed energy change deviates from correct value.");
            }
#endif      // End of debugging code


            if(kappa != 0.0 && serialMode == false) {

                // What follows is the second rejection test for the variance-constrained
                // semi-grandcanonical method.

                // MPI sum of total change in number of particles.
                MPI_Allreduce(&deltaN.front(), &deltaNGlobal.front(), deltaN.size(), MPI_INT, MPI_SUM, world);

                // Perform outer MC acceptance test.
                // This is done in sync by all processors.
                double A = 0.0;
                for(int i = 1; i <= atom->ntypes; i++) {
                    A += deltaNGlobal[i] * deltaNGlobal[i];
                    A += 2.0 * deltaNGlobal[i] * (speciesCounts[i] - (int)(targetConcentration[i] * atom->natoms));
                }
                double deltaB = -(kappa / atom->natoms) * A;
                if(deltaB < 0.0) {
                    if(deltaB < log(random->uniform())) {
                        std::fill(deltaN.begin(), deltaN.end(), 0);
                        std::fill(deltaNGlobal.begin(), deltaNGlobal.end(), 0);
                        selectedAtom = -1;
                    }
                }

                // Update global species counters.
                for(int i = 1; i <= atom->ntypes; i++)
                    speciesCounts[i] += deltaNGlobal[i];
            }
            else if(serialMode) {
                // Update the local species counters.
                for(int i = 1; i <= atom->ntypes; i++)
                    speciesCounts[i] += deltaN[i];
            }

            // Make accepted atom swap permanent.
            if(selectedAtom >= 0) {
                if(pairEAM) {
#if CDEAM_MC_SUPPORT
                    if(pairCDEAM) // Concentration dependent EAM case:
                        flipAtomCDEAM(selectedAtom, selectedAtomNL, oldSpecies, newSpecies);
                    else   // Standard EAM case:
#endif
                        flipAtomEAM(selectedAtom, selectedAtomNL, oldSpecies, newSpecies);
                }
                else
                    flipAtomGeneric(selectedAtom, oldSpecies, newSpecies);
                nAcceptedSwapsLocal++;
            }
            else {
                nRejectedSwapsLocal++;
            }

            // Update variable that keeps track of the current total energy.
            totalPotentialEnergy += deltaE;

            if(oversizeWindow) {
                // In case of an oversized sampling window we have to exchange the atom types and all other
                // per-atom quantities after each and every swap step. This is very slow and should only be used
                // for debugging purposes.
                communicateRhoAndTypes();
            }
        }

        // Finally the changed electron densities and atom types must be exchanged before
        // the sampling window is moved.
        if(!oversizeWindow)
            communicateRhoAndTypes();
    }

    // MPI sum total number of accepted/rejected swaps.
    MPI_Allreduce(&nAcceptedSwapsLocal, &nAcceptedSwaps, 1, MPI_INT, MPI_SUM, world);
    MPI_Allreduce(&nRejectedSwapsLocal, &nRejectedSwaps, 1, MPI_INT, MPI_SUM, world);

    // For (parallelized) semi-grandcanonical MC we have to determine the current concentrations now.
    // For the serial version and variance-constrained MC it has already been done in the loop.
    if(kappa == 0.0 && serialMode == false) {
        const int *type = atom->type;
        std::vector<int> localSpeciesCounts(atom->ntypes+1, 0);
        for(int i = 0; i < atom->nlocal; i++, ++type) {
            ASSERT(*type >= 1 && *type <= atom->ntypes);
            if(mask[i] & groupbit)
                localSpeciesCounts[*type]++;
        }
        MPI_Allreduce(&localSpeciesCounts.front(), &speciesCounts.front(), localSpeciesCounts.size(), MPI_INT, MPI_SUM, world);
    }
}

/*********************************************************************
 * Fetches the electron densities for the local ghost atoms
 * from the neighbor nodes.
 *********************************************************************/
void FixSemiGrandCanonicalMC::fetchGhostAtomElectronDensities()
{
    if(pairEAM) {
        // Transfer original EAM rho values.
        communicationStage = 1;
        comm->forward_comm(this);
    }
}

/*********************************************************************
 * Transfers the locally changed electron densities and atom
 * types to the neighbors.
 *********************************************************************/
void FixSemiGrandCanonicalMC::communicateRhoAndTypes()
{
    // Electron densities can have changed for real atoms as well as ghost atoms during the last MC step.
    // So we have to perform a forward and a reverse communication to keep everything in sync.
    // In the array changedAtoms we kept track of which rhos have been changed by the MC. This helps us
    // here to not overwrite values when doing the bidirectional exchange.

    if(pairEAM) {
        // Transfer changed electron densities of ghost atoms to the real atoms.
        communicationStage = 2;
        comm->reverse_comm(this);
    }

    // Transfer changed atom types and electron densities of the real atoms to the ghost atoms.
    communicationStage = 3;
    comm->forward_comm(this);
}

/*********************************************************************
 * This is for MPI communication with neighbor nodes.
 *********************************************************************/
int FixSemiGrandCanonicalMC::pack_forward_comm(int n, int* list, double* buf, int pbc_flag, int* pbc)
{
    int m = 0;
    if(communicationStage == 1) {
        // Send electron densities of local atoms to neighbors.
        ASSERT(pairEAM->rho != nullptr);
#if CDEAM_MC_SUPPORT
        if(pairCDEAM == nullptr) {
#endif
            for(int i = 0; i < n; i++)
                buf[m++] = pairEAM->rho[list[i]];
#if CDEAM_MC_SUPPORT
        }
        else {
            // In case of the CD-EAM model we have to send the RhoB values and D values as well.
            for(int i = 0; i < n; i++) {
                buf[m++] = pairCDEAM->rho[list[i]];
                buf[m++] = pairCDEAM->rhoB[list[i]];
                buf[m++] = pairCDEAM->D_values[list[i]];
            }
        }
#endif
    }
    else if(communicationStage == 3) {
        if(pairEAM) {
            // Send types and rhos of real atoms to the ghost atoms of the neighbor proc.
#if CDEAM_MC_SUPPORT
            if(pairCDEAM == nullptr) {
#endif
                for(int i = 0; i < n; i++) {
                    buf[m++] = atom->type[list[i]];
                    buf[m++] = pairEAM->rho[list[i]];
                }
#if CDEAM_MC_SUPPORT
            }
            else {
                // In case of the CD-EAM model we have to send the RhoB values and D values as well.
                for(int i = 0; i < n; i++) {
                    buf[m++] = atom->type[list[i]];
                    buf[m++] = pairCDEAM->rho[list[i]];
                    buf[m++] = pairCDEAM->rhoB[list[i]];
                    buf[m++] = pairCDEAM->D_values[list[i]];
                }
            }
#endif
        }
        else {
            // Generic potential case:
            for(int i = 0; i < n; i++) {
                buf[m++] = atom->type[list[i]];
            }
        }
    }
    else {
        ASSERT(false);
    }
    return m;
}

/*********************************************************************
 * This is for MPI communication with neighbor nodes.
 *********************************************************************/
void FixSemiGrandCanonicalMC::unpack_forward_comm(int n, int first, double* buf)
{
    if(communicationStage == 1) {
        // Receive electron densities of ghost atoms from neighbors.
        int last = first + n;
#if CDEAM_MC_SUPPORT
        if(pairCDEAM == nullptr) {
#endif
            for(int i = first; i < last; i++)
                pairEAM->rho[i] = *buf++;
#if CDEAM_MC_SUPPORT
        }
        else {
            // Also receive the partial densities and D values when using the CD-EAM model.
            for(int i = first; i < last; i++) {
                pairCDEAM->rho[i] = *buf++;
                pairCDEAM->rhoB[i] = *buf++;
                pairCDEAM->D_values[i] = *buf++;
            }
        }
#endif
    }
    else if(communicationStage == 3) {
        int last = first + n;
        if(pairEAM) {
            // Receive types and rhos of real atoms of the neighbor proc and assign them
            // to the local ghost atoms.
#if CDEAM_MC_SUPPORT
            if(pairCDEAM == nullptr) {
#endif
                for(int i = first; i < last; i++, buf += 2) {
                    atom->type[i] = (int)buf[0];
                    // We have to make sure that rhos changed locally do not get overridden by the rhos
                    // sent by the neighbor procs.
                    ASSERT(i < (int)changedAtoms.size());
                    if(!changedAtoms[i])
                        pairEAM->rho[i] = buf[1];
                }
#if CDEAM_MC_SUPPORT
            }
            else {
                // Also receive the partial densities and D values when using the CD-EAM model.
                for(int i = first; i < last; i++, buf += 4) {
                    atom->type[i] = (int)buf[0];
                    // We have to make sure that values changed locally do not get overridden by the values
                    // sent by the neighbor procs.
                    ASSERT(i < (int)changedAtoms.size());
                    if(!changedAtoms[i]) {
                        pairCDEAM->rho[i] = buf[1];
                        pairCDEAM->rhoB[i] = buf[2];
                        pairCDEAM->D_values[i] = buf[3];
                    }
                }
            }
#endif
        }
        else {
            // Generic potential case:
            for(int i = first; i < last; i++, buf += 1) {
                atom->type[i] = (int)buf[0];
            }
        }
    }
    else {
        ASSERT(false);
    }
}

/*********************************************************************
 * This is for MPI communication with neighbor nodes.
 *********************************************************************/
int FixSemiGrandCanonicalMC::pack_reverse_comm(int n, int first, double* buf)
{
    ASSERT(communicationStage == 2);
    int m = 0;

    // Send changed electron densities of ghost atoms to the real atoms of neighbor procs.
    ASSERT(pairEAM->rho != nullptr);
    int last = first + n;
#if CDEAM_MC_SUPPORT
    if(pairCDEAM == nullptr) {
#endif
        for(int i = first; i < last; i++)
            buf[m++] = pairEAM->rho[i];
#if CDEAM_MC_SUPPORT
    }
    else {
        // In case of the CD-EAM model we have to send the RhoB values and D values as well.
        for(int i = first; i < last; i++) {
            buf[m++] = pairCDEAM->rho[i];
            buf[m++] = pairCDEAM->rhoB[i];
            buf[m++] = pairCDEAM->D_values[i];
        }
    }
#endif
    return m;
}

/*********************************************************************
 * This is for MPI communication with neighbor nodes.
 *********************************************************************/
void FixSemiGrandCanonicalMC::unpack_reverse_comm(int n, int *list, double* buf)
{
    ASSERT(communicationStage == 2);

    // Received changed electron densities of ghost atoms of neighbor procs and assign them to our
    // real atoms.
    ASSERT(pairEAM->rho != nullptr);
#if CDEAM_MC_SUPPORT
    if(pairCDEAM == nullptr) {
#endif
        for(int i = 0; i < n; i++, buf++) {
            // We have to make sure that rhos changed locally do not get overridden by the rhos
            // sent by the neighbor procs.
            ASSERT(list[i] < (int)changedAtoms.size());
            if(!changedAtoms[list[i]])
                pairEAM->rho[list[i]] = *buf;
        }
#if CDEAM_MC_SUPPORT
    }
    else {
        for(int i = 0; i < n; i++, buf += 3) {
            // We have to make sure that rhos changed locally do not get overridden by the rhos
            // sent by the neighbor procs.
            ASSERT(list[i] < (int)changedAtoms.size());
            if(!changedAtoms[list[i]]) {
                pairCDEAM->rho[list[i]] = buf[0];
                pairCDEAM->rhoB[list[i]] = buf[1];
                pairCDEAM->D_values[list[i]] = buf[2];
            }
        }
    }
#endif
}

/*********************************************************************
 * Positions the sampling window inside the node's bounding box.
 *********************************************************************/
bool FixSemiGrandCanonicalMC::placeSamplingWindow()
{
    ASSERT(neighborList != nullptr);
    ASSERT(neighborList->inum == atom->nlocal);

    // By default the size of the sampling window is the size of the processor bounds minus two cutoff radii.
    // This ensures that changing atoms in the sampling windows of two adjacent processors cannot affect
    // the same atoms in the region between the two sampling windows.
    // For debugging purposes the sampling window can be chosen larger than the default size. Then it is
    // considered an 'oversize' window.
    bool oversizeWindow = false;

    // Align the sampling window to one of the 8 corners of the processor cell.
    double samplingWindowLo[3];
    double samplingWindowHi[3];
    double margin[3];
    for(int i = 0; i < 3; i++) {

        margin[i] = interactionRadius * 2.0;
        if(samplingWindowUserSize > 0.0) {
            margin[i] = (domain->subhi[i] - domain->sublo[i]) * (1.0 - samplingWindowUserSize);
            if(margin[i] < interactionRadius * 2.0)
                oversizeWindow = true;
        }

        double shift = (double)((samplingWindowPosition >> i) & 1) * margin[i];
        samplingWindowLo[i] = domain->sublo[i] + shift;
        samplingWindowHi[i] = domain->subhi[i] + shift - margin[i];

        // Check if processor cells are large enough.
        // Node bounds must be at least four times as large as the atom interaction radius.
        // Sampling window must be at least half as wise as the processor cell to cover the cell completely.
        if(samplingWindowHi[i] - samplingWindowLo[i] + 1e-6 < (domain->subhi[i] - domain->sublo[i]) * 0.5) {
            error->one(FLERR, "Per-node simulation cell is too small for fix sgcmc. Processor cell size must be at least 4 times cutoff radius.");
        }
    }
    // Increase counter by one.
    // Since we are only using the lower 3 bits of the integer value the alignment will
    // be the same after 8 iterations.
    samplingWindowPosition += 1;

    // Compile a list of atoms that are inside the sampling window.
    samplingWindowAtoms.resize(0);
    samplingWindowAtoms.reserve(atom->nlocal);
    numSamplingWindowAtoms = 0;
    numFixAtomsLocal = 0;

    ASSERT(atom->nlocal == neighborList->inum);
    const int *mask = atom->mask;
    for(int ii = 0; ii < neighborList->inum; ii++) {
        int i = neighborList->ilist[ii];
        if(mask[i] & groupbit) {
            numFixAtomsLocal++;
            const double* x = atom->x[i];
            // Is atom inside window region?
            if(x[0] >= samplingWindowLo[0] && x[0] < samplingWindowHi[0] &&
                x[1] >= samplingWindowLo[1] && x[1] < samplingWindowHi[1] &&
                x[2] >= samplingWindowLo[2] && x[2] < samplingWindowHi[2])
            {
                // Atoms within a distance of two times the interaction radius from the cell border
                // are less often inside the sampling window than atoms in the center of the node cell,
                // which are always inside the window.
                // We therefore have to increase their probability here to make them chosen
                // as often as the core atoms.
                int multiplicity = 1;
                for(int k=0; k < 3; k++) {
                    if(x[k] < domain->sublo[k] + margin[k] ||
                       x[k] > domain->subhi[k] - margin[k])
                        multiplicity *= 2;
                }

                for(int m = 0; m < multiplicity; m++)
                    samplingWindowAtoms.push_back(ii);

                numSamplingWindowAtoms++;
            }
        }
    }

    return oversizeWindow;
}

/*********************************************************************
 * Calculates the change in energy that swapping the given
 * atom would produce. This routine is for the standard EAM potential.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * flipAtomNL [in]
 *   This specifies the atom to be swapped. It's an index into the neighbor list.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new species of the atom. The atom's type is not changed by this routine. It only computes the induced energy change.
 *
 * Return value:
 *   The expected change in total potential energy.
 *********************************************************************/
double FixSemiGrandCanonicalMC::computeEnergyChangeEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies)
{
    double p;
    int m;
    double const* rho = pairEAM->rho;
    double* coeff;
    double new_total_rho_i = 0.0;
    double deltaE = 0.0;

    // Calculate change of electron density at the surrounding
    // sites induced by the swapped atom. Then calculate the change of embedding energy for each neighbor atom.
    // Also recalculate the total electron density at the site of the swapped atom.

    double xi = atom->x[flipAtom][0];
    double yi = atom->x[flipAtom][1];
    double zi = atom->x[flipAtom][2];

    // Loop over all neighbors of the selected atom.
    ASSERT(flipAtomNL < neighborList->inum);
    int* jlist = neighborList->firstneigh[flipAtomNL];
    int jnum = neighborList->numneigh[flipAtomNL];
    for(int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];

        double delx = xi - atom->x[j][0];
        double dely = yi - atom->x[j][1];
        double delz = zi - atom->x[j][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        if(rsq >= pairEAM->cutforcesq) continue;

        int jtype = atom->type[j];
        double r = sqrt(rsq);

        p = r * pairEAM->rdr + 1.0;
        m = static_cast<int>(p);
        m = MIN(m, pairEAM->nr - 1);
        p -= m;
        p = MIN(p, 1.0);

        // Calculate change of pair energy ij.
        coeff = pairEAM->z2r_spline[pairEAM->type2z2r[oldSpecies][jtype]][m];
        double oldz2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = pairEAM->z2r_spline[pairEAM->type2z2r[newSpecies][jtype]][m];
        double newz2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        deltaE += (newz2 - oldz2) / r;

        // Calculate change of electron density at site j.
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[oldSpecies][jtype]][m];
        double oldrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[newSpecies][jtype]][m];
        double newrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        double delta_rho = newrho_contr - oldrho_contr;

        // Sum total rho at site of swapped atom.
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[jtype][newSpecies]][m];
        new_total_rho_i += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        // Calculate old embedding energy of atom j.
        p = rho[j] * pairEAM->rdrho + 1.0;
        m = static_cast<int>(p);
        m = MAX(1, MIN(m, pairEAM->nrho - 1));
        p -= m;
        p = MIN(p, 1.0);
        coeff = pairEAM->frho_spline[pairEAM->type2frho[jtype]][m];
        double oldF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        // Calculate new embedding energy of atom j.
        p = (rho[j] + delta_rho) * pairEAM->rdrho + 1.0;
        m = static_cast<int>(p);
        m = MAX(1, MIN(m, pairEAM->nrho - 1));
        p -= m;
        p = MIN(p, 1.0);
        coeff = pairEAM->frho_spline[pairEAM->type2frho[jtype]][m];
        double newF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        deltaE += newF - oldF;
    }

    // Compute the change in embedding energy of the changing atom.
    p = rho[flipAtom] * pairEAM->rdrho + 1.0;
    m = static_cast<int>(p);
    m = MAX(1, MIN(m, pairEAM->nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    coeff = pairEAM->frho_spline[pairEAM->type2frho[oldSpecies]][m];
    double oldF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

    p = new_total_rho_i * pairEAM->rdrho + 1.0;
    m = static_cast<int>(p);
    m = MAX(1, MIN(m, pairEAM->nrho - 1));
    p -= m;
    p = MIN(p, 1.0);
    coeff = pairEAM->frho_spline[pairEAM->type2frho[newSpecies]][m];
    double newF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

    deltaE += newF - oldF;

    return deltaE;
}

#if CDEAM_MC_SUPPORT

/*********************************************************************
 * Calculates the change in energy that swapping the given
 * atom would produce.
 * This routine is for the concentration-dependent EAM potential.
 * It has to do some extra work in comparison to the standard EAM routine above.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * flipAtomNL [in]
 *   This specifies the atom to be swapped. It's an index into the neighbor list.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new species of the atom. The atom's type is not changed by this routine. It only computes the induced energy change.
 *
 * Return value:
 *   The expected change in total potential energy.
 *********************************************************************/
double FixSemiGrandCanonicalMC::computeEnergyChangeCDEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies)
{
    ASSERT(pairCDEAM != nullptr);  // Make sure we have a CD-EAM potential in use.

    double p;
    int m;
    double* rho = pairEAM->rho;
    double* coeff;
    double new_total_rho_i = 0.0;
    double new_total_rhoB_i = 0.0;

    // The energy change to calculate.
    double deltaE = 0.0;

    // Calculate each change of electron density (and partial density) at the
    // surrounding sites induced by the swapped atom. Also calculate the change of pair interaction energy.
    // Then calculate the change of embedding energy for each neighbor atom.

    double xi = atom->x[flipAtom][0];
    double yi = atom->x[flipAtom][1];
    double zi = atom->x[flipAtom][2];

    /// The change in atom i's D value.
    double deltaD_i = 0.0;

    // Loop over all neighbors of the selected atom.
    int* jlist = neighborList->firstneigh[flipAtomNL];
    int jnum = neighborList->numneigh[flipAtomNL];
    for(int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];

        double delx = xi - atom->x[j][0];
        double dely = yi - atom->x[j][1];
        double delz = zi - atom->x[j][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        if(rsq >= pairEAM->cutforcesq) continue;

        int jtype = atom->type[j];
        double r = sqrt(rsq);

        p = r * pairEAM->rdr + 1.0;
        m = static_cast<int>(p);
        m = MIN(m, pairEAM->nr - 1);
        p -= m;
        p = MIN(p, 1.0);

        // Calculate change of electron density at site j.
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[oldSpecies][jtype]][m];
        double oldrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[newSpecies][jtype]][m];
        double newrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        double delta_rho = newrho_contr - oldrho_contr;

        // Sum total rho at site of swapped atom.
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[jtype][newSpecies]][m];
        double rho_ji = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        new_total_rho_i += rho_ji;
        if(jtype == pairCDEAM->speciesB)
            new_total_rhoB_i += rho_ji;

        // Determine change of partial electron density at site j.
        // It increases if the atom i becomes a B atom and it decreases if it was an B atom.
        double delta_rhoB_j = 0;
        if(newSpecies == pairCDEAM->speciesB) delta_rhoB_j = newrho_contr;
        else if(oldSpecies == pairCDEAM->speciesB) delta_rhoB_j = -oldrho_contr;

        // Now we can calculate the new concentration x_j at site j.
        double new_x_j = (pairCDEAM->rhoB[j] + delta_rhoB_j) / (rho[j] + delta_rho);
        // Calculate the old concentration x_j at site j as well.
        double old_x_j = pairCDEAM->rhoB[j] / rho[j];

        // Calculate change of pair energy ij.
        coeff = pairEAM->z2r_spline[pairEAM->type2z2r[oldSpecies][jtype]][m];
        double oldz2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = pairEAM->z2r_spline[pairEAM->type2z2r[newSpecies][jtype]][m];
        double newz2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        /// Was the old i-j interaction concentration dependent?
        if((jtype == pairCDEAM->speciesA && oldSpecies == pairCDEAM->speciesB)
                || (jtype == pairCDEAM->speciesB && oldSpecies == pairCDEAM->speciesA)) {

            // The old pair interaction was concentration dependent.
            // Please note that it must now become concentration independent that
            // either the A or the B atom has changed to another species.
            // We here require that only the A-B interactions are concentration dependent.

            // Add the new static pair interaction term.
            deltaE += newz2 / r;

            // The D values of i and j decrease due to the removed concentration dependent cross interaction.
            double deltaD_j = -oldz2 / r;
            deltaD_i += deltaD_j;

            // Since the concentration x_j at site j changes, the mixed pair interaction between atom j and
            // other atoms k is also affected.
            // The energy change of site j due to the cross term is: (new_h*new_D - old_h*old_D)/2
            double old_h_j = pairCDEAM->evalH(old_x_j);
            double new_h_j = pairCDEAM->evalH(new_x_j);
            deltaE += 0.5 * (new_h_j * (pairCDEAM->D_values[j] + deltaD_j) - old_h_j * pairCDEAM->D_values[j]);
        }
        else {

            // The old pair interaction was not concentration dependent. Now it might have
            // become dependent. Do check:
            if((jtype == pairCDEAM->speciesA && newSpecies == pairCDEAM->speciesB)
                || (jtype == pairCDEAM->speciesB && newSpecies == pairCDEAM->speciesA)) {

                // The new pair interaction is concentration dependent. It's an AB interaction.

                // Subtract old static pair interaction.
                deltaE -= oldz2 / r;

                // The D values of i and j increase due to the created AB cross interaction.
                double deltaD_j = newz2 / r;
                deltaD_i += deltaD_j;

                // Since the concentration x_j at site j changes, the mixed pair interaction between atom j and
                // other atoms k is also affected.
                // The energy change of site j due to the cross term is: (new_h*new_D - old_h*old_D)/2
                double old_h_j = pairCDEAM->evalH(old_x_j);
                double new_h_j = pairCDEAM->evalH(new_x_j);
                deltaE += 0.5 * (new_h_j * (pairCDEAM->D_values[j] + deltaD_j) - old_h_j * pairCDEAM->D_values[j]);
            }
            else {

                // The pair interaction stays concentration independent.
                // This is like standard EAM.
                deltaE += (newz2 - oldz2) / r;

                double old_h_j = pairCDEAM->evalH(old_x_j);
                double new_h_j = pairCDEAM->evalH(new_x_j);
                deltaE += 0.5 * (new_h_j * (pairCDEAM->D_values[j]) - old_h_j * pairCDEAM->D_values[j]);
            }
        }

        // Calculate old embedding energy of atom j.
        p = rho[j] * pairEAM->rdrho + 1.0;
        m = static_cast<int>(p);
        m = MAX(1, MIN(m, pairEAM->nrho - 1));
        p -= m;
        p = MIN(p, 1.0);
        coeff = pairEAM->frho_spline[pairEAM->type2frho[jtype]][m];
        double oldF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        // Calculate new embedding energy of atom j.
        p = (rho[j] + delta_rho) * pairEAM->rdrho + 1.0;
        m = static_cast<int>(p);
        m = MAX(1, MIN(m, pairEAM->nrho - 1));
        p -= m;
        p = MIN(p, 1.0);
        coeff = pairEAM->frho_spline[pairEAM->type2frho[jtype]][m];
        double newF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        deltaE += newF - oldF;
    }

    ASSERT(rho[flipAtom] > 0.0);
    ASSERT(new_total_rho_i > 0.0);

    // Things are easier if the rho(r) functional does not depend on the type of both atoms I and J (as for Finnis/Sinclair type potentials).
    if(rho[flipAtom] == new_total_rho_i && pairCDEAM->rhoB[flipAtom] == new_total_rhoB_i) {
        // Calculate local concentration at site i. This did not change.
        double x_i = pairCDEAM->rhoB[flipAtom] / rho[flipAtom];

        // Calculate h(x_i) polynomial function.
        double h_i = pairCDEAM->evalH(x_i);

        // This is the energy change at site i due to the cross term:
        deltaE += 0.5 * h_i * deltaD_i;

        // Compute the change in embedding energy of the swapping atom.
        p = rho[flipAtom] * pairEAM->rdrho + 1.0;
        m = static_cast<int>(p);
        m = MAX(1, MIN(m, pairEAM->nrho - 1));
        p -= m;
        p = MIN(p, 1.0);
        coeff = pairEAM->frho_spline[pairEAM->type2frho[oldSpecies]][m];
        double oldF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = pairEAM->frho_spline[pairEAM->type2frho[newSpecies]][m];
        double newF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        deltaE += newF - oldF;
    }
    else {
        // Calculate the new and old concentration at site i.
        double x_i_old = pairCDEAM->rhoB[flipAtom] / rho[flipAtom];
        double x_i_new = new_total_rhoB_i / new_total_rho_i;

        // Calculate h(x_i) polynomial function.
        double old_h_i = pairCDEAM->evalH(x_i_old);
        double new_h_i = pairCDEAM->evalH(x_i_new);

        // This is the energy change at site i due to the cross term:
        deltaE += 0.5 * (new_h_i * (pairCDEAM->D_values[flipAtom] + deltaD_i) - old_h_i * pairCDEAM->D_values[flipAtom]);

        // Compute the change in embedding energy of the swapping atom.
        p = rho[flipAtom] * pairEAM->rdrho + 1.0;
        m = static_cast<int>(p);
        m = MAX(1, MIN(m, pairEAM->nrho - 1));
        p -= m;
        p = MIN(p, 1.0);
        coeff = pairEAM->frho_spline[pairEAM->type2frho[oldSpecies]][m];
        double oldF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        p = new_total_rho_i * pairEAM->rdrho + 1.0;
        m = static_cast<int>(p);
        m = MAX(1, MIN(m, pairEAM->nrho - 1));
        p -= m;
        p = MIN(p, 1.0);
        coeff = pairEAM->frho_spline[pairEAM->type2frho[newSpecies]][m];
        double newF = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        deltaE += newF - oldF;
    }

    return deltaE;
}

#endif

#if TERSOFF_MC_SUPPORT

/*********************************************************************
 * Calculates the change in energy that swapping the given
 * atom would produce. This routine is for the Tersoff potential.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * flipAtomNL [in]
 *   This specifies the atom to be swapped. It's an index into the neighbor list.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new species of the atom. The atom's type is not changed by this routine. It only computes the induced energy change.
 *
 * Return value:
 *   The expected change in total potential energy.
 *********************************************************************/
double FixSemiGrandCanonicalMC::computeEnergyChangeTersoff(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies)
{
    // This routine is called even when no trial move is being performed during the
    // the current iteration to keep the parallel processors in sync. If no trial
    // move is performed then the energy is calculated twice for the same state of the system.
    if(flipAtom >= 0) {
        // Change system. Perform trial move.
        atom->type[flipAtom] = newSpecies;
    }
    // Transfer changed atom types of the real atoms to the ghost atoms.
    communicationStage = 3;
    comm->forward_comm(this);

    // Calculate new total energy.
    double newEnergy = 0;
    if(flipAtom >= 0) {
        newEnergy = computeEnergyTersoff(flipAtom);
    }

    // Undo trial move. Restore old system state.
    if(flipAtom >= 0) {
        atom->type[flipAtom] = oldSpecies;
    }
    // Transfer changed atom types of the real atoms to the ghost atoms.
    communicationStage = 3;
    comm->forward_comm(this);

    // Calculate old total energy.
    double oldEnergy = 0;
    if(flipAtom >= 0 || totalPotentialEnergy == 0) {
        totalPotentialEnergy = oldEnergy = computeEnergyTersoff(flipAtom);
    }
    else oldEnergy = totalPotentialEnergy;

    return newEnergy - oldEnergy;
}

/// Computes the energy of the atom group around the flipped atom using the Tersoff potential.
double FixSemiGrandCanonicalMC::computeEnergyTersoff(int flipAtom)
{
    double totalEnergy = 0;
    for(int ii = 0; ii < atom->nlocal; ii++)
        totalEnergy += computeAtomicEnergyTersoff(ii);
    return totalEnergy;
}

/// Computes the energy of an atom using the Tersoff potential.
double FixSemiGrandCanonicalMC::computeAtomicEnergyTersoff(int i)
{
    double **x = atom->x;
    int *tag = atom->tag;
    int *type = atom->type;

    int* numneigh = neighborList->numneigh;
    int** firstneigh = neighborList->firstneigh;

    double atomicEnergy = 0;

    int itype = pairTersoff->map[type[i]];
    int itag = tag[i];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];

    // two-body interactions, skip half of them
    int* jlist = firstneigh[i];
    int jnum = numneigh[i];
    for(int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;
        int jtype = pairTersoff->map[type[j]];
        int jtag = tag[j];

        if (itag > jtag) {
            if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
            if ((itag+jtag) % 2 == 1) continue;
        } else {
            if (x[j][2] < x[i][2]) continue;
            if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
            if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }

        double delx = xtmp - x[j][0];
        double dely = ytmp - x[j][1];
        double delz = ztmp - x[j][2];
        double rsq = delx*delx + dely*dely + delz*delz;

        int iparam_ij = pairTersoff->elem2param[itype][jtype][jtype];
        if (rsq > pairTersoff->params[iparam_ij].cutsq) continue;

        double r = sqrt(rsq);
        double tmp_fc = pairTersoff->ters_fc(r, &pairTersoff->params[iparam_ij]);
        double tmp_exp = exp(-pairTersoff->params[iparam_ij].lam1 * r);
        atomicEnergy += tmp_fc * pairTersoff->params[iparam_ij].biga * tmp_exp;
    }

    // three-body interactions
    // skip immediately if I-J is not within cutoff
    for(int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;
        int jtype = pairTersoff->map[type[j]];
        int iparam_ij = pairTersoff->elem2param[itype][jtype][jtype];

        double delr1[3];
        delr1[0] = x[j][0] - xtmp;
        delr1[1] = x[j][1] - ytmp;
        delr1[2] = x[j][2] - ztmp;
        double rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
        if(rsq1 > pairTersoff->params[iparam_ij].cutsq) continue;

        // accumulate bondorder zeta for each i-j interaction via loop over k
        double zeta_ij = 0.0;
        for (int kk = 0; kk < jnum; kk++) {
            if (jj == kk) continue;
            int k = jlist[kk];
            k &= NEIGHMASK;
            int ktype = pairTersoff->map[type[k]];
            int iparam_ijk = pairTersoff->elem2param[itype][jtype][ktype];
            double delr2[3];
            delr2[0] = x[k][0] - xtmp;
            delr2[1] = x[k][1] - ytmp;
            delr2[2] = x[k][2] - ztmp;
            double rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
            if(rsq2 > pairTersoff->params[iparam_ijk].cutsq) continue;
            zeta_ij += pairTersoff->zeta(&pairTersoff->params[iparam_ijk],rsq1,rsq2,delr1,delr2);
        }

        double r = sqrt(rsq1);
        double fa = pairTersoff->ters_fa(r,&pairTersoff->params[iparam_ij]);
        double bij = pairTersoff->ters_bij(zeta_ij,&pairTersoff->params[iparam_ij]);
        atomicEnergy += 0.5*bij*fa;
    }

    return atomicEnergy;
}

#endif

/*********************************************************************
 * Calculates the change in energy that swapping the given atom would produce.
 * This routine is for the general case of an arbitrary potential and
 * IS VERY SLOW! It computes the total energies of the system for the unmodified state
 * and for the modified state and then returns the difference of both values.
 * This routine should only be used for debugging purposes.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new species of the atom. The atom's type is not changed by this method. It only computes the induced energy change.
 *
 * Return value:
 *   The expected change in total potential energy.
 *********************************************************************/
double FixSemiGrandCanonicalMC::computeEnergyChangeGeneric(int flipAtom, int oldSpecies, int newSpecies)
{
    // This routine is called even when no trial move is being performed during the
    // the current iteration to keep the parallel processors in sync. If no trial
    // move is performed then the energy is calculated twice for the same state of the system.
    if(flipAtom >= 0) {
        // Change system. Perform trial move.
        atom->type[flipAtom] = newSpecies;
    }
    // Transfer changed atom types of the real atoms to the ghost atoms.
    communicationStage = 3;
    comm->forward_comm(this);

    // Calculate new total energy.
    double newEnergy = computeTotalEnergy();

    // Undo trial move. Restore old system state.
    if(flipAtom >= 0) {
        atom->type[flipAtom] = oldSpecies;
    }
    // Transfer changed atom types of the real atoms to the ghost atoms.
    communicationStage = 3;
    comm->forward_comm(this);

    // Calculate old total energy.
    double oldEnergy = computeTotalEnergy();

    // Restore the correct electron densities and forces.
    update->integrate->setup_minimal(0);
    fetchGhostAtomElectronDensities();

    return newEnergy - oldEnergy;
}

/*********************************************************************
 * Lets LAMMPS calculate the total potential energy of the system.
 *********************************************************************/
double FixSemiGrandCanonicalMC::computeTotalEnergy()
{
    ASSERT(compute_pe != nullptr);

    int eflag = 1;
    int vflag = 0;

    if(force->pair) force->pair->compute(eflag,vflag);

    if(atom->molecular) {
        if(force->bond) force->bond->compute(eflag,vflag);
        if(force->angle) force->angle->compute(eflag,vflag);
        if(force->dihedral) force->dihedral->compute(eflag,vflag);
        if(force->improper) force->improper->compute(eflag,vflag);
    }

    if(force->kspace) force->kspace->compute(eflag,vflag);

    update->eflag_global = update->ntimestep;
    return compute_pe->compute_scalar();
}

/*********************************************************************
 * Flips the type of one atom and changes the electron densities
 * of nearby atoms accordingly.
 * This routine is for the case of a standard EAM potential.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * flipAtomNL [in]
 *   This specifies the atom to be swapped. It's an index into the neighbor list.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new type to be assigned to the atom.
 *********************************************************************/
void FixSemiGrandCanonicalMC::flipAtomEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies)
{
    double p;
    int m;
    double* rho = pairEAM->rho;
    double* coeff;
    double new_total_rho_i = 0.0;

    // Change atom's type and mark it for exchange.
    atom->type[flipAtom] = newSpecies;
    changedAtoms[flipAtom] = true;

    // Rescale particle velocity vector to conserve kinetic energy.
    double vScaleFactor = sqrt(atom->mass[oldSpecies] / atom->mass[newSpecies]);
    atom->v[flipAtom][0] *= vScaleFactor;
    atom->v[flipAtom][1] *= vScaleFactor;
    atom->v[flipAtom][2] *= vScaleFactor;

    double xi = atom->x[flipAtom][0];
    double yi = atom->x[flipAtom][1];
    double zi = atom->x[flipAtom][2];

    // Loop over all neighbors of the selected atom.
    int* jlist = neighborList->firstneigh[flipAtomNL];
    int jnum = neighborList->numneigh[flipAtomNL];
    for(int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        ASSERT(j < (int)changedAtoms.size());

        double delx = xi - atom->x[j][0];
        double dely = yi - atom->x[j][1];
        double delz = zi - atom->x[j][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        if(rsq >= pairEAM->cutforcesq) continue;

        int jtype = atom->type[j];
        double r = sqrt(rsq);
        p = r * pairEAM->rdr + 1.0;
        m = static_cast<int>(p);
        m = MIN(m, pairEAM->nr - 1);
        p -= m;
        p = MIN(p, 1.0);

        // Calculate change of electron density at site j.
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[oldSpecies][jtype]][m];
        double oldrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[newSpecies][jtype]][m];
        double newrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        double delta_rho = newrho_contr - oldrho_contr;

        rho[j] += delta_rho;

        // Sum total rho at site of swapped atom.
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[jtype][newSpecies]][m];
        new_total_rho_i += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        // Set the flag for this atom to indicate that its rho has changed and needs
        // to be transfered at end of MC step.
        changedAtoms[j] = true;
    }

    // Store newly calculated electron density at swapped atom site.
    rho[flipAtom] = new_total_rho_i;
}

#if CDEAM_MC_SUPPORT

/*********************************************************************
 * Flips the type of one atom and changes the electron densities and
 * D values of nearby atoms accordingly.
 * This routine is for the case of the concentration dependent CD-EAM potential.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * flipAtomNL [in]
 *   This specifies the atom to be swapped. It's an index into the neighbor list.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new type to be assigned to the atom.
 *********************************************************************/
void FixSemiGrandCanonicalMC::flipAtomCDEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies)
{
    double p;
    int m;
    double* rho = pairEAM->rho;
    double* coeff;
    double new_total_rho_i = 0.0;
    double new_total_rhoB_i = 0.0;

    atom->type[flipAtom] = newSpecies;

    // Rescale particle velocity vector to conserve kinetic energy.
    double vScaleFactor = sqrt(atom->mass[oldSpecies] / atom->mass[newSpecies]);
    atom->v[flipAtom][0] *= vScaleFactor;
    atom->v[flipAtom][1] *= vScaleFactor;
    atom->v[flipAtom][2] *= vScaleFactor;

    double xi = atom->x[flipAtom][0];
    double yi = atom->x[flipAtom][1];
    double zi = atom->x[flipAtom][2];

    /// The change in atom i's D value.
    double deltaD_i = 0.0;

    // Loop over all neighbors of the selected atom.
    int* jlist = neighborList->firstneigh[flipAtomNL];
    int jnum = neighborList->numneigh[flipAtomNL];
    for(int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        ASSERT(j < (int)changedAtoms.size());

        double delx = xi - atom->x[j][0];
        double dely = yi - atom->x[j][1];
        double delz = zi - atom->x[j][2];
        double rsq = delx*delx + dely*dely + delz*delz;
        if(rsq >= pairEAM->cutforcesq) continue;

        int jtype = atom->type[j];
        double r = sqrt(rsq);
        p = r * pairEAM->rdr + 1.0;
        m = static_cast<int>(p);
        m = MIN(m, pairEAM->nr - 1);
        p -= m;
        p = MIN(p, 1.0);

        // Calculate change of electron density at site j.
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[oldSpecies][jtype]][m];
        double oldrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[newSpecies][jtype]][m];
        double newrho_contr = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        double delta_rho = newrho_contr - oldrho_contr;

        rho[j] += delta_rho;
        // Sum total rho at site of swapped atom.
        coeff = pairEAM->rhor_spline[pairEAM->type2rhor[jtype][newSpecies]][m];
        double rho_ji = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        new_total_rho_i += rho_ji;
        if(jtype == pairCDEAM->speciesB)
            new_total_rhoB_i += rho_ji;

        // Determine change of partial electron density at site j.
        // It increases if the atom i becomes a B atom and it decreases if it was a B atom.
        if(newSpecies == pairCDEAM->speciesB)
            pairCDEAM->rhoB[j] += newrho_contr;
        else if(oldSpecies == pairCDEAM->speciesB)
            pairCDEAM->rhoB[j] -= oldrho_contr;

        /// The change in atom j's D value.
        double deltaD_j;

        /// Was the old i-j interaction concentration dependent?
        if((jtype == pairCDEAM->speciesA && oldSpecies == pairCDEAM->speciesB)
                || (jtype == pairCDEAM->speciesB && oldSpecies == pairCDEAM->speciesA)) {

            // The old pair interaction was concentration dependent.
            // Please note that it must now become concentration independent that
            // either the A or the B atom has changed to another species.
            // We here require that only the A-B interactions are concentration dependent.

            // Calculate change of pair energy ij.
            coeff = pairEAM->z2r_spline[pairEAM->type2z2r[oldSpecies][jtype]][m];
            double oldz2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

            // The D value of j decreases due to the removed cross interaction.
            deltaD_j = -oldz2 / r;
        }
        else {

            // The old pair interaction was not concentration dependent. Now it might have
            // become dependent. Do check:
            if((jtype == pairCDEAM->speciesA && newSpecies == pairCDEAM->speciesB)
                || (jtype == pairCDEAM->speciesB && newSpecies == pairCDEAM->speciesA)) {

                // The new pair interaction is concentration dependent. It's an AB interaction.

                // Calculate change of pair energy ij.
                coeff = pairEAM->z2r_spline[pairEAM->type2z2r[newSpecies][jtype]][m];
                double newz2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

                // The D value of j increases due to the created cross interaction.
                deltaD_j = newz2 / r;
            }
            else deltaD_j = 0.0;
        }
        pairCDEAM->D_values[j] += deltaD_j;

        // D_i changes in the same way as D_j changes.
        deltaD_i += deltaD_j;

        // Set the flag for this atom to indicate that its rho has changed and needs
        // to be transfered at end of MC step.
        changedAtoms[j] = true;
    }

    // Store newly calculated electron densities and D value at swapped atom site.
    rho[flipAtom] = new_total_rho_i;
    pairCDEAM->rhoB[flipAtom] = new_total_rhoB_i;
    pairCDEAM->D_values[flipAtom] += deltaD_i;

    changedAtoms[flipAtom] = true;
}

#endif

/*********************************************************************
 * Flips the type of one atom.
 * This routine is for the generic case.
 *
 * Parameters:
 *
 * flipAtom [in]
 *   This specifies the atom to be swapped. It's an index into the local list of atoms.
 *
 * oldSpecies [in]
 *   The current species of the atom before the routine is called.
 *
 * newSpecies [in]
 *   The new type to be assigned to the atom.
 *********************************************************************/
void FixSemiGrandCanonicalMC::flipAtomGeneric(int flipAtom, int oldSpecies, int newSpecies)
{
    atom->type[flipAtom] = newSpecies;

    // Rescale particle velocity vector to conserve kinetic energy.
    double vScaleFactor = sqrt(atom->mass[oldSpecies] / atom->mass[newSpecies]);
    atom->v[flipAtom][0] *= vScaleFactor;
    atom->v[flipAtom][1] *= vScaleFactor;
    atom->v[flipAtom][2] *= vScaleFactor;

    changedAtoms[flipAtom] = true;
}

/*********************************************************************
 * Lets the fix report one of its internal state variables to LAMMPS.
 *********************************************************************/
double FixSemiGrandCanonicalMC::compute_vector(int index)
{
    if(index == 0) return nAcceptedSwaps;
    if(index == 1) return nRejectedSwaps;
    index -= 1;
    ASSERT(index >= 1 && index < (int)speciesCounts.size());
    int totalAtoms = 0;
    for(int i = 0; i < (int)speciesCounts.size(); i++)
        totalAtoms += speciesCounts[i];
    if(index  <= atom->ntypes) return (double)speciesCounts[index] / (totalAtoms > 0 ? totalAtoms : 1);
    ASSERT(false); return 0.0;
}

/*********************************************************************
 * Reports the memory usage of this fix to LAMMPS.
 *********************************************************************/
double FixSemiGrandCanonicalMC::memory_usage()
{
    return (changedAtoms.size() * sizeof(bool)) +
        (samplingWindowAtoms.size() * sizeof(int));
}

