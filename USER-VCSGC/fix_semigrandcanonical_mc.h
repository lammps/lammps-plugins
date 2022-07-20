/* ----------------------------------------------------------------------
 * Parallel Monte-Carlo code for the semi-grandcanonical ensemble (SGC)
 * and the variance-constrained semi-grandcanonical ensemble (VC-SGC).
 *
 * See Sadigh et al., Phys. Rev. B 85, 184203 (2012) for a 
 * description of the algorithm.
 *
 * Code author: Alexander Stukowski (stukowski@mm.tu-darmstadt.de)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sgcmc,FixSemiGrandCanonicalMC)

#else

#ifndef FIX_SEMIGRANDCANONICAL_MC_H
#define FIX_SEMIGRANDCANONICAL_MC_H

#include "fix.h"

// Setting this to 1 enables support for the concentration-dependent EAM potential (pair_style eam/cd) 
// in the Monte Carlo routine. Setting to 0 limits support to standard EAM only and removes all dependencies
// on the CD-EAM potential code.
#ifndef CDEAM_MC_SUPPORT
#define CDEAM_MC_SUPPORT                0
#endif

// Setting this to 1 enables support for Tersoff-like potentials (pair_style tersoff)
// in the Monte Carlo routine.
#ifndef TERSOFF_MC_SUPPORT
#define TERSOFF_MC_SUPPORT              0
#endif
// Setting this to 1 enables additional debugging/sanity checks (with a small performance penalty).
#ifndef SGCMC_DEBUG
#define SGCMC_DEBUG                             0
#endif

#include <vector>
#include <math.h>

namespace LAMMPS_NS {

class FixSemiGrandCanonicalMC : public Fix
{
public:

	/// Fix class constructor.
	FixSemiGrandCanonicalMC(class LAMMPS*, int, char **);

	/// Fix class destructor.
	virtual ~FixSemiGrandCanonicalMC();

	/******************** Virtual methods from Fix base class ************************/

	/// The return value of this method specifies at which points the fix is invoked during the simulation.
	virtual int setmask();

	/// This gets called by the system before the simulation starts.
	virtual void init();

	/// Assigns the requested neighbor list to the fix.
	virtual void init_list(int id, NeighList *ptr);

	/// Called after the EAM force calculation during each timestep.
	/// This method triggers the MC routine from time to time.
	virtual void post_force(int vflag);

	/// Lets the fix report one of its internal state variables to LAMMPS.
	virtual double compute_vector(int index);

	/// This is for MPI communication with neighbor nodes.
	virtual int pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc);
	virtual void unpack_forward_comm(int n, int first, double *buf);
	virtual int pack_reverse_comm(int n, int first, double* buf);
	virtual void unpack_reverse_comm(int n, int *list, double* buf);

	/// Reports the memory usage of this fix to LAMMPS.
	virtual double memory_usage();

	/******************** Monte-Carlo routines ************************/

	/// This routine does one full MC step.
	void doMC();

	/// Fetches the electron densities for the local ghost atoms from the neighbor nodes.
	void fetchGhostAtomElectronDensities();

	/// Positions the sampling window inside the node's bounding box.
	bool placeSamplingWindow();

	/// Calculates the change in energy that swapping the given atom would produce.
	/// This routine is for the case of a standard EAM potential.
	double computeEnergyChangeEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies);

#if CDEAM_MC_SUPPORT
	/// Calculates the change in energy that swapping the given atom would produce.
	/// This routine is for the case of the concentration dependent CD-EAM potential.
	double computeEnergyChangeCDEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies);
#endif

#if TERSOFF_MC_SUPPORT
	/// Calculates the change in energy that swapping the given atom would produce.
	/// This routine is for the Tersoff potential.
	double computeEnergyChangeTersoff(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies);

	/// Computes the energy of the atom group around the flipped atom using the Tersoff potential.
	double computeEnergyTersoff(int flipAtom);

	/// Computes the energy of an atom using the Tersoff potential.
	double computeAtomicEnergyTersoff(int i);
#endif

	/// Calculates the change in energy that swapping the given atom would produce.
	/// This routine is for the general case of an arbitrary potential and
	/// IS VERY SLOW! It computes the total energies of the system for the unmodified state
	/// and for the modified state and then returns the difference of both values.
	/// This routine should only be used for debugging purposes.
	double computeEnergyChangeGeneric(int flipAtom, int oldSpecies, int newSpecies);

	/// Lets LAMMPS calculate the total potential energy of the system.
	double computeTotalEnergy();

	/// Flips the type of one atom and changes the electron densities of nearby atoms accordingly.
	/// This routine is for the case of a standard EAM potential.
	void flipAtomEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies);

	/// Flips the type of one atom and changes the electron densities and D values of nearby atoms accordingly.
	/// This routine is for the case of the concentration dependent CD-EAM potential.
	void flipAtomCDEAM(int flipAtom, int flipAtomNL, int oldSpecies, int newSpecies);

	/// Flips the type of one atom.
	/// This routine is for the generic case.
	void flipAtomGeneric(int flipAtom, int oldSpecies, int newSpecies);

	/// Transfers the locally changed electron densities and atom types to the neighbors.
	void communicateRhoAndTypes();

#if SGCMC_DEBUG
	/// Allocate atom-based array.
	void grow_arrays(int);
	/// Copy values within local atom-based array.
	void copy_arrays(int, int);
	/// Initialize one atom's array values, called when atom is created.
	void set_arrays(int);
	/// Pack values in local atom-based array for exchange with another proc.
	int pack_exchange(int, double *);
	/// Unpack values in local atom-based array from exchange with another proc.
	int unpack_exchange(int, double *);
#endif

private:

	/// Sends the given formatted string to the log file and stdout.
	void printLog(const char* format, ...);

private:

	/// The number of MD steps between each MC cycle.
	int nevery_mdsteps;

	/// The number of times the sampling window should be repositioned during
	/// one MC cycle.
	int numSamplingWindowMoves;

	/// The fraction of atoms that should be swapped per MC step.
	double swap_fraction;

	/// The maximum interaction radius of all potentials.
	double interactionRadius;

	/// The inverse MC temperature.
	double beta;

	/// Chemical potential differences for all species. The differences are relative to the chemical
	/// potential of the first species. Note that this array is based on index 1 (not 0 as normal C arrays).
	/// This means the first two elements of this vector are always zero.
	std::vector<double> deltamu;

	/// Enables serial implementation without second rejection in the VCSGC ensemble.
	bool serialMode;

	/// The MC variance constraint parameter.
	double kappa;

	/// The target concentration values for each species. The concentration of first species is
	/// implicitely defined as one minues all other concentrations. Please note that this vector
	/// is based on index 1. The first element at index 0 is not used.
	std::vector<double> targetConcentration;

	/// The master seed value for the random number generators on all nodes.
	int seed;

        /// The random number generator that is in sync with all other nodes.
        class RanPark* random;

        /// The local random number generator for this proc only.
        class RanPark* localRandom;

	/// The total number of atoms of the different species in the whole system.
	/// Divide this by the total number of atoms to get the global concentration.
	/// Since LAMMPS atom types start at index 1 this array is also based on index 1.
	/// The first array element at index 0 is not used.
	std::vector<int> speciesCounts;

	/// The full neighbor list used by this fix.
	class NeighList* neighborList;

	/// The user-defined size of the sampling window. It is specified as a fraction
	/// of the processor's cell boundaries.
	/// If this parameter is 0 then the default sampling window size is used.
	double samplingWindowUserSize;

	/// This counter is increased each time the sampling window is repositioned.
	/// The lowest 3 bits of this integer value specify the alignment of the sampling window in
	/// the processor cell. That means that after 8 iterations the 8 corners have been sampled
	/// and it starts at the first corner again.
	int samplingWindowPosition;

	/// Array with indices of all atoms that are within the sampling window.
	/// Note 1: that an atom can be more than once in this array if its multiplicity is greater than one.
	/// Note 2: Indices are into the I array of the neighbor list and not into the atoms array.
	std::vector<int> samplingWindowAtoms;

	/// The number of atoms inside the sampling window. Counting each atom only once.
	int numSamplingWindowAtoms;

	/// The number of local atoms that are in the fix group.
	int numFixAtomsLocal;

        /// Pointer to the EAM potential class.
        /// This is required to access the Rho arrays calculated by the potential class and its potential tables.
        class PairEAM* pairEAM;

#if CDEAM_MC_SUPPORT
        /// Pointer to the CD-EAM potential class.
        /// This is required to access the RhoB arrays calculated by the potential class.
        /// The pointer is NULL if only the standard EAM model is used in the simulation.
        class PairEAMCD* pairCDEAM;
#endif

#if TERSOFF_MC_SUPPORT
        /// Pointer to the Tersoff potential class.
        /// This is required to access the parameters of the potential when computing the
        /// change in energy.
        class PairTersoff* pairTersoff;
#endif

	/// This array contains a boolean value per atom (real and ghosts) that indicates whether
	/// the electron density or another property at that site has been affected by one of the accepted MC swaps.
	std::vector<bool> changedAtoms;

	/// This counter indicates the current MPI communication stage to let the
	/// pack/unpack routines know which data is being transmitted.
	int communicationStage;

	/// The total number of accepted swaps during the last MC step.
	int nAcceptedSwaps;

	/// The total number of rejected swaps during the last MC step.
	int nRejectedSwaps;

	/// Keeps track of the current total potential energy.
	/// This is only used when no routine is available that can efficiently calculate the
	/// local energy change due to an atom swap.
	double totalPotentialEnergy;

	/// A compute used to compute the total potential energy of the system.
	class Compute* compute_pe;

#if SGCMC_DEBUG
	/// This per-atom array counts how often each atom is picked for a trial move.
	/// This is only used for debugging purposes.
	double* trialCounters;
#endif
};

};


#endif	// FIX_SEMIGRANDCANONICAL_MC_H

#endif	// FIX_CLASS
