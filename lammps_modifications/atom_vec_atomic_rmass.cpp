/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_atomic_rmass.h"
#include "atom.h"
#include "error.h"
#include "fix.h"
#include "math_const.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecAtomicRmass::AtomVecAtomicRmass(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_ATOM;

  // enable rmass flag in Atom class
  atom->rmass_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"rmass"};
  fields_copy = {"rmass"};
  fields_border = {"rmass"};
  fields_border_vel = {"rmass"};
  fields_exchange = {"rmass"};
  fields_restart = {"rmass"};
  fields_create = {"rmass"};
  fields_data_atom = {"id", "type", "rmass", "x"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecAtomicRmass::grow_pointers()
{
  rmass = atom->rmass;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecAtomicRmass::create_atom_post(int ilocal)
{
  rmass[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecAtomicRmass::data_atom_post(int ilocal)
{
  if (rmass[ilocal] <= 0.0) error->one(FLERR, "Invalid mass in Atoms section of data file");
} 

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecAtomicRmass::pack_data_pre(int ilocal)
{
  rmass_one = rmass[ilocal];
}


/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecAtomicRmass::pack_data_post(int ilocal)
{
  rmass[ilocal] = rmass_one;
}
