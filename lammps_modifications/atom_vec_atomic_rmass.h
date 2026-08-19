#ifdef ATOM_CLASS
// clang-format off
AtomStyle(atomic/rmass,AtomVecAtomicRmass);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_ATOMIC_RMASS_H
#define LMP_ATOM_VEC_ATOMIC_RMASS_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecAtomicRmass : virtual public AtomVec {
 public:
  AtomVecAtomicRmass(class LAMMPS *);
  void grow_pointers() override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;
  void pack_data_pre(int) override;
  void pack_data_post(int) override;

 protected:
  double *rmass;
  double rmass_one;
};

}    // namespace LAMMPS_NS

#endif
#endif 