/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef PSI4_CORE_GIAO_H
#define PSI4_CORE_GIAO_H

#include "onebody.h"
#include "osrecur.h"

namespace psi {

/**
\begin{equation}
  S^{B_i}_{\mu\nu}=\left[\left({\rm M}_j-{\rm N}_j\right)
  \left\langle\mu \left| {\rm k}\right|\nu\right\rangle-\left({\rm M}_k-{\rm N}_k\right)
  \left\langle\mu \left| {\rm j}\right|\nu\right\rangle\right]
\end{equation}
*/

struct GIAO_OverlapInt : public OneBodyAOInt
{
private:
    //! Obara an Saika recursion object to be used.
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    //! Computes the GIAO overlap between two gaussian shells
    void compute_pair(const GaussianShell&, const GaussianShell&);

public:

    //! Constructor. Use an IntegralFactory
    GIAO_OverlapInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);

    //! Virtual destructor
    virtual ~GIAO_OverlapInt();
};

} // namespace psi

#endif //PSI4_CORE_GIAO_H
