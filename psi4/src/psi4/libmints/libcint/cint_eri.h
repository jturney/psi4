/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef PSI4_CORE_CINT_ERI_H
#define PSI4_CORE_CINT_ERI_H

#include "psi4/libmints/twobody.h"

#include <array>
#include <vector>
#include <memory>

namespace psi {

class BasisSet;

class CINTERI : public TwoBodyAOInt
{
    std::vector<int> atm_;
    std::vector<int> bas_;
    std::vector<double> env_;

    std::array<int, 4> basis_shell_start_;
    int nshell_;
    int natom_;

public:
    CINTERI(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
            std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4,
            int deriv = 0);

    virtual ~CINTERI()
    { }

    virtual size_t compute_shell(int, int, int, int) override;

    virtual size_t compute_shell(const AOShellCombinationsIterator&) override;

    virtual size_t compute_shell_deriv1(int, int, int, int) override;

    virtual size_t compute_shell_deriv2(int, int, int, int) override;
};

} // psi namespace

#endif //PSI4_CORE_CINT_ERI_H
