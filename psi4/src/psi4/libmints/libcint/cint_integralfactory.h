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

#ifndef PSI4_CORE_CINT_INTEGRALFACTORY_H
#define PSI4_CORE_CINT_INTEGRALFACTORY_H

#include "psi4/libmints/default/default_integralfactory.h"

namespace psi {

class CINTIntegralFactory : public DefaultIntegralFactory
{
public:
    /** Initialize IntegralFactory object given a BasisSet for each center. */
    CINTIntegralFactory(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
            std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4);
    /** Initialize IntegralFactory object given a BasisSet for two centers. Becomes (bs1 bs1 | bs1 bs1). */
    explicit CINTIntegralFactory(std::shared_ptr<BasisSet> bs1);

    /// Returns an ERI integral object
    std::unique_ptr<TwoBodyAOInt> eri(int deriv=0, bool use_shell_pairs=true) override;
};

} // psi namespace
#endif //PSI4_CORE_CINT_INTEGRALFACTORY_H
