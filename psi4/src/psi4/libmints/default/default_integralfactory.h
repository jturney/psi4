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

#ifndef PSI4_CORE_DEFAULT_INTEGRALFACTORY_HH
#define PSI4_CORE_DEFAULT_INTEGRALFACTORY_HH

#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include <vector>

#include "psi4/libmints/onebody.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"

namespace psi {

class DefaultIntegralFactory
{
protected:

    /// Center 1 basis set
    std::shared_ptr<BasisSet> bs1_;
    /// Center 2 basis set
    std::shared_ptr<BasisSet> bs2_;
    /// Center 3 basis set
    std::shared_ptr<BasisSet> bs3_;
    /// Center 4 basis set
    std::shared_ptr<BasisSet> bs4_;

    /// Provides ability to transform to sphericals (d=0, f=1, g=2)
    std::vector<SphericalTransform> spherical_transforms_;
    /// Provides ability to transform from sphericals (d=0, f=1, g=2)
    std::vector<ISphericalTransform> ispherical_transforms_;

    /// Initializes spherical harmonic transformations
    void init_spherical_harmonics(int max_am);

public:

    /** Initialize IntegralFactory object given a BasisSet for each center. */
    DefaultIntegralFactory(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
            std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4);
    /** Initialize IntegralFactory object given a BasisSet for two centers. Becomes (bs1 bs1 | bs1 bs1). */
    explicit DefaultIntegralFactory(std::shared_ptr<BasisSet> bs1);

    /// Return the basis set on center 1.
    std::shared_ptr<BasisSet> basis1() const;
    /// Return the basis set on center 2.
    std::shared_ptr<BasisSet> basis2() const;
    /// Return the basis set on center 3.
    std::shared_ptr<BasisSet> basis3() const;
    /// Return the basis set on center 4.
    std::shared_ptr<BasisSet> basis4() const;

    /// Set the basis set for each center.
    virtual void set_basis(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                           std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4);

    /// Returns an OneBodyInt that computes the overlap integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_overlap(int deriv=0);

    /// Returns an OneBodyInt that computes the overlap integral.
    virtual std::unique_ptr<OneBodySOInt> so_overlap(int deriv=0);

    /// Returns a ThreeCenterOverlapINt that computes the overlap between three centers
    virtual std::unique_ptr<ThreeCenterOverlapInt> overlap_3c();

    /// Returns an OneBodyInt that computes the kinetic energy integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_kinetic(int deriv=0);
    virtual std::unique_ptr<OneBodySOInt> so_kinetic(int deriv=0);

    /// Returns an OneBodyInt that computes the nuclear attraction integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_potential(int deriv=0);
    virtual std::unique_ptr<OneBodySOInt> so_potential(int deriv=0);

    /// Returns an OneBodyInt that computes the ECP integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_ecp(int deriv=0);
    virtual std::unique_ptr<OneBodySOInt> so_ecp(int deriv=0);

    /// Returns an OneBodyInt that computes the relativistic nuclear attraction integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_rel_potential(int deriv=0);
    virtual std::unique_ptr<OneBodySOInt> so_rel_potential(int deriv=0);

    /// Returns the OneBodyInt that computes the pseudospectral grid integrals
    virtual std::unique_ptr<OneBodyAOInt> ao_pseudospectral(int deriv = 0);
    virtual std::unique_ptr<OneBodySOInt> so_pseudospectral(int deriv = 0);

    /// Returns an OneBodyInt that computes the dipole integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_dipole(int deriv=0);
    virtual std::unique_ptr<OneBodySOInt> so_dipole(int deriv=0);

    /// Returns an OneBodyInt that computes the quadrupole integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_quadrupole();
    virtual std::unique_ptr<OneBodySOInt> so_quadrupole();

    /// Returns an OneBodyInt that computes arbitrary-order multipole integrals.
    virtual std::unique_ptr<OneBodyAOInt> ao_multipoles(int order);
    virtual std::unique_ptr<OneBodySOInt> so_multipoles(int order);

    /// Returns an OneBodyInt that computes the traceless quadrupole integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_traceless_quadrupole();
    virtual std::unique_ptr<OneBodySOInt> so_traceless_quadrupole();

    /// Returns an OneBodyInt that computes the nabla integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_nabla(int deriv=0);
    virtual std::unique_ptr<OneBodySOInt> so_nabla(int deriv=0);

    /// Returns an OneBodyInt that computes the nabla integral.
    virtual std::unique_ptr<OneBodyAOInt> ao_angular_momentum(int deriv=0);
    virtual std::unique_ptr<OneBodySOInt> so_angular_momentum(int deriv=0);

    /// Returns a OneBodyInt that computes the multipole potential integrals for EFP
    virtual std::unique_ptr<OneBodyAOInt> ao_efp_multipole_potential(int deriv=0);
    virtual std::unique_ptr<OneBodySOInt> so_efp_multipole_potential(int deriv=0);

    /// Returns an OneBodyInt that computes the electric field
    virtual std::unique_ptr<OneBodyAOInt> electric_field();

    /// Returns an OneBodyInt that computes the point electrostatic potential
    virtual std::unique_ptr<OneBodyAOInt> electrostatic();

    /// Returns an OneBodyInt that computes the electrostatic potential at desired points
    /// Want to change the name of this after the PCM dust settles
    virtual std::unique_ptr<OneBodyAOInt> pcm_potentialint();

    /// Returns an ERI integral object
    virtual std::unique_ptr<TwoBodyAOInt> eri(int deriv=0, bool use_shell_pairs=true);

    /// Returns an ERD ERI integral object, if available.  Otherwise returns a libint integral object
    virtual std::unique_ptr<TwoBodyAOInt> erd_eri(int deriv=0, bool use_shell_pairs=true);

    /// Returns an erf ERI integral object (omega integral)
    virtual std::unique_ptr<TwoBodyAOInt> erf_eri(double omega, int deriv=0, bool use_shell_pairs=true);

    /// Returns an erf complement ERI integral object (omega integral)
    virtual std::unique_ptr<TwoBodyAOInt> erf_complement_eri(double omega, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12 integral object
    virtual std::unique_ptr<TwoBodyAOInt> f12(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12Scaled integral object
    virtual std::unique_ptr<TwoBodyAOInt> f12_scaled(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12 squared integral object
    virtual std::unique_ptr<TwoBodyAOInt> f12_squared(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12G12 integral object
    virtual std::unique_ptr<TwoBodyAOInt> f12g12(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

    /// Returns an F12 double commutator integral object
    virtual std::unique_ptr<TwoBodyAOInt> f12_double_commutator(std::shared_ptr<CorrelationFactor> cf, int deriv=0, bool use_shell_pairs=true);

};

}

#endif //PSI4_CORE_DEFAULT_INTEGRALFACTORY_HH
