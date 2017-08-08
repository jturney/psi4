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

#include "default_integralfactory.h"

#include "psi4/libmints/default/rel_potential.h"
#include "psi4/libmints/default/electricfield.h"
#include "psi4/libmints/default/tracelessquadrupole.h"
#include "psi4/libmints/default/efpmultipolepotential.h"
#include "psi4/libmints/default/eri.h"
#include "psi4/libmints/default/multipoles.h"
#include "psi4/libmints/default/quadrupole.h"
#include "psi4/libmints/default/angularmomentum.h"
#include "psi4/libmints/default/nabla.h"
#include "psi4/libmints/default/dipole.h"
#include "psi4/libmints/default/electrostatic.h"
#include "pseudospectral.h"
#include "psi4/libmints/default/kinetic.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libmints/default/overlap.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/default/potentialint.h"
#include "psi4/libmints/ecpint.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/liberd/erd_eri.h"

#ifdef USING_simint
#include "psi4/libmints/simint/siminteri.h"
#endif

#include <libint/libint.h>

namespace psi {

namespace {
template<typename BaseT, typename T, typename... Ts>
std::unique_ptr<BaseT> make_unique(Ts&&... params)
{
    return std::unique_ptr<BaseT>(new T(std::forward<Ts>(params)...));
}
}

DefaultIntegralFactory::DefaultIntegralFactory(std::shared_ptr<BasisSet> bs1,
                                               std::shared_ptr<BasisSet> bs2,
                                               std::shared_ptr<BasisSet> bs3,
                                               std::shared_ptr<BasisSet> bs4)
{

    set_basis(bs1, bs2, bs3, bs4);
}

DefaultIntegralFactory::DefaultIntegralFactory(std::shared_ptr<BasisSet> bs1)
{
    set_basis(bs1, bs1, bs1, bs1);
}

std::shared_ptr<BasisSet> DefaultIntegralFactory::basis1() const
{
    return bs1_;
}

std::shared_ptr<BasisSet> DefaultIntegralFactory::basis2() const
{
    return bs2_;
}

std::shared_ptr<BasisSet> DefaultIntegralFactory::basis3() const
{
    return bs3_;
}

std::shared_ptr<BasisSet> DefaultIntegralFactory::basis4() const
{
    return bs4_;
}

void DefaultIntegralFactory::set_basis(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                       std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4)
{
    bs1_ = bs1;
    bs2_ = bs2;
    bs3_ = bs3;
    bs4_ = bs4;

    // Use the max am from libint
    init_spherical_harmonics(LIBINT_MAX_AM + 1);
}

void DefaultIntegralFactory::init_spherical_harmonics(int max_am)
{
    spherical_transforms_.clear();
    ispherical_transforms_.clear();

    for (int i = 0; i <= max_am; ++i) {
        spherical_transforms_.emplace_back(SphericalTransform(i));
        ispherical_transforms_.emplace_back(ISphericalTransform(i));
    }
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_overlap(int deriv)
{
    return make_unique<OneBodyAOInt, OverlapInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_overlap(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_overlap(deriv));
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<ThreeCenterOverlapInt> DefaultIntegralFactory::overlap_3c()
{
    return make_unique<ThreeCenterOverlapInt, ThreeCenterOverlapInt>(bs1_, bs2_, bs3_);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_kinetic(int deriv)
{
    return make_unique<OneBodyAOInt, KineticInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_kinetic(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_kinetic(deriv));
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_potential(int deriv)
{
    return make_unique<OneBodyAOInt, PotentialInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_potential(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_potential(deriv));
    return make_unique<OneBodySOInt, PotentialSOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_ecp(int deriv)
{
    return make_unique<OneBodyAOInt, ECPInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_ecp(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_ecp(deriv));
    return make_unique<OneBodySOInt, ECPSOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_rel_potential(int deriv)
{
    return make_unique<OneBodyAOInt, RelPotentialInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_rel_potential(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_rel_potential(deriv));
    return make_unique<OneBodySOInt, RelPotentialSOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_pseudospectral(int deriv)
{
    return make_unique<OneBodyAOInt, PseudospectralInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_pseudospectral(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_pseudospectral(deriv));
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::electrostatic()
{
    return make_unique<OneBodyAOInt, ElectrostaticInt>(bs1_, bs2_, 0);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::pcm_potentialint()
{
    return make_unique<OneBodyAOInt, PCMPotentialInt>(bs1_, bs2_, 0);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_dipole(int deriv)
{
    return make_unique<OneBodyAOInt, DipoleInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_dipole(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_dipole(deriv));
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_nabla(int deriv)
{
    return make_unique<OneBodyAOInt, NablaInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_nabla(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_nabla(deriv));
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_angular_momentum(int deriv)
{
    return make_unique<OneBodyAOInt, AngularMomentumInt>(bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_angular_momentum(int deriv)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_angular_momentum(deriv));
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_quadrupole()
{
    return make_unique<OneBodyAOInt, QuadrupoleInt>(bs1_, bs2_);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_quadrupole()
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_quadrupole());
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_multipoles(int order)
{
    return make_unique<OneBodyAOInt, MultipoleInt>(bs1_, bs2_, order);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_efp_multipole_potential(int order)
{
    return make_unique<OneBodyAOInt, EFPMultipolePotentialInt>(bs1_, bs2_, order);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_efp_multipole_potential(int order)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_efp_multipole_potential(order));
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_multipoles(int order)
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_multipoles(order));
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::ao_traceless_quadrupole()
{
    return make_unique<OneBodyAOInt, TracelessQuadrupoleInt>(bs1_, bs2_);
}

std::unique_ptr<OneBodySOInt> DefaultIntegralFactory::so_traceless_quadrupole()
{
    std::shared_ptr<OneBodyAOInt> ao_int(ao_traceless_quadrupole());
    return make_unique<OneBodySOInt, OneBodySOInt>(ao_int);
}

std::unique_ptr<OneBodyAOInt> DefaultIntegralFactory::electric_field()
{
    return make_unique<OneBodyAOInt, ElectricFieldInt>(bs1_, bs2_);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::erd_eri(int deriv, bool use_shell_pairs)
{
#ifdef USING_simint
    if (deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "SIMINT")
        return make_unique<TwoBodyAOInt, SimintERI>(bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
#elif defined USING_erd
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "ERD")
        return new ERDERI(deriv, use_shell_pairs);
#endif
    return eri(deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::eri(int deriv, bool use_shell_pairs)
{
#ifdef USING_simint
    if (deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "SIMINT")
        return make_unique<TwoBodyAOInt, SimintERI>(bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
#elif defined USING_erd
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "ERD")
        return make_unique<TwoBodyAOInt, ERDERI>(deriv, use_shell_pairs);
#endif
    return make_unique<TwoBodyAOInt, ERI>(bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::erf_eri(double omega, int deriv, bool use_shell_pairs)
{
    return make_unique<TwoBodyAOInt, ErfERI>(omega, bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::erf_complement_eri(double omega, int deriv, bool use_shell_pairs)
{
    return make_unique<TwoBodyAOInt, ErfComplementERI>(omega, bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::f12(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<TwoBodyAOInt, F12>(cf, bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::f12_scaled(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<TwoBodyAOInt, F12Scaled>(cf, bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::f12_squared(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<TwoBodyAOInt, F12Squared>(cf, bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::f12g12(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<TwoBodyAOInt, F12G12>(cf, bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> DefaultIntegralFactory::f12_double_commutator(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<TwoBodyAOInt, F12DoubleCommutator>(cf, bs1_, bs2_, bs3_, bs4_, deriv, use_shell_pairs);
}

}
