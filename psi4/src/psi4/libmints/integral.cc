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
#include "psi4/libmints/integral.h"
#include "psi4/libmints/shellrotation.h"
#include "psi4/libmints/cartesianiter.h"
#include "psi4/libmints/rel_potential.h"
#include "psi4/libmints/electricfield.h"
#include "psi4/libmints/tracelessquadrupole.h"
#include "psi4/libmints/efpmultipolepotential.h"
#include "psi4/libmints/eri.h"
#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/quadrupole.h"
#include "psi4/libmints/angularmomentum.h"
#include "psi4/libmints/nabla.h"
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/pseudospectral.h"
#include "psi4/libmints/kinetic.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libmints/overlap.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/potentialint.h"
#include "psi4/libmints/ecpint.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/erd_eri.h"

#ifdef USING_simint
#include "psi4/libmints/siminteri.h"
#endif

#include <libint/libint.h>

namespace psi {

namespace {
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}
} // namespace

IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2,
                                 std::shared_ptr<BasisSet> bs3,
                                 std::shared_ptr<BasisSet> bs4)
{

    set_basis(bs1, bs2, bs3, bs4);
}

IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1)
{
    set_basis(bs1, bs1, bs1, bs1);
}

std::shared_ptr<BasisSet> IntegralFactory::basis1() const
{
    return bs1_;
}

std::shared_ptr<BasisSet> IntegralFactory::basis2() const
{
    return bs2_;
}

std::shared_ptr<BasisSet> IntegralFactory::basis3() const
{
    return bs3_;
}

std::shared_ptr<BasisSet> IntegralFactory::basis4() const
{
    return bs4_;
}

void IntegralFactory::set_basis(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4)
{
    bs1_ = bs1;
    bs2_ = bs2;
    bs3_ = bs3;
    bs4_ = bs4;

    // Use the max am from libint
    init_spherical_harmonics(LIBINT_MAX_AM + 1);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_overlap(int deriv)
{
    return make_unique<OverlapInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_overlap(int deriv)
{
    return make_unique<OneBodySOInt>(ao_overlap(deriv), this);
}

std::unique_ptr<ThreeCenterOverlapInt> IntegralFactory::overlap_3c()
{
    return make_unique<ThreeCenterOverlapInt>(spherical_transforms_, bs1_, bs2_, bs3_);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_kinetic(int deriv)
{
    return make_unique<KineticInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_kinetic(int deriv)
{
    return make_unique<OneBodySOInt>(ao_kinetic(deriv), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_potential(int deriv)
{
    return make_unique<PotentialInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_potential(int deriv)
{
    return make_unique<PotentialSOInt>(ao_potential(deriv), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_ecp(int deriv)
{
    return make_unique<ECPInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_ecp(int deriv)
{
    return make_unique<ECPSOInt>(ao_ecp(deriv), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_rel_potential(int deriv)
{
    return make_unique<RelPotentialInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_rel_potential(int deriv)
{
    return make_unique<RelPotentialSOInt>(ao_rel_potential(deriv), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_pseudospectral(int deriv)
{
    return make_unique<PseudospectralInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_pseudospectral(int deriv)
{
    return make_unique<OneBodySOInt>(ao_pseudospectral(deriv), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::electrostatic()
{
    return make_unique<ElectrostaticInt>(spherical_transforms_, bs1_, bs2_, 0);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::pcm_potentialint()
{
    return make_unique<PCMPotentialInt>(spherical_transforms_, bs1_, bs2_, 0);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_dipole(int deriv)
{
    return make_unique<DipoleInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_dipole(int deriv)
{
    return make_unique<OneBodySOInt>(ao_dipole(deriv), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_nabla(int deriv)
{
    return make_unique<NablaInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_nabla(int deriv)
{
    return make_unique<OneBodySOInt>(ao_nabla(deriv), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_angular_momentum(int deriv)
{
    return make_unique<AngularMomentumInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_angular_momentum(int deriv)
{
    return make_unique<OneBodySOInt>(ao_angular_momentum(deriv), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_quadrupole()
{
    return make_unique<QuadrupoleInt>(spherical_transforms_, bs1_, bs2_);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_quadrupole()
{
    return make_unique<OneBodySOInt>(ao_quadrupole(), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_multipoles(int order)
{
    return make_unique<MultipoleInt>(spherical_transforms_, bs1_, bs2_, order);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_efp_multipole_potential(int order)
{
    return make_unique<EFPMultipolePotentialInt>(spherical_transforms_, bs1_, bs2_, order);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_efp_multipole_potential(int order)
{
    return make_unique<OneBodySOInt>(ao_efp_multipole_potential(order), this);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_multipoles(int order)
{
    return make_unique<OneBodySOInt>(ao_multipoles(order), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_traceless_quadrupole()
{
    return make_unique<TracelessQuadrupoleInt>(spherical_transforms_, bs1_, bs2_);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_traceless_quadrupole()
{
    return make_unique<OneBodySOInt>(ao_traceless_quadrupole(), this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::electric_field(int deriv)
{
    return make_unique<ElectricFieldInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erd_eri(int deriv, bool use_shell_pairs)
{
#ifdef USING_simint
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "SIMINT")
        return make_unique<SimintERI>(this, deriv, use_shell_pairs);
#elif defined USING_erd
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "ERD")
        return make_unique<ERDERI>(this, deriv, use_shell_pairs);
#endif
    return eri(deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::eri(int deriv, bool use_shell_pairs)
{
#ifdef USING_simint
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "SIMINT")
        return make_unique<SimintERI>(this, deriv, use_shell_pairs);
#elif defined USING_erd
    if(deriv == 0 && Process::environment.options.get_str("INTEGRAL_PACKAGE") == "ERD")
        return make_unique<ERDERI>(this, deriv, use_shell_pairs);
#endif
    return make_unique<ERI>(this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_eri(double omega, int deriv, bool use_shell_pairs)
{
    return make_unique<ErfERI>(omega, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_complement_eri(double omega, int deriv, bool use_shell_pairs)
{
    return make_unique<ErfComplementERI>(omega, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<F12>(cf, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12_scaled(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<F12Scaled>(cf, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12_squared(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<F12Squared>(cf, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12g12(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<F12G12>(cf, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12_double_commutator(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return make_unique<F12DoubleCommutator>(cf, this, deriv, use_shell_pairs);
}

void IntegralFactory::init_spherical_harmonics(int max_am)
{
    spherical_transforms_.clear();
    ispherical_transforms_.clear();

    for (int i = 0; i <= max_am; ++i) {
        spherical_transforms_.push_back(SphericalTransform(i));
        ispherical_transforms_.push_back(ISphericalTransform(i));
    }
}

AOShellCombinationsIterator IntegralFactory::shells_iterator()
{
    return AOShellCombinationsIterator(bs1_, bs2_, bs3_, bs4_);
}

AOShellCombinationsIterator *IntegralFactory::shells_iterator_ptr()
{
    return new AOShellCombinationsIterator(bs1_, bs2_, bs3_, bs4_);
}

AOIntegralsIterator IntegralFactory::integrals_iterator(int p, int q, int r, int s)
{
    return AOIntegralsIterator(bs1_->shell(p), bs2_->shell(q), bs3_->shell(r), bs4_->shell(s));
}

CartesianIter *IntegralFactory::cartesian_iter(int l) const
{
    return new CartesianIter(l);
}

RedundantCartesianIter *IntegralFactory::redundant_cartesian_iter(int l) const
{
    return new RedundantCartesianIter(l);
}

RedundantCartesianSubIter *IntegralFactory::redundant_cartesian_sub_iter(int l) const
{
    return new RedundantCartesianSubIter(l);
}

ShellRotation IntegralFactory::shell_rotation(int am, SymmetryOperation &so, int pure) const
{
    ShellRotation r(am, so, this, pure);
    return r;
}

SphericalTransformIter *IntegralFactory::spherical_transform_iter(int am, int inv, int subl) const
{
    if (subl != -1)
        throw NOT_IMPLEMENTED_EXCEPTION();

    if (inv) {
        return new SphericalTransformIter(ispherical_transforms_[am]);
    }
    return new SphericalTransformIter(spherical_transforms_[am]);
}

} // namespace psi
