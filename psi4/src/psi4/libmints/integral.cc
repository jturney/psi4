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
#include "psi4/libmints/default/pseudospectral.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/default/ecpint.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/liberd/erd_eri.h"

#ifdef USING_simint
#include "psi4/libmints/simint/siminteri.h"
#endif

#include <libint/libint.h>

#include "psi4/libmints/default/default_integralfactory.h"

using namespace psi;

IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2,
                                 std::shared_ptr<BasisSet> bs3,
                                 std::shared_ptr<BasisSet> bs4)
{
    pImpl_ = std::unique_ptr<DefaultIntegralFactory>(new DefaultIntegralFactory(bs1, bs2, bs3, bs4));
}

IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1)
{
    pImpl_ = std::unique_ptr<DefaultIntegralFactory>(new DefaultIntegralFactory(bs1, bs1, bs1, bs1));
}

IntegralFactory::~IntegralFactory()
{

}

std::shared_ptr<BasisSet> IntegralFactory::basis1() const
{
    return pImpl_->basis1();
}

std::shared_ptr<BasisSet> IntegralFactory::basis2() const
{
    return pImpl_->basis2();
}

std::shared_ptr<BasisSet> IntegralFactory::basis3() const
{
    return pImpl_->basis3();
}

std::shared_ptr<BasisSet> IntegralFactory::basis4() const
{
    return pImpl_->basis4();
}

void IntegralFactory::set_basis(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4)
{
    pImpl_->set_basis(bs1, bs2, bs3, bs4);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_overlap(int deriv)
{
    return pImpl_->ao_overlap(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_overlap(int deriv)
{
    return pImpl_->so_overlap(deriv);
}

std::unique_ptr<ThreeCenterOverlapInt> IntegralFactory::overlap_3c()
{
    return pImpl_->overlap_3c();
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_kinetic(int deriv)
{
    return pImpl_->ao_kinetic(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_kinetic(int deriv)
{
    return pImpl_->so_kinetic(deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_potential(int deriv)
{
    return pImpl_->ao_potential(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_potential(int deriv)
{
    return pImpl_->so_potential(deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_ecp(int deriv)
{
    return pImpl_->ao_ecp(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_ecp(int deriv)
{
	return pImpl_->so_ecp(deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_rel_potential(int deriv)
{
    return pImpl_->ao_rel_potential(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_rel_potential(int deriv)
{
    return pImpl_->so_rel_potential(deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_pseudospectral(int deriv)
{
    return pImpl_->ao_pseudospectral(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_pseudospectral(int deriv)
{
    return pImpl_->so_pseudospectral(deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::electrostatic()
{
    return pImpl_->electrostatic();
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::pcm_potentialint()
{
    return pImpl_->pcm_potentialint();
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_dipole(int deriv)
{
    return pImpl_->ao_dipole(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_dipole(int deriv)
{
    return pImpl_->so_dipole(deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_nabla(int deriv)
{
    return pImpl_->ao_nabla(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_nabla(int deriv)
{
    return pImpl_->so_nabla(deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_angular_momentum(int deriv)
{
    return pImpl_->ao_angular_momentum(deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_angular_momentum(int deriv)
{
    return pImpl_->so_angular_momentum(deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_quadrupole()
{
    return pImpl_->ao_quadrupole();
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_quadrupole()
{
    return pImpl_->so_quadrupole();
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_multipoles(int order)
{
    return pImpl_->ao_multipoles(order);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_efp_multipole_potential(int order)
{
    return pImpl_->ao_efp_multipole_potential(order);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_efp_multipole_potential(int order)
{
    return pImpl_->so_efp_multipole_potential(order);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_multipoles(int order)
{
    return pImpl_->so_multipoles(order);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_traceless_quadrupole()
{
    return pImpl_->ao_traceless_quadrupole();
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_traceless_quadrupole()
{
    return pImpl_->so_traceless_quadrupole();
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::electric_field()
{
    return pImpl_->electric_field();
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erd_eri(int deriv, bool use_shell_pairs)
{
    return pImpl_->erd_eri(deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::eri(int deriv, bool use_shell_pairs)
{
    return pImpl_->eri(deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_eri(double omega, int deriv, bool use_shell_pairs)
{
    return pImpl_->erf_eri(omega, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_complement_eri(double omega, int deriv, bool use_shell_pairs)
{
    return pImpl_->erf_complement_eri(omega, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return pImpl_->f12(cf, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12_scaled(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return pImpl_->f12_scaled(cf, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12_squared(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return pImpl_->f12_squared(cf, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12g12(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return pImpl_->f12g12(cf, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12_double_commutator(std::shared_ptr<CorrelationFactor> cf, int deriv, bool use_shell_pairs)
{
    return pImpl_->f12_double_commutator(cf, deriv, use_shell_pairs);
}

AOShellCombinationsIterator IntegralFactory::shells_iterator()
{
    return AOShellCombinationsIterator(pImpl_->basis1(),
                                       pImpl_->basis2(),
                                       pImpl_->basis3(),
                                       pImpl_->basis4());
}

AOShellCombinationsIterator* IntegralFactory::shells_iterator_ptr()
{
    return new AOShellCombinationsIterator(pImpl_->basis1(),
                                           pImpl_->basis2(),
                                           pImpl_->basis3(),
                                           pImpl_->basis4());
}

AOIntegralsIterator IntegralFactory::integrals_iterator(int p, int q, int r, int s)
{
    return AOIntegralsIterator(pImpl_->basis1()->shell(p),
                               pImpl_->basis2()->shell(q),
                               pImpl_->basis3()->shell(r),
                               pImpl_->basis4()->shell(s));
}

CartesianIter* IntegralFactory::cartesian_iter(int l) const
{
    return new CartesianIter(l);
}

RedundantCartesianIter* IntegralFactory::redundant_cartesian_iter(int l) const
{
    return new RedundantCartesianIter(l);
}

RedundantCartesianSubIter* IntegralFactory::redundant_cartesian_sub_iter(int l) const
{
    return new RedundantCartesianSubIter(l);
}

ShellRotation IntegralFactory::shell_rotation(int am, SymmetryOperation &so, int pure) const
{
    ShellRotation r(am, so, pure);
    return r;
}

SphericalTransformIter* IntegralFactory::spherical_transform_iter(int am, int inv, int subl) const
{
    if (subl != -1)
        throw NOT_IMPLEMENTED_EXCEPTION();

    if (inv) {
        return new SphericalTransformIter(ISphericalTransform::transforms[am]);
    }
    return new SphericalTransformIter(SphericalTransform::transforms[am]);
}
