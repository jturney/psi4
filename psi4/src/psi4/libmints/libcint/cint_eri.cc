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

#include "cint_eri.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "cint/cint.h"

namespace psi {

namespace {

int load_basis_info(std::shared_ptr<BasisSet> bs, int& shell_off, std::vector<int>& bas, int& env_off, std::vector<double>& env)
{
    for (int n = 0; n < bs->nshell(); n++) {
        bas[(shell_off + n) * BAS_SLOTS + ATOM_OF] = bs->shell_to_center(n);

        const GaussianShell& shell = bs->shell(n);
        bas[(shell_off + n) * BAS_SLOTS + ANG_OF] = shell.am();
        bas[(shell_off + n) * BAS_SLOTS + NPRIM_OF] = shell.nprimitive();
        bas[(shell_off + n) * BAS_SLOTS + NCTR_OF] = 1;


        // load in exponents
        bas[(shell_off + n) * BAS_SLOTS + PTR_EXP] = env_off;
        for (int exp = 0; exp < shell.nprimitive(); exp++) {
            env[env_off++] = shell.exp(exp);
        }

        // load in coefficients
        bas[(shell_off + n) * BAS_SLOTS + PTR_COEFF] = env_off;
        for (int coef = 0; coef < shell.nprimitive(); coef++) {
            env[env_off++] = shell.original_coef(coef) * CINTgto_norm(bas[(shell_off + n) * BAS_SLOTS + ANG_OF], env[bas[(shell_off + n) * BAS_SLOTS + PTR_EXP] + coef]);
        }
    }

    return bs->nshell();
}

}

CINTERI::CINTERI(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4, int deriv)
        : TwoBodyAOInt(bs1, bs2, bs3, bs4, deriv)
{
    size_t max_cart = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                      INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    try {
        tformbuf_ = new double[max_cart];
    }
    catch (std::bad_alloc& e) {
        outfile->Printf("Error allocating tformbuf_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(tformbuf_, 0, sizeof(double) * max_cart);


    try {
        target_full_ = new double[max_cart];
	target_ = target_full_;
    }
    catch (std::bad_alloc& e) {
        outfile->Printf("Error allocating target_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(target_, 0, sizeof(double) * max_cart);

    try {
        source_full_ = new double[max_cart];
	source_ = source_full_;
    }
    catch (std::bad_alloc& e) {
        outfile->Printf("Error allocating source_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(source_, 0, sizeof(double) * max_cart);

    std::shared_ptr<Molecule> molecule = bs1->molecule();

    // Figure out how many unique basis sets we have.
    int nshell = bs1->nshell();
    int envsize = bs1->nprimitive();
    if (bs2 != bs1) {
        nshell += bs2->nshell();
        envsize += bs2->nprimitive();
    }
    if (bs3 != bs2 && bs3 != bs1) {
        nshell += bs3->nshell();
        envsize += bs3->nprimitive();
    }
    if (bs4 != bs3 && bs4 != bs2 && bs4 != bs1) {
        nshell += bs4->nshell();
        envsize += bs4->nprimitive();
    }
    envsize = envsize * 2 + PTR_ENV_START + molecule->natom()*3;
    // Populate the data structures for cint
    atm_.resize(molecule->natom() * ATM_SLOTS, 0);
    bas_.resize(nshell * BAS_SLOTS, 0);
    env_.resize(envsize);

    int off = PTR_ENV_START;

    // Load in atom information
    for (int i = 0; i < molecule->natom(); i++) {
        atm_[i * ATM_SLOTS + CHARGE_OF] = static_cast<int>(molecule->Z(i));
        atm_[i * ATM_SLOTS + PTR_COORD] = off;

        env_[off + 0] = molecule->x(i);
        env_[off + 1] = molecule->y(i);
        env_[off + 2] = molecule->z(i);

        off += 3;
    }

    std::shared_ptr<BasisSet> sets[] = {bs1, bs2, bs3, bs4};
    basis_shell_start_[0] = 0;

    // load bs1
    nshell = load_basis_info(bs1, basis_shell_start_[0], bas_, off, env_);

    // check if bs2 is equivalent to bs1, if so, just use bs1, otherwise load it
    if (bs2 == bs1) {
        basis_shell_start_[1] = basis_shell_start_[0];
    } else {
        basis_shell_start_[1] = basis_shell_start_[0] + nshell;
        nshell = load_basis_info(bs2, basis_shell_start_[1], bas_, off, env_);
    }

    // check if bs3 is equivalent to something
    if (bs3 == bs2) {
        basis_shell_start_[2] = basis_shell_start_[1];
    } else if (bs3 == bs1) {
        basis_shell_start_[2] = basis_shell_start_[0];
    } else {
        basis_shell_start_[2] = basis_shell_start_[1] + nshell;
        nshell = load_basis_info(bs3, basis_shell_start_[2], bas_, off, env_);
    }

    // check if bs4 is equivalent to something
    if (bs4 == bs3) {
        basis_shell_start_[3] = basis_shell_start_[2];
    } else if (bs4 == bs2) {
        basis_shell_start_[3] = basis_shell_start_[1];
    } else if (bs4 == bs1) {
        basis_shell_start_[3] = basis_shell_start_[0];
    } else {
        basis_shell_start_[3] = basis_shell_start_[2] + nshell;
        nshell = load_basis_info(bs4, basis_shell_start_[3], bas_, off, env_);
    }

    nshell_ = basis_shell_start_[3] + nshell;
    natom_ = molecule->natom();

    // All data is loaded.
}

size_t CINTERI::compute_shell(int p, int q, int r, int s)
{
    int shells[4] = {p+basis_shell_start_[0],
                     q+basis_shell_start_[1],
                     r+basis_shell_start_[2],
                     s+basis_shell_start_[3]};

    int np = INT_NCART(bas_[shells[0] + ANG_OF]);
    int nq = INT_NCART(bas_[shells[1] + ANG_OF]);
    int nr = INT_NCART(bas_[shells[2] + ANG_OF]);
    int ns = INT_NCART(bas_[shells[3] + ANG_OF]);

    const auto& shell1 = original_bs1_->shell(p);
    const auto& shell2 = original_bs2_->shell(q);
    const auto& shell3 = original_bs3_->shell(r);
    const auto& shell4 = original_bs4_->shell(s);

    bool do_cart = force_cartesian_ || (shell1.is_cartesian() &&
                                        shell2.is_cartesian() &&
                                        shell3.is_cartesian() &&
                                        shell4.is_cartesian());

    int retval = 0;
    if (do_cart) {
        retval = cint2e_cart(target_, shells, atm_.data(), natom_, bas_.data(), nshell_, env_.data(), nullptr);
        return retval;
    }

    retval = cint2e_cart(source_, shells, atm_.data(), natom_, bas_.data(), nshell_, env_.data(), nullptr);
    if (0 != retval) {
        pure_transform(p, q, r, s, 1, false);
        return retval;
    } else {
        // zero memory
        ::memset(target_, 0, sizeof(double) * np * nq * nr * ns);
        return 0;
    }
}

size_t CINTERI::compute_shell(const AOShellCombinationsIterator& shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t CINTERI::compute_shell_deriv1(int, int, int, int)
{
    throw PSIEXCEPTION("ERROR - SIMINT CANNOT HANDLE DERIVATIVES\n");
}

size_t CINTERI::compute_shell_deriv2(int, int, int, int)
{
    throw PSIEXCEPTION("ERROR - SIMINT CANNOT HANDLE DERIVATIVES\n");
}

} // psi namespace
