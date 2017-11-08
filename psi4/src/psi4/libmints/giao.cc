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

#include "giao.h"
#include "molecule.h"
#include "basisset.h"
#include "integral.h"
#include "vector.h"

namespace psi {

GIAO_OverlapInt::GIAO_OverlapInt(std::vector<SphericalTransform> &spherical_transforms,
                                 std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2,
                                 int nderiv)
        : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv),
          overlap_recur_(bs1->max_am() + 2, bs2->max_am() + 2)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1 + 1) * (maxam1 + 2) / 2;
    int maxnao2 = (maxam2 + 1) * (maxam2 + 2) / 2;

    // Increase buffer size to handle x, y, and z components
    if (deriv_ == 0) {
        buffer_ = new double[3 * maxnao1 * maxnao2];
        set_chunks(3);
    }
}

GIAO_OverlapInt::~GIAO_OverlapInt()
{
    delete[] buffer_;
}

void GIAO_OverlapInt::compute_pair(const GaussianShell &s1, const GaussianShell &s2)
{
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    int ydisp = INT_NCART(am1) * INT_NCART(am2);
    int zdisp = ydisp + INT_NCART(am1) * INT_NCART(am2);

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 3 * INT_NCART(am1) * INT_NCART(am2) * sizeof(double));

    double **x = overlap_recur_.x();
    double **y = overlap_recur_.y();
    double **z = overlap_recur_.z();

    for (int p1 = 0; p1 < nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2 = 0; p2 < nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0 / gamma;

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1 * A[0] + a2 * B[0]) * oog;
            P[1] = (a1 * A[1] + a2 * B[1]) * oog;
            P[2] = (a1 * A[2] + a2 * B[2]) * oog;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double over_pf = exp(-a1 * a2 * AB2 * oog) * sqrt(M_PI * oog) * M_PI * oog * c1 * c2;

            // Do recursion
            overlap_recur_.compute(PA, PB, gamma, am1 + 1, am2 + 1);

            ao12 = 0;
            for (int ii = 0; ii <= am1; ii++) {
                int l1 = am1 - ii;
                for (int jj = 0; jj <= ii; jj++) {
                    int m1 = ii - jj;
                    int n1 = jj;

                    for (int kk = 0; kk <= am2; kk++) {
                        int l2 = am2 - kk;
                        for (int ll = 0; ll <= kk; ll++) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            double x00 = x[l1][l2], y00 = y[m1][m2], z00 = z[n1][n2];
                            double x10 = x[l1 + 1][l2], y10 = y[m1 + 1][m2], z10 = z[n1 + 1][n2];

                            // Compute dipole moment components
                            double DAx = (x10 + x00 * (A[0] - origin_[0])) * y00 * z00 * over_pf;
                            double DAy = x00 * (y10 + y00 * (A[1] - origin_[1])) * z00 * over_pf;
                            double DAz = x00 * y00 * (z10 + z00 * (A[2] - origin_[2])) * over_pf;

                            double SBx = (A[1] - B[1]) * DAz - (A[2] - B[2]) * DAy;
                            double SBy = (A[2] - B[2]) * DAx - (A[0] - B[0]) * DAz;
                            double SBz = (A[0] - B[0]) * DAy - (A[1] - B[1]) * DAx;

                            buffer_[ao12] += SBx;
                            buffer_[ao12+ydisp] += SBy;
                            buffer_[ao12+zdisp] += SBz;

                            ao12++;
                        }
                    }
                }
            }
        }
    }
}

} // namespace psi
