set(sources_list local.cc
        onebody.cc
        orbitalspace.cc
        osrecur.cc
        maketab.cc
        oeprop.cc
        writer.cc
        transform.cc
        sieve.cc
        multipolesymmetry.cc
        shellrotation.cc
        deriv.cc
        integralparameters.cc
        twobody.cc
        vector.cc
        sobasis.cc
        cartesianiter.cc
        basisset.cc
        wavefunction.cc
        irrep.cc
        eribase.cc
        fjt.cc
        chartab.cc
        corrtab.cc
        symop.cc
        benchmark.cc
        get_writer_file_prefix.cc
        3coverlap.cc
        petitelist.cc
        solidharmonics.cc
        #orthog.cc
        sointegral.cc
        extern.cc
        factory.cc
        psimath.cc
        dimension.cc
        molecule.cc
        intvector.cc
        mintshelper.cc
        coordentry.cc
        pseudospectral.cc
        integral.cc
        matrix.cc
        gshell.cc
        integraliter.cc
        pointgrp.cc
        rep.cc
        cdsalclist.cc
        bessel.cc
        gaussquad.cc
        ecpint.cc

        default/angularmomentum.cc
        default/dipole.cc
        default/efpmultipolepotential.cc
        default/electricfield.cc
        default/electrostatic.cc
        default/eri.cc
        default/kinetic.cc
        default/multipoles.cc
        default/nabla.cc
        default/overlap.cc
        default/potential.cc
        default/potentialint.cc
        default/quadrupole.cc
        default/rel_potential.cc
        default/tracelessquadrupole.cc
        default/x2cint.cc

        liberd/erd_eri.cc
        )

if (TARGET simint::simint)
    list(APPEND sources_list simint/siminteri.cc)
    add_definitions("-DUSING_simint")
endif ()
psi4_add_module(lib mints sources_list iwl options psi4util trans efp_solver)

target_link_libraries(mints PRIVATE Libint::libint)
if (TARGET dkh::dkh)
    target_link_libraries(mints PRIVATE dkh::dkh)
endif ()
if (TARGET libefp::efp)
    target_link_libraries(mints PRIVATE libefp::efp)
endif ()
if (TARGET erd::erd)
    target_link_libraries(mints PRIVATE erd::erd)
endif ()
if (TARGET simint::simint)
    target_link_libraries(mints PRIVATE simint::simint)
endif ()
