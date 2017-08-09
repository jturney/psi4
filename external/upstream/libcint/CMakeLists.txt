set(LIBCINT_VER "2.8")

if (${ENABLE_libcint})
    find_package(libcint ${LIBCINT_VER} CONFIG QUIET)

    if (${libcint_FOUND})
        get_property(_loc TARGET libcint::libcint PROPERTY LOCATION)
        message(STATUS "${Cyab}Found libcint${ColourReset}: ${_loc} (found version ${libcint_VERSION})")
        add_library(libcint_external INTERFACE)
    else ()
        include(ExternalProject)
        message(STATUS "Suitable libcint could not be located, ${Magenta}Building libcint${ColourReset} instead.")
        ExternalProject_Add(libcint_external
            GIT_REPOSITORY https://github.com/jturney/libcint
            GIT_TAG master
            UPDATE_COMMAND ""
            CMAKE_ARGS 
                -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
                -DNAMESPACE_INSTALL_INCLUDEDIR=/libxc
                -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                # OpenMP irrelevant
                -DENABLE_XHOST=${ENABLE_XHOST}
                -DBUILD_FPIC=${BUILD_FPIC}
                -DENABLE_GENERIC=${ENABLE_GENERIC}
                -DLIBC_INTERJECT=${LIBC_INTERJECT}
                -DBUILD_TESTING=OFF
            CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            DESTDIR=${CMAKE_BINARY_DIR}/stage)
        set(libcint_DIR ${STAGED_INSTALL_PREFIX}/share/cmake/libcint CACHE PATH "path to internally build libcintConfig.cmake" FORCE)
    endif ()
endif ()
