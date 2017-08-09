set(CINT_VER "2.8")

if (${ENABLE_cint})
    find_package(cint ${LIBCINT_VER} CONFIG QUIET)

    if (${cint_FOUND})
        get_property(_loc TARGET cint::cint PROPERTY LOCATION)
        message(STATUS "${Cyab}Found libcint${ColourReset}: ${_loc} (found version ${cint_VERSION})")
        add_library(cint_external INTERFACE)
    else ()
        include(ExternalProject)
        message(STATUS "Suitable cint could not be located, ${Magenta}Building libcint${ColourReset} instead.")
        ExternalProject_Add(cint_external
            GIT_REPOSITORY https://github.com/jturney/libcint
            GIT_TAG updates
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
