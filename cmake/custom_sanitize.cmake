if(ENABLE_ASAN)
    if((CMAKE_CXX_COMPILER_ID MATCHES Clang) OR
        (CMAKE_CXX_COMPILER_ID MATCHES AppleClang) OR
        (CMAKE_CXX_COMPILER_ID MATCHES GNU))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    endif()
endif()
