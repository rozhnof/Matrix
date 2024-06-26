set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -Wpedantic -g -fno-omit-frame-pointer)
add_compile_options(-Wno-language-extension-token)
add_compile_options(-Wno-error=unused-command-line-argument)
add_compile_options(-gdwarf-4)
# add_compile_options(-Werror)

message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")