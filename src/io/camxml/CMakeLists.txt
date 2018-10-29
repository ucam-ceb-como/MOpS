# Specification file for the CMake cross platform build system

project (CAMXML)

file(GLOB_RECURSE INCS "include/*.h")
source_group("Header Files" FILES ${INCS})

add_library(camxml source/stringf.cpp
                   source/unicodef.cpp
                   source/xmlattribute.cpp
                   source/xmldocument.cpp
                   source/xmlelement.cpp
				   ${INCS})

include_directories (include ../comostrings/include)

# Anyone linking against this library will also need to link against the comostrings
# library, because it is references from camxml.
target_link_libraries(camxml comostrings)

install (TARGETS camxml EXPORT camxml DESTINATION lib)
install (DIRECTORY include DESTINATION .)
