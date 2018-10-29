# Specification file for the CMake build system
#
# Laurence R. McGLashan
# 22 June 2011

project (CHEMKINREADER)

file(GLOB_RECURSE INCS "include/*.h")
source_group("Header Files" FILES ${INCS})

add_library(chemkinReader 
            source/stringFunctions.cpp
            source/element.cpp 
            source/transport.cpp
            source/species.cpp
            source/thermo.cpp 
            source/reaction.cpp
            source/thermoParser.cpp 
            source/transportParser.cpp
            source/reactionParser.cpp
            source/chemkinReader.cpp
			source/phase.cpp
			${INCS})

include_directories(include
                    ${Boost_INCLUDE_DIRS}
                    ../../sprogc/include
                    ../../utils/include)

target_link_libraries(chemkinReader ${Boost_LIBRARIES})

install (TARGETS chemkinReader EXPORT chemkinReader DESTINATION lib)
install (DIRECTORY include DESTINATION .)
