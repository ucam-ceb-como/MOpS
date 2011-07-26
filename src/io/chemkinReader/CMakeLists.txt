# Specification file for the CMake build system
#
# Laurence R. McGLashan
# 22 June 2011

project (CHEMKINREADER)

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
            source/chemkinReader.cpp)

include_directories(include
                    ${Boost_INCLUDE_DIRS}
                    ../../sprogc/include
                    ../../utils/include)

target_link_libraries(chemkinReader ${Boost_LIBRARIES})
