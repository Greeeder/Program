# Install script for directory: D:/JOSEP/OpenWAM/Source

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files/Project")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("D:/JOSEP/OpenWAM/Units/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/CheckXML/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Fluids/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Solids/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/1DPipes/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Act/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Boundaries/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Concentric Pipe/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Connections/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Control/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/DataIO/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Labels/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Integrable/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Math_wam/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/DPF/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Engine/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Extern/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/ODModels/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Output/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Q2DTurbo/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Turbocompressor/cmake_install.cmake")
  include("D:/JOSEP/OpenWAM/Wrappers/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "D:/JOSEP/OpenWAM/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
