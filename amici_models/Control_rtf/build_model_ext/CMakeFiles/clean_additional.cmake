# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Control_rtf.dir/Control_rtfPYTHON_wrap.cxx"
  "swig/Control_rtf.py"
  )
endif()
