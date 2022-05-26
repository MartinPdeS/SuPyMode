#Building EigenSolver---------------------------------------------------------------------
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/SuPyMode/bin)
pybind11_add_module(EigenSolver MODULE SuPyMode/includes/interface.cpp )
target_link_libraries(EigenSolver PRIVATE ${PYTHON_LIBRARIES})
set_target_properties(EigenSolver PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(EigenSolver PROPERTIES OUTPUT_NAME "EigenSolver")
target_compile_options (EigenSolver PRIVATE -O3)
