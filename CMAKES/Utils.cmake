#Upload Pypi package-----------------------------------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Upload.txt
  COMMAND rm -rf dist/*
  COMMAND python3 setup.py bdist_wheel --NewMinor
  COMMAND python3 -m twine upload --password $ENV{PyPiPassword} --username $ENV{PyPiToken} --repository pypi dist/*
  COMMENT "Upload on Pypi"
  COMMAND rm -rf dist *.egg* CMakeCache.txt CMakeFiles *.cmake build)
add_custom_target(Upload DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Upload.txt)



#Cleaning command---------------------------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt
  COMMAND rm -rf CMakeCache.txt cmake_install.cmake CMakeFiles
  COMMENT "Cleaning cmake output files")

add_custom_target(Clean DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt)
