#Upload Pypi package------------------------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/UploadPypi.txt
  COMMAND python3.8 -m twine upload --password ${Password} --username ${Token} --repository pypi /Project/SuPyMode/output/*
  COMMENT "Upload on Pypi")

add_custom_target(UploadPypi DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/UploadPypi.txt)


#Cleaning command---------------------------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt
  COMMAND rm -rf CMakeCache.txt cmake_install.cmake CMakeFiles
  COMMENT "Cleaning cmake output files")

add_custom_target(Clean DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt)
