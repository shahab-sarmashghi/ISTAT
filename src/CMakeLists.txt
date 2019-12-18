add_subdirectory(lib)
set(SOURCE_FILES main_chr.cpp)

add_executable(istat ${SOURCE_FILES})
target_link_libraries(istat DP PB FUNC PREP)
install(TARGETS istat DESTINATION ${ISTAT_INSTALL_BIN_DIR})
