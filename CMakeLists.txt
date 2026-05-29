cmake_minimum_required (VERSION 3.15)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("CMakeSimply1")

# Определяем файлы исходного кода
set(SOURCES
    main.cpp
    NoShared/rotate.cpp
)

set(HEADERS
  NoShared/rotate.h
  NoShared/emu_vector_types.h
)

add_executable (${PROJECT_NAME} ${SOURCES} ${HEADERS})

