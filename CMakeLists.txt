cmake_minimum_required(VERSION 3.10)

# Nombre del proyecto y lenguaje
project(c-psd-estimator C)

# Estándar de C
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Nombre del ejecutable
set(EXECUTABLE_NAME "c-psd-estimator")

# Archivos fuente
set(SOURCE_FILES
    main.c
    Modules/utils.c
    Modules/psd-estimators.c
)

# Archivos de cabecera (opcional, solo para IDEs)
set(HEADER_FILES
    Modules/utils.h
    Modules/psd-estimators.h
    Modules/datatypes.h
)

# Crear ejecutable
add_executable(${EXECUTABLE_NAME} ${SOURCE_FILES} ${HEADER_FILES})

# Incluir directorio de headers locales
target_include_directories(${EXECUTABLE_NAME} PUBLIC
    ${PROJECT_SOURCE_DIR}/Modules
)

# -----------------------------
# Dependencias externas
# -----------------------------

# FFTW3
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${FFTW_LIBRARIES} m)

# HackRF (opcional, se activa con -DHACKRF=ON)
option(HACKRF "Enable HackRF support" ON)
if(HACKRF)
    pkg_check_modules(LIBHACKRF REQUIRED libhackrf)
    target_include_directories(${EXECUTABLE_NAME} PRIVATE ${LIBHACKRF_INCLUDE_DIRS})
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${LIBHACKRF_LIBRARIES})
endif()

# -----------------------------
# Configuración de compilación
# -----------------------------
# Guardar el binario en la raíz del proyecto
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# Flags extra (warnings + optimización por defecto)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(${EXECUTABLE_NAME} PRIVATE -Wall -Wextra -g)
else()
    target_compile_options(${EXECUTABLE_NAME} PRIVATE -Wall -O2)
endif()
