# Meson file for building the preprocessor in a standalone fashion.
# It is not used when building Dynare as a whole.

project('dynare-preprocessor', 'cpp',
        version : '6-unstable',
        default_options : [ 'cpp_std=gnu++20', 'warning_level=2' ])

add_global_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language : 'cpp')

cpp_compiler = meson.get_compiler('cpp')

subdir('src')
subdir('doc')
