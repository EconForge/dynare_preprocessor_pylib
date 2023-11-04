# Meson file for building the preprocessor in a standalone fashion.
# It is not used when building Dynare as a whole.

project('dynare-preprocessor', 'cpp',
        version : '6-unstable',
        default_options : [ 'cpp_std=c++20', 'warning_level=0' ],
        meson_version : '>=0.64.0')

add_global_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language : 'cpp')

if get_option('warning_level').to_int() >= 2
  add_global_arguments('-Wold-style-cast', language : 'cpp')
endif

cpp_compiler = meson.get_compiler('cpp')

subdir('src')
subdir('doc')
