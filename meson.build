# Meson file for building the preprocessor in a standalone fashion.
# It is not used when building Dynare as a whole.

project('dynare-preprocessor', 'cpp',
        version : '8-unstable',
        license : 'GPL-3.0-or-later',
        license_files : 'COPYING',
        # NB: update C++ standard in .clang-format whenever the following is modified
        default_options : [ 'cpp_std=gnu++20,c++20', 'warning_level=2' ],
        meson_version : '>=1.3.0')

add_global_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language : 'cpp')

cpp_compiler = meson.get_compiler('cpp')

if get_option('warning_level').to_int() >= 2 and cpp_compiler.get_id() != 'msvc' and cpp_compiler.get_id() != 'clang-cl'
  add_global_arguments('-Wold-style-cast', language : 'cpp')
endif

if cpp_compiler.get_id() == 'msvc'
  add_global_arguments('/FI' + (meson.project_source_root() / 'src' / 'Portability.hh').replace('\\', '/'), language : 'cpp')
  add_global_arguments('/DYY_NO_UNISTD_H', language : 'cpp')
  add_global_arguments('/DNOMINMAX', language : 'cpp')
  add_global_arguments('/DWIN32_LEAN_AND_MEAN', language : 'cpp')
endif

subdir('src')
subdir('doc')
subdir('tests')
