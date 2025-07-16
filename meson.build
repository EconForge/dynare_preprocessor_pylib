project('dynare-preprocessor', 'cpp',
		version : '7-unstable',
		# NB: update C++ standard in .clang-format whenever the following is modified
		default_options : [ 'cpp_std=gnu++20', 'warning_level=2' ],
		meson_version : '>=1.1')

build_library = get_option('build_library')

# if build_library.disabled()

	# Meson file for building the preprocessor in a standalone fashion.
	# It is not used when building Dynare as a whole.

add_global_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language : 'cpp')

if get_option('warning_level').to_int() >= 2
add_global_arguments('-Wold-style-cast', language : 'cpp')
endif

cpp_compiler = meson.get_compiler('cpp')

subdir('src')
subdir('doc')

# endif