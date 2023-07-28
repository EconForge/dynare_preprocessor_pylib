latexmk = find_program('latexmk')

# We have to set TEXINPUTS because the current directory is not the source
# directory when latexmk is invoked (and using the -cd option in combination
# with -outdir/-auxdir does not work in all cases because @OUTDIR@ and
# @PRIVATE_DIR@ can be relative paths)

# Without the -g flag, latexmk remembers a previous build failure and will refuse
# to recompile even if the error has been fixed in the TeX source

# The \graphicspath{} command does not compute directories relative to TEXINPUTS,
# so add these manually
logos_dir = meson.current_source_dir() / 'logos'

custom_target('macroprocessor.pdf',
              output : 'macroprocessor.pdf',
              input : 'macroprocessor/macroprocessor.tex',
              command : [ latexmk, '-pdf', '-g', '-outdir=@OUTDIR@', '-auxdir=@PRIVATE_DIR@', '@INPUT@'],
              env : { 'TEXINPUTS': meson.current_source_dir() + '/macroprocessor:' + logos_dir + ':' },
              install : true,
              install_dir : 'share/doc/dynare')

custom_target('preprocessor.pdf',
              output : 'preprocessor.pdf',
              input : 'preprocessor/preprocessor.tex',
              command : [ latexmk, '-pdf', '-g', '-outdir=@OUTDIR@', '-auxdir=@PRIVATE_DIR@', '@INPUT@'],
              env : { 'TEXINPUTS': meson.current_source_dir() + '/preprocessor:' + logos_dir + ':' },
              install : true,
              install_dir : 'share/doc/dynare')