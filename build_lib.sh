#!/bin/bash

export PREFIX=$(python -c "import sys;print(sys.prefix)") 
export BOOST_ROOT=$PREFIX

#cp $RECIPE_DIR/meson.build src/meson.build

meson setup --prefix=$PREFIX --bindir=$PREFIX/bin --libdir=$PREFIX/lib --includedir=$PREFIX/include \
    --buildtype=release build_lib \
    -Dcpp_args="-w  -Wno-enum-constexpr-conversion -I${PREFIX}/include/pybind11"  \
    -Dcpp_link_args="-w  -Wno-enum-constexpr-conversion -I${PREFIX}/include/pybind11" \
    -Dcpp_link_args="-pthread" \
    -Dbuild_library="enabled"
    
meson compile -C build_lib
meson install -C build_lib #--destdir="../

#rm $PREFIX/bin/python
