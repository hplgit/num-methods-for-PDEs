#!/bin/sh
# Compile extension modules for the loop


# Cython (easier with pyximport)
module=wave2D_u0_loop_cy
rm -f $module.so
python setup_${module}.py build_ext --inplace   # compile
python -c "import $module"                      # test
if [ $? -eq 0 ]; then                           # success?
  echo "Cython module $module successfully built"
else
  echo "Building Cython module $module failed"
  exit 1
fi

# Fortran
module=wave2D_u0_loop_f77
rm -f $module.so
f2py -m $module -h ${module}.pyf --overwrite-signature $module.f
f2py -c $module.pyf --fcompiler=gfortran --build-dir tmp_build_f77 \
     -DF2PY_REPORT_ON_ARRAY_COPY=1 $module.f
python -c "
import $module as m
print m.__doc__
print m.advance.__doc__"
if [ $? -eq 0 ]; then    # success?
  echo "Fortran module $module successfully built"
else
  echo "Building Fortran module $module failed"
  exit 1
fi

# C via f2py
module=wave2D_u0_loop_c_f2py
rm -f $module.so
f2py -m $module -h ${module}.pyf --overwrite-signature \
      ${module}_signature.f
f2py -c $module.pyf --build-dir tmp_build_c \
     -DF2PY_REPORT_ON_ARRAY_COPY=1 wave2D_u0_loop_c.c
python -c "
import $module as m
print m.__doc__
print m.advance.__doc__"
if [ $? -eq 0 ]; then    # success?
  echo "C module $module successfully built"
else
  echo "Building C module $module failed"
  exit 1
fi

# Cython interface to C code
module=wave2D_u0_loop_c_cy
rm -f $module.so
python setup_${module}.py build_ext --inplace   # compile
python -c "import $module"                      # test
if [ $? -eq 0 ]; then                           # success?
  echo "Cython module $module successfully built"
else
  echo "Building Cython module $module failed"
  exit 1
fi

