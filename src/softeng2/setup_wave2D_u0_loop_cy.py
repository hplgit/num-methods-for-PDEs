from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

cymodule = 'wave2D_u0_loop_cy'
setup(
  name=cymodule,
  ext_modules=[Extension(cymodule, [cymodule + '.pyx'],
                         libraries=[], # C libs to link with
                         )],
  cmdclass={'build_ext': build_ext},
)
