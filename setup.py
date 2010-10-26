import os, sys

import build_system

SAGE_ROOT = os.environ['SAGE_ROOT']

INCLUDES = ['%s/%s/'%(SAGE_ROOT,x) for x in
            ['local/include/csage', 'local/include', 'local/include/python2.6/',
             'devel/sage/sage/ext', 'devel/sage', 'devel/sage/sage/gsl']]

if '-ba' in sys.argv:
    print "Rebuilding all Cython extensions."
    FORCE = True
else:
    FORCE = False

def Extension(*args, **kwds):
    if not kwds.has_key('include_dirs'):
        kwds['include_dirs'] = INCLUDES
    if not kwds.has_key('force'):
        kwds['force'] = FORCE
    return build_system.Extension(*args, **kwds)

ext_modules = [
    Extension("psage.ellff.ellff",
              ["psage/ellff/ellff.pyx",
               "psage/ellff/ell.cpp",
               "psage/ellff/ell_surface.cpp",
               "psage/ellff/euler.cpp",
               "psage/ellff/helper.cpp",
               "psage/ellff/jacobi.cpp",
               "psage/ellff/lzz_pEExtra.cpp",
               "psage/ellff/lzz_pEratX.cpp"],
              language = 'c++'),
    
    Extension("psage.function_fields.function_field_element",
              ["psage/function_fields/function_field_element.pyx"]),
    

]

build_system.cythonize(ext_modules)

build_system.setup(
    name = 'psage',
    version = "10.10.26",
    description = "PSAGE: Software for Arithmetic Geometry",
    author = 'William Stein',
    author_email = 'wstein@gmail.com',
    url = 'http://purple.sagemath.org',
    license = 'GPL v3+',
    packages = ['psage'],
    platforms = ['any'],
    download_url = 'NA',
    ext_modules = ext_modules,
)

