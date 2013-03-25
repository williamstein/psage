#################################################################################
#
# (c) Copyright 2010 William Stein
#
#  This file is part of PSAGE
#
#  PSAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PSAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################


import os, sys


if sys.maxint != 2**63 - 1:
    print "*"*70
    print "The PSAGE library only works on 64-bit computers.  Terminating build."
    print "*"*70
    sys.exit(1)


import build_system

SAGE_ROOT = os.environ['SAGE_ROOT']
SAGE_LOCAL = os.environ['SAGE_LOCAL']

INCLUDES = ['%s/%s/'%(SAGE_ROOT,x) for x in
             ('devel/sage/sage/ext', 'devel/sage', 'devel/sage/sage/gsl')] \
         + ['%s/%s/'%(SAGE_LOCAL,x) for x in
             ('include/csage', 'include', 'include/python')]
if '-ba' in sys.argv:
    print "Rebuilding all Cython extensions."
    sys.argv.remove('-ba')
    FORCE = True
else:
    FORCE = False

def Extension(*args, **kwds):
    if not kwds.has_key('include_dirs'):
        kwds['include_dirs'] = INCLUDES
    else:
        kwds['include_dirs'] += INCLUDES
    if not kwds.has_key('force'):
        kwds['force'] = FORCE

    # Disable warnings when running GCC step -- cython has already parsed the code and
    # generated any warnings; the GCC ones are noise.
    if not kwds.has_key('extra_compile_args'):
        kwds['extra_compile_args'] = ['-w']
    else:
        kwds['extra_compile_args'].append('-w')

    E = build_system.Extension(*args, **kwds)
    E.libraries = ['csage'] + E.libraries
    return E


numpy_include_dirs = [os.path.join(SAGE_LOCAL,
                                   'lib/python/site-packages/numpy/core/include')]

ext_modules = [
# Remove until the database is rewritten to not use ZODB (which was removed from Sage 5.8)
#    Extension("psage.ellff.ellff",
#              ["psage/ellff/ellff.pyx",
#               "psage/ellff/ell.cpp",
#               "psage/ellff/ell_surface.cpp",
#               "psage/ellff/euler.cpp",
#               "psage/ellff/helper.cpp",
#               "psage/ellff/jacobi.cpp",
#               "psage/ellff/lzz_pEExtra.cpp",
#               "psage/ellff/lzz_pEratX.cpp"],
#              language = 'c++'),

    Extension("psage.function_fields.function_field_element",
              ["psage/function_fields/function_field_element.pyx"]),

    Extension("psage.modform.jacobiforms.jacobiformd1nn_fourierexpansion_cython",
              ["psage/modform/jacobiforms/jacobiformd1nn_fourierexpansion_cython.pyx"]),
    
    Extension("psage.modform.paramodularforms.siegelmodularformg2_misc_cython",
              ["psage/modform/paramodularforms/siegelmodularformg2_misc_cython.pyx"]),

    Extension("psage.modform.paramodularforms.siegelmodularformg2_fourierexpansion_cython",
              ["psage/modform/paramodularforms/siegelmodularformg2_fourierexpansion_cython.pyx"]),

    Extension("psage.modform.paramodularforms.siegelmodularformg2vv_fegenerators_cython",
              ["psage/modform/paramodularforms/siegelmodularformg2vv_fegenerators_cython.pyx"]),

    Extension("psage.modform.paramodularforms.paramodularformd2_fourierexpansion_cython",
              ["psage/modform/paramodularforms/paramodularformd2_fourierexpansion_cython.pyx"]),

    Extension("psage.modform.siegel.fastmult",
              ["psage/modform/siegel/fastmult.pyx"]),

    Extension('psage.modform.maass.mysubgroups_alg',
              ['psage/modform/maass/mysubgroups_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),

    Extension('psage.modform.maass.maass_forms_alg',
              ['psage/modform/maass/maass_forms_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),

    Extension('psage.modform.maass.lpkbessel',
              ['psage/modform/maass/lpkbessel.pyx']),

    Extension("psage.modform.rational.modular_symbol_map",
              ["psage/modform/rational/modular_symbol_map.pyx"]),

    Extension("psage.modform.rational.padic_elliptic_lseries_fast",
              ["psage/modform/rational/padic_elliptic_lseries_fast.pyx"]),

    Extension("psage.modform.hilbert.sqrt5.sqrt5_fast",
              ["psage/modform/hilbert/sqrt5/sqrt5_fast.pyx"],
              libraries = ['ntl', 'gmp'],
              language = 'c++'),

    Extension("psage.ellcurve.lseries.sqrt5",
              ["psage/ellcurve/lseries/sqrt5.pyx"],
              libraries = ['ntl', 'gmp'],
              language = 'c++'),

    Extension("psage.ellcurve.lseries.helper",
              ["psage/ellcurve/lseries/helper.pyx"]),

    Extension('psage.ellcurve.galrep.wrapper',
              sources = ['psage/ellcurve/galrep/wrapper.pyx', 'psage/ellcurve/galrep/galrep.c'],
              libraries = ['gmp']),

    Extension('psage.ellcurve.minmodel.sqrt5',
              sources = ['psage/ellcurve/minmodel/sqrt5.pyx'],
              libraries = ['gmp']),

    Extension('psage.rh.mazur_stein.game',
              sources = ['psage/rh/mazur_stein/game.pyx']),

    Extension('psage.rh.mazur_stein.book_cython',
              sources = ['psage/rh/mazur_stein/book_cython.pyx']),

    Extension("psage.ellcurve.lseries.fast_twist",
              ["psage/ellcurve/lseries/fast_twist.pyx"],
              libraries = ['gsl']),

    Extension("psage.ellcurve.lseries.aplist_sqrt5",
              ["psage/ellcurve/lseries/aplist_sqrt5.pyx"],
              language = 'c++'),

    Extension("psage.number_fields.sqrt5.prime",
              ["psage/number_fields/sqrt5/prime.pyx"],
              libraries = ['pari']),

    Extension("psage.modform.rational.special_fast",
              ["psage/modform/rational/special_fast.pyx", SAGE_ROOT + "/devel/sage/sage/libs/flint/fmpq_poly.c"],
              libraries = ['gmp', 'flint'],
              language = 'c++',
              include_dirs = [SAGE_LOCAL + '/include/FLINT/', SAGE_ROOT + '/devel/sage/sage/libs/flint/'],
              extra_compile_args = ['-std=c99']),

    Extension("psage.ellcurve.xxx.rankbound",
              sources = [   'psage/ellcurve/xxx/rankbound.pyx',
                            'psage/ellcurve/xxx/rankbound_.cc',
                            'psage/ellcurve/xxx/mathlib.cc',
                            'psage/libs/smalljac/wrapper_g1.c'],
              libraries = ['gmp', 'm'],
              include_dirs = ['psage/libs/smalljac/'],
              language = 'c'
              )
]

for g in [1, 2]:
    e = Extension('psage.libs.smalljac.wrapper%s'%g,
                  sources = ['psage/libs/smalljac/wrapper%s.pyx'%g,
                             'psage/libs/smalljac/wrapper_g%s.c'%g],
                  libraries = ['gmp', 'm'])
    ext_modules.append(e)
## Fredrik Stroemberg: my additional modules.
my_extensions = [
    Extension('psage.rings.mpc_extras',
              sources = ['psage/rings/mpc_extras.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp']),

    Extension('psage.modform.maass.inc_gamma',
              ['psage/modform/maass/inc_gamma.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp']),
    Extension('psage.modform.maass.mysubgroups_alg',
              ['psage/modform/maass/mysubgroups_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc']),

    Extension('psage.modform.maass.permutation_alg',
              ['psage/modform/maass/permutation_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc']),
    
    Extension('psage.modform.maass.maass_forms_parallel_alg',
              ['psage/modform/maass/maass_forms_parallel_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp']),
#              include_dirs = numpy_include_dirs),

    Extension('psage.modform.maass.pullback_algorithms',
              ['psage/modform/maass/pullback_algorithms.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp']),
#              include_dirs = numpy_include_dirs),
    Extension('psage.modform.maass.linear_system',
              ['psage/modform/maass/linear_system.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),
    Extension('psage.modform.maass.automorphic_forms_alg',
              ['psage/modform/maass/automorphic_forms_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
#              include_dirs = numpy_include_dirs,
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp']),
<<<<<<< Updated upstream
    Extension('psage.modform.maass.hilbert_modular_group_alg',
              ['psage/modform/maass/hilbert_modular_group_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),
=======

    Extension('psage.modform.hilbert.hn_class',
              ['psage/modform/hilbert/hn_class.pyx'],
              libraries = ['m']),
#              include_dirs = numpy_include_dirs),
    Extension('psage.modform.hilbert.test_nf',
              sources=['psage/modform/hilbert/test_nf.pyx'],
              language="c++"),
    
    Extension('psage.modform.hilbert.hilbert_modular_group_alg',
              ['psage/modform/hilbert/hilbert_modular_group_alg.pyx'],
              language="c++"),
#              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++", "ntl"],
              #libraries = ['ntl','gmp','gmpxx','mpfr','mpc','csage','m','pari','stdc++','flint'],
#              extra_compile_args = ['-std=c99']),


#              language="c++",
#              libraries = ["csage", "flint", "gmp", "gmpxx", "m", "stdc++", "ntl"],
#              include_dirs = numpy_include_dirs),

>>>>>>> Stashed changes
    Extension('psage.zfunctions.selberg_z_alg',
              ['psage/zfunctions/selberg_z_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
#              include_dirs = numpy_include_dirs,
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp']),
              
    
    Extension('psage.modform.maass.vv_harmonic_weak_maass_forms_alg',
              ['psage/modform/maass/vv_harmonic_weak_maass_forms_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),
    Extension('psage.modform.maass.maass_forms_phase2',
                  sources=['psage/modform/maass/maass_forms_phase2.pyx'],
                  libraries = ['m','gmp','mpfr','mpc'],
                  include_dirs = numpy_include_dirs),
    Extension('psage.modform.maass.lpkbessel',
              ['psage/modform/maass/lpkbessel.pyx']),
    Extension('psage.modules.vector_complex_dense',
              sources = ['psage/modules/vector_complex_dense.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),
    Extension('psage.modules.vector_real_mpfr_dense',
              sources = ['psage/modules/vector_real_mpfr_dense.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),
    
    Extension('psage.modules.weil_module_alg',
              sources = ['psage/modules/weil_module_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),
    Extension('psage.matrix.matrix_complex_dense',
              sources = ['psage/matrix/matrix_complex_dense.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp'],
              include_dirs = numpy_include_dirs),

    Extension('psage.matrix.linalg_complex_dense',
              sources = ['psage/matrix/linalg_complex_dense.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp']),
    Extension('psage.modform.maass.poincare_series_alg',
                  ['psage/modform/maass/poincare_series_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc']),
    
    Extension('psage.modform.maass.poincare_series_alg_vv',
          ['psage/modform/maass/poincare_series_alg_vv.pyx'],
          libraries = ['m','gmp','mpfr','mpc']),
        
    Extension('psage.modform.maass.eisenstein_series',
              ['psage/modform/maass/eisenstein_series.pyx'],
              libraries = ['m','gmp','mpfr','mpc']),
#              include_dirs = numpy_include_dirs),

    Extension("psage.modform.maass.test_parallel",
              sources = [ 'psage/modform/maass/test_parallel.pyx' ],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp']) 
]

ext_modules.extend(my_extensions)

# I just had a long chat with Robert Bradshaw (a Cython dev), and he
# told me the following functionality -- turning an Extension with
# Cython code into one without -- along with proper dependency
# checking, is now included in the latest development version of
# Cython (Nov 2, 2010).  It's supposed to be a rewrite he did of the
# code in the Sage library.  Hence once that gets released, we should
# switch to using it here.

build_system.cythonize(ext_modules)

build_system.setup(
    name = 'psage',
    version = "2011.01.06",
    description = "PSAGE: Software for Arithmetic Geometry",
    author = 'William Stein',
    author_email = 'wstein@gmail.com',
    url = 'http://purple.sagemath.org',
    license = 'GPL v2+',
    packages = ['psage',
                'psage.ellcurve',
                'psage.ellcurve.lseries',

#                'psage.ellff',

                'psage.function_fields',

                'psage.lmfdb',
                'psage.lmfdb.ellcurves',
                'psage.lmfdb.ellcurves.sqrt5',

         		'psage.matrix',

                'psage.modform',

                'psage.modform.fourier_expansion_framework',
                'psage.modform.fourier_expansion_framework.gradedexpansions',
                'psage.modform.fourier_expansion_framework.modularforms',
                'psage.modform.fourier_expansion_framework.monoidpowerseries',

                'psage.modform.hilbert',
                'psage.modform.hilbert.sqrt5',

                'psage.modform.rational',

                'psage.modform.siegel',
                'psage.modform.jacobiforms',
                'psage.modform.maass',

        		'psage.modules',

                'psage.number_fields',
                'psage.number_fields.sqrt5',

                'psage.rh',
                'psage.rh.mazur_stein',
                'psage.rings',
                'psage.zfunctions'
                ],
    platforms = ['any'],
    download_url = 'NA',
    ext_modules = ext_modules
)

