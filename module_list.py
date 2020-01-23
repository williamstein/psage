r"""
List of extension modules
"""
#from distutils.extension import Extension
import os
import pkgconfig

# CBLAS can be one of multiple implementations
cblas_pc = pkgconfig.parse('cblas')
cblas_libs = list(cblas_pc['libraries'])
cblas_library_dirs = list(cblas_pc['library_dirs'])
cblas_include_dirs = list(cblas_pc['include_dirs'])
# TODO: Remove Cygwin hack by installing a suitable cblas.pc
if os.path.exists('/usr/lib/libblas.dll.a'):
    cblas_libs = ['gslcblas']
# GNU Scientific Library
# Note we replace the built-in gslcblas with the above cblas
gsl_pc = pkgconfig.parse('gsl')
gsl_libs = list(set(gsl_pc['libraries']).difference(set(['gslcblas'])).union(set(cblas_libs)))
gsl_library_dirs = list(gsl_pc['library_dirs'])
gsl_include_dirs = list(gsl_pc['include_dirs'])

aliases = dict(
    GSL_LIBRARIES=gsl_libs,
    GSL_LIBDIR=gsl_library_dirs,
    GSL_INCDIR=gsl_include_dirs,
)
import numpy
numpy_include_dirs = [numpy.get_include()]
import setuptools
class Extension(setuptools.extension.Extension):
    def __init__(self, name, sources, include_dirs=[],
                  language="c", force=False, **kwds):
        #print "kwds=",kwds
        #print "module=",module
        setuptools.Extension.__init__(self, name, sources, language=language,
                                       include_dirs=include_dirs, **kwds)         
#     def __init__(self, module, sources, include_dirs,
#                  language="c", force=False, **kwds):
#         self.cython_cmds = []
#         for i in range(len(sources)):
#             f = sources[i]
#             if f.endswith('.pyx'):
#                 sources[i], cmds = cython(f) #, language, include_dirs, force)
#                 for c in cmds:
#                     self.cython_cmds.append(c)
#         setuptools.Extension.__init__(self, module, sources, language=language,
#                                       include_dirs=include_dirs, **kwds)


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

#    Extension("psage.function_fields.function_field_element",
#              ["psage/function_fields/function_field_element.pyx"]),

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

    Extension('psage.modform.arithgroup.mysubgroups_alg',
              ['psage/modform/arithgroup/mysubgroups_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),

    Extension('psage.modform.maass.maass_forms_alg',
              ['psage/modform/maass/maass_forms_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc','ntl'],
              include_dirs = numpy_include_dirs),

    Extension('psage.modform.maass.lpkbessel',
              ['psage/modform/maass/lpkbessel.pyx'],
              libraries = ['m', 'gmp','mpfr','mpc','ntl'],
              include_dirs = numpy_include_dirs),

              
    Extension("psage.modform.rational.modular_symbol_map",
              ["psage/modform/rational/modular_symbol_map.pyx"]),

    Extension("psage.modform.rational.padic_elliptic_lseries_fast",
              ["psage/modform/rational/padic_elliptic_lseries_fast.pyx"]),

#     Extension("psage.modform.hilbert.sqrt5.sqrt5_fast",
#               ["psage/modform/hilbert/sqrt5/sqrt5_fast.pyx"],
# #              libraries = ['ntl', 'gmp'],
#               libraries = ['gmp']),
# #              language = 'c++'),

    # Extension("psage.ellcurve.lseries.sqrt5",
    #           ["psage/ellcurve/lseries/sqrt5.pyx"],
    #           libraries = ['ntl', 'gmp'],
    #           language = 'c++'),

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
#              ["psage/modform/rational/special_fast.pyx", SAGE_ROOT + "/devel/sage/sage/libs/flint/fmpq_poly.c"],
              ["psage/modform/rational/special_fast.pyx"],
              libraries = ['gmp', 'flint'],
              language = 'c++'),
#              include_dirs = [SAGE_LOCAL + '/include/FLINT/', SAGE_ROOT + '/devel/sage/sage/libs/flint/'],
#              extra_compile_args = ['-std=c99']),

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
# Compile using openmp and hide deprecatiopn warning for numpy API since we don't use it here
extra_compile_args = ['-fopenmp','-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION']
## Fredrik Stroemberg: my additional modules.
my_extensions = [
    Extension('psage.groups.permutation_alg',
              ['psage/groups/permutation_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc']),
    Extension('psage.rings.mp_cimports',
              sources= ['psage/rings/mp_cimports.pyx'],
              libraries = ['gmp','mpfr','mpc']),
    Extension('psage.rings.mpc_extras',
              sources = ['psage/rings/mpc_extras.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=extra_compile_args,
              extra_link_args=['-fopenmp']),

    Extension('psage.modform.periods.period_polynomials_algs',
              ['psage/modform/periods/period_polynomials_algs.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=extra_compile_args,
              extra_link_args=['-fopenmp']),

    Extension('psage.functions.inc_gamma',
              ['psage/functions/inc_gamma.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=extra_compile_args,
              extra_link_args=['-fopenmp']),
    
    Extension('psage.modform.arithgroup.mysubgroups_alg',
              ['psage/modform/arithgroup/mysubgroups_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc']),

    Extension('psage.modform.arithgroup.sl2z_subgroups_alg',
              ['psage/modform/arithgroup/sl2z_subgroups_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc']),

    Extension('psage.modform.maass.maass_forms_parallel_alg',
              ['psage/modform/maass/maass_forms_parallel_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=extra_compile_args,
              extra_link_args=['-fopenmp']),

    Extension('psage.modform.maass.pullback_algorithms',
              ['psage/modform/maass/pullback_algorithms.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=extra_compile_args,
              include_dirs = numpy_include_dirs,
              extra_link_args=['-fopenmp']),

    Extension('psage.modform.maass.linear_system',
              ['psage/modform/maass/linear_system.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),
    Extension('psage.modform.maass.automorphic_forms_alg',
              ['psage/modform/maass/automorphic_forms_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs,
              extra_compile_args=extra_compile_args,
              extra_link_args=['-fopenmp']),

    Extension('psage.modform.hilbert.hn_class',
              ['psage/modform/hilbert/hn_class.pyx'],
              libraries = ['m']),
    
    Extension('psage.modform.hilbert.hilbert_modular_group_alg',
              ['psage/modform/hilbert/hilbert_modular_group_alg.pyx'],
              libraries = ["flint", "gmp", "gmpxx", "m","ntl"],
              language="c"),

    Extension('psage.zfunctions.selberg_z_alg',
              ['psage/zfunctions/selberg_z_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs,
              extra_compile_args=extra_compile_args,
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
              extra_compile_args=extra_compile_args,
              extra_link_args=['-fopenmp'],
              include_dirs = numpy_include_dirs),

    Extension('psage.matrix.linalg_complex_dense',
              sources = ['psage/matrix/linalg_complex_dense.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=extra_compile_args,
              extra_link_args=['-fopenmp']),
    Extension('psage.modform.maass.poincare_series_alg',
                  ['psage/modform/maass/poincare_series_alg.pyx'],
              libraries = ['m','gmp','mpfr','mpc']),
    
    Extension('psage.modform.maass.poincare_series_alg_vv',
          ['psage/modform/maass/poincare_series_alg_vv.pyx'],
          libraries = ['m','gmp','mpfr','mpc']),
        
    Extension('psage.modform.maass.eisenstein_series',
              ['psage/modform/maass/eisenstein_series.pyx'],
              libraries = ['m','gmp','mpfr','mpc'],
              include_dirs = numpy_include_dirs),

    Extension("psage.modform.maass.test_parallel",
              sources = [ 'psage/modform/maass/test_parallel.pyx' ],
              libraries = ['m','gmp','mpfr','mpc'],
              extra_compile_args=extra_compile_args,
              extra_link_args=['-fopenmp']),


    Extension("psage.groups.dirichlet_conrey",
              ['psage/groups/dirichlet_conrey.pyx'],
              extra_compile_args = ['-w','-O2'])
]

ext_modules.extend(my_extensions)

## Stephan Ehlen's additional modules

sehlen_extensions = [
      Extension('psage.modules.weil_invariants',
              sources = ['psage/modules/weil_invariants.pyx'],
              libraries = ['m']
#              include_dirs = numpy_include_dirs),
     )
]

ext_modules.extend(sehlen_extensions)
