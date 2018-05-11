#################################################################################
#
# (c) Copyright 2010 William Stein
#
#   2016 - 06 - 22: Largely reqritten by Fredrik Stromberg to work with Sage 7.2 and make use of
#                   more standard install routines (based on how Sage setup.py currently works)
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
from sage.env import sage_include_directories,SAGE_INC,SAGE_LIB,SAGE_LOCAL
import subprocess 

if sys.maxint != 2**63 - 1:
    print "*"*70
    print "The PSAGE library only works on 64-bit computers.  Terminating build."
    print "*"*70
    sys.exit(1)

if '-ba' in sys.argv:
    print "Rebuilding all Cython extensions."
    sys.argv.remove('-ba')
    FORCE = True
else:
    FORCE = False

from module_list import ext_modules,aliases

include_dirs = sage_include_directories(use_sources=True)
include_dirs = include_dirs + [SAGE_LIB]
include_dirs = include_dirs + [os.path.join(SAGE_LIB,"cysignals")]
include_dirs = include_dirs + [os.path.join(SAGE_LIB,"sage/ext/")]

extra_compile_args = [ "-fno-strict-aliasing" ]
extra_link_args = [ ]

DEVEL = False
if DEVEL:
    extra_compile_args.append('-ggdb')

# Work around GCC-4.8.0 bug which miscompiles some sig_on() statements,
# as witnessed by a doctest in sage/libs/gap/element.pyx if the
# compiler flag -Og is used. See also
# * http://trac.sagemath.org/sage_trac/ticket/14460
# * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56982
if subprocess.call("""$CC --version | grep -i 'gcc.* 4[.]8' >/dev/null """, shell=True) == 0:
    extra_compile_args.append('-fno-tree-dominator-opts')

    
lib_headers = { "gmp":     [ os.path.join(SAGE_INC, 'gmp.h') ],   # cf. #8664, #9896
                "gmpxx":   [ os.path.join(SAGE_INC, 'gmpxx.h') ],
                "ntl":     [ os.path.join(SAGE_INC, 'NTL', 'config.h') ]
            }
for m in ext_modules:
    m.depends = m.depends + [__file__]

    # Add dependencies for the libraries
    for lib in lib_headers:
        if lib in m.libraries:
            m.depends += lib_headers[lib]

    m.extra_compile_args = m.extra_compile_args + extra_compile_args
    m.extra_link_args = m.extra_link_args + extra_link_args
    m.library_dirs = m.library_dirs + [os.path.join(SAGE_LOCAL, "lib")]
    m.include_dirs = m.include_dirs + include_dirs

#print "include_dirs=",include_dirs

def run_cythonize():
    from Cython.Build import cythonize
    import Cython.Compiler.Options
    import Cython.Compiler.Main
    debug = False
    gdb_debug = True
    if os.environ.get('SAGE_DEBUG', None) == 'yes':
        print('Enabling Cython debugging support')
        debug = True
        Cython.Compiler.Main.default_options['gdb_debug'] = True
        Cython.Compiler.Main.default_options['output_dir'] = 'build'
        gdb_debug=True

    profile = False    
    if os.environ.get('SAGE_PROFILE', None) == 'yes':
        print('Enabling Cython profiling support')
        profile = True
   
    # Sage uses these directives (mostly for historical reasons).
    Cython.Compiler.Options.embed_pos_in_docstring = True
    Cython.Compiler.Options.get_directive_defaults()['autotestdict'] = False
    Cython.Compiler.Options.get_directive_defaults()['cdivision'] = True
    Cython.Compiler.Options.get_directive_defaults()['fast_getattr'] = True
    # The globals() builtin in Cython was fixed to return to the current scope,
    # but Sage relies on the broken behavior of returning to the nearest
    # enclosing Python scope (e.g. to perform variable injection).
    Cython.Compiler.Options.old_style_globals = True
    Cython.Compiler.Main.default_options['cache'] = False
    #print "include_dirs1=",include_dirs
    global ext_modules
    ext_modules = cythonize(
        ext_modules,
        gdb_debug=gdb_debug,
        nthreads=int(os.environ.get('SAGE_NUM_THREADS', 0)),
        #    build_dir=SAGE_CYTHONIZED,
        force=FORCE,
        include_path = include_dirs,
        aliases=aliases,
        compiler_directives={
            'embedsignature': True,
            'profile': profile,
        })
print("Updating Cython code....")
import time
t = time.time()
run_cythonize()
print("Finished Cythonizing, time: %.2f seconds." % (time.time() - t))
import distutils
#for m in ext_modules:
#    print m,isinstance(m,distutils.extension.Extension)
#from distutils.core import setup
from setuptools import setup
code = setup(
    name = 'psage',
    version = "2016.1.0",
    description = "PSAGE: Software for Arithmetic Geometry",
    author = 'William Stein',
    author_email = 'wstein@gmail.com',
    url = 'http://purple.sagemath.org',
    license = 'GPL v2+',
    packages = ['psage',
                'psage.ellcurve',
                'psage.ellcurve.lseries',
                'psage.functions',
#                'psage.ellff',

#                'psage.function_fields',
                'psage.groups',
                
                
                'psage.lmfdb',
                'psage.lmfdb.ellcurves',
                'psage.lmfdb.ellcurves.sqrt5',

         		'psage.matrix',

                'psage.modform',
                'psage.modform.arithgroup',

                'psage.modform.fourier_expansion_framework',
                'psage.modform.fourier_expansion_framework.gradedexpansions',
                'psage.modform.fourier_expansion_framework.modularforms',
                'psage.modform.fourier_expansion_framework.monoidpowerseries',

                'psage.modform.hilbert',
                'psage.modform.hilbert.sqrt5',

                'psage.modform.rational',

                'psage.modform.siegel',
                'psage.modform.jacobi',
                'psage.modform.vector_valued',
                'psage.modform.jacobiforms',
                'psage.modform.weilrep_tools',
                'psage.modform.maass',
                'psage.modform.periods',

        		'psage.modules',

                'psage.number_fields',
                'psage.number_fields.sqrt5',

                'psage.rh',
                'psage.rh.mazur_stein',
                'psage.rings',
                'psage.zfunctions'
                ],
    platforms = ['any'],
    download_url = 'N/A',
    ext_modules = ext_modules
)

