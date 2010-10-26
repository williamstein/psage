import os, sys

from setuptools import setup, Extension

SAGE_ROOT = os.environ['SAGE_ROOT']
INCLUDES = ['%s/%s/'%(SAGE_ROOT,x) for x in
            ['local/include/csage', 'local/include', 'local/include/python2.6/',  'devel/sage/sage/ext',
             'devel/sage', 'devel/sage/sage/gsl']]

def cython(f, m):
    """
    Given a pair p = (f, m), with a .pyx file f which is a part the
    module m, call Cython on f

    INPUT:
         p -- a 2-tuple f, m

    copy the file to SITE_PACKAGES, and return a string
    which will call Cython on it.
    """
    assert f.endswith('.pyx')
    # output filename
    dir, f = os.path.split(f)
    ext = 'cpp' if m.language == 'c++' else 'c'
    outfile = os.path.splitext(f)[0] + '.' + ext
    includes = ''.join(["-I '%s' "%x for x in INCLUDES])

    # call cython
    cmd = "cd %s && python `which cython` --embed-positions --directive cdivision=True %s -o %s %s"%(
        dir, includes, outfile, f)
                       
    print cmd
    return os.system(cmd)

E = Extension("psage.ellff.ellff",
              ["psage/ellff/ellff.cpp"],
              language='c++',
              include_dirs = INCLUDES
              )

cython('psage/ellff/ellff.pyx', E)

ext_modules = [
    E
]



setup(
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

