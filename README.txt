This is the PSAGE library:

    http://code.google.com/p/purplesage/

The target audience of PSAGE is research mathematicians, with an
emphasis on arithmetic geometry. PSAGE is closely related to Sage, but
is aimed squarely at research mathematicians.

UNIT TESTS:

  Install nose: 
    sage -sh; easy_install nose

  Then from the root of the psage install:
    nosetests -v --processes=8   # run 8 tests in parallel
  
DOCTESTS:

  Do this to run 8 tests in parallel:

    sage -tp 8 --force_lib psage/

  A lot of doctests in certain parts of psage will fail.   This is because various authors wrote doctests for code
  as if it were in Sage (or as if the file is pre-imported), but for the tests to pass in Sage one must explicitly
  import functions from psage. 
 
