from __future__ import absolute_import
# import modules 
from . import automorphic_forms
from . import eisenstein_series 
from . import maass_forms
from . import multiplier_systems
from . import poincare_series
from . import maass_forms_alg
from . import maass_forms_phase2 
from . import find_maass_forms

# import specific functionality that we think the end user wants

# spaces 
from .automorphic_forms import (AutomorphicFormSpace,
                               HalfIntegralWeightForms,
                               HarmonicWeakMaassForm,
                               HarmonicWeakMaassFormSpace,
                               HolomorphicModularForms)


from .maass_forms import MaassWaveForms, EisensteinSeries, Maasswaveform
from .weil_rep_simple import WeilRepDiscriminantForm
from .multiplier_systems import (MultiplierSystem,
                                TrivialMultiplier,
                                ThetaMultiplier,
                                EtaMultiplier,
                                TestMultiplier,
                                MultiplierByGenerator,
                                InducedRepresentationMultiplier,
                                WeilRepMultiplier,EtaQuotientMultiplier)
from .poincare_series import PoincareSeries
from .vv_harmonic_weak_maass_forms import VVHarmonicWeakMaassForms
#from maass_forms_alg import get_Y_from_M,get_M_from_Y,get_M_and_Y
#from lpkbessel import besselk_dp
#from pullback_algorithms import *
#from maass_forms_alg import *
#from maass_forms_phase2 import *
#from automorphic_forms_alg import *
#from linear_system import *
#from poincare_series_vv import *
#from harmonic_test import *
