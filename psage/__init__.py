# Import Python code
import modform

# Import Cython code
try:
    from function_fields import FunctionField
    from ellff import ellff_EllipticCurve
except ImportError, msg:
    print msg
    print "WARNING: Error importing Cython modules."
    pass
