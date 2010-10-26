import function_field

def FunctionField(X, names=None):
    """
    Return the function field defined by X.

    INPUT:

        - `X` -- a field; return the function field in one variable over X.

        - ``names`` -- name of variable as a string
    
    EXAMPLES::

        sage: FunctionField(QQ,'alpha')
        Rational function field in alpha over Rational Field
        sage: K.<alpha> = FunctionField(GF(7)); K
        Rational function field in alpha over Finite Field of size 7
    """
    return function_field.RationalFunctionField(X, names=names)
    
