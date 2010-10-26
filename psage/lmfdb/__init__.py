from mfdb import MFDB

def reload():
    reload = __builtins__['reload']
    import mfdb; reload(mfdb)
    import collection; reload(collection)
    import converter; reload(converter)
    import newforms; reload(newforms)
    import objectdb; reload(objectdb)
