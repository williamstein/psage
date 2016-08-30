__all__ = ['all']
from all import *

import logging
import os

LOG_LEVEL = os.getenv('LOG_LEVEL',logging.WARN)
try:
    import colorlog
    from colorlog import ColoredFormatter
    LOGFORMAT = " %(asctime)s %(log_color)s%(levelname)-10s%(filename)s:%(lineno)d%(reset)s | %(log_color)s%(message)s%(reset)s"
    formatter = ColoredFormatter(LOGFORMAT)
    stream = logging.StreamHandler()
    stream.setLevel(LOG_LEVEL)
    stream.setFormatter(formatter)    
except ImportError:
    LOGFORMAT = "  %(asctime)s %(levelname)-10s%(filename)s:%(lineno)d | %(message)s"
    LOG_LEVEL = logging.DEBUG
    #logging.root.setLevel(LOG_LEVEL)
    logging.basicConfig(
        level=LOG_LEVEL,
        format=LOGFORMAT,
#        format='%(asctime)s %(levelname)s: %(message)s '
#        '[in %(pathname)s:%(lineno)d]',
        datefmt='%Y/%m/%d-%H:%M%p',
        )

        
