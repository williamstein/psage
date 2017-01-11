__all__ = ['all']

import logging
import os
logging.basicConfig(
        level=int(os.getenv('LOG_LEVEL',40)),
        format='%(asctime)s %(levelname)s: %(message)s '
        '[in %(pathname)s:%(lineno)d]',
        datefmt='%Y%m%d-%H:%M%p',
        )
