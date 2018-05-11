import logging
import sys


def setup_logging(name, level="INFO"):
    """Create an instance of log and preconfigure it for logging to stdout (useful
    when working within jupyter notebook)

    :param str name: name of the module for logging (usually __name__ should be given)
    :param str level: log level, one of [INFO, DEBUG, CRITICAL, WARN, ERROR, FATAL, NOTSET]
    :returns <logging.Logger> object.
    """

    if level not in logging._nameToLevel.keys():
        raise ValueError('level {} not recognized! Please use one of {}'.format(
            level, logging._nameToLevel.keys()
        ))

    logging.basicConfig(level=level, stream=sys.stdout,
                        format='%(asctime)s %(levelname)-4s : %(message)s',
                        datefmt='%m-%d %H:%M')

    return logging.getLogger(name)
