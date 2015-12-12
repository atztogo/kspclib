#!/usr/local/bin/python
import sys
import logging
import logging.config
import os
import yaml

def config_logger(filename="/logging.yaml",level=logging.INFO):
        """
        Load logging configuration from file, if file is not found 
        set the logging level to DEBUG
        
        """
        path_to_script=os.path.dirname(os.path.realpath(__file__))
        filename=path_to_script+filename
        if os.path.exists(filename):
                print "Setting logger from logging.yaml"
                sys.stdout.flush()
                with open(filename,'rt') as f:
                        config=yaml.load(f.read())
                logging.config.dictConfig(config)
        else:
                print "Setting default logger"
                sys.stdout.flush()
                logging.basicConfig(level=level)
