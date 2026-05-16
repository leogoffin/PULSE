import os
import sys
ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

import numpy as np
from findCp import findCp
from atmosphere_isa import Air
from functions import *
from pystdatm import temperature, pressure, density, speed_of_sound, viscosity
from unit_converter import *
import combustion_chambers as cc
import compressors as compr
import turbines as turb
from controllers import Pcontrol