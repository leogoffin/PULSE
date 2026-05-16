import os
import sys
ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

import numpy as np
from modules.findCp import findCp
from modules.atmosphere_isa import Air
from modules.functions import *
from pystdatm import temperature, pressure, density, speed_of_sound, viscosity
from modules.unit_converter import *
import modules.combustion_chambers as cc
import modules.compressors as compr
import modules.turbines as turb
from modules.controllers import Pcontrol