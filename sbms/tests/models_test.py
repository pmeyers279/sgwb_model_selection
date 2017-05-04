import unittest
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from sbms.io_sbms import read_ini
from sbms.models import schumann, power_law, broken_power_law
import numpy.testing as npt
import numpy as np

class TestModel(unittest.TestCase):
    def test_power_law(self):
        omgw_f,f = power_law.omega_gw_spectrum(-9, alpha=0)
        npt.assert_array_almost_equal(omgw_f, 1e-9*np.ones(omgw_f.size))
    def test_broken_power_law(self):
        omgw_f,f = broken_power_law.omega_gw_spectrum(-9, f_break=50,
                alpha1=0, alpha2=0)
        npt.assert_array_almost_equal(omgw_f, 1e-9*np.ones(omgw_f.size))
    def test_schumann_model(self):
        # let's at least run it and see if it fails.
        omgw_f,f = schumann.omega_gw_spectrum()


if __name__=="__main__":
    unittest.main()
