import unittest
from ..io import read_ini
from ..models import *

class TestModel(unittest.TestCase):
    def test_power_law(self):
        params = read_ini('sbms/tests/data/test_ini.ini')
        print params
        omgw_f,f = power_law.omega_gw_spectrum(1e-9, alpha=0)
        omgw_f2,f = omgwf(params)
        print omgw_f
        print omgw_f2
if __name__=="__main__":
    unittest.main()
