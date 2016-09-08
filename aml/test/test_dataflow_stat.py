from openalea.core.alea import *

pm = PackageManager()
pm.init(verbose=False)

# These tests use gnuplot interface, which requires human interaction
# Consequently, they cannot be used within builbot (which hangs forever)
# We added a flags inside aml/src/aml/wralea/py_stat.py to prevent gnuplot
# to be launched if these tests are run with nosetests. The remaining of the
# nodes are run.
# In order to have the gnuplot interface, run this script with python instead of nosetests 
def test_demo_corsican():
    """ Test changepoint demo corsican  """
    res = run(('demo.changepoint','Corsican pine change point'),{},pm=pm)
    assert res == []

def test_demo_mtg_multiplemtg():
    """ Test dataflow demo dycorinia """
    res = run(('demo.changepoint','Dycorinia change point'),{},pm=pm)
    assert res == []


if __name__ == "__main__":
    test_demo_corsican()
    test_demo_mtg_multiplemtg()

