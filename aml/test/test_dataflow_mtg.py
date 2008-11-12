from openalea.core.alea import *

pm = PackageManager()
pm.init(verbose=False)

def test_demo_mtg_PlantFrame2():
    """ Test dataflow demo PlantFrame2 """
    res = noserun(('demo.mtg','PlantFrame2'),{},pm=pm)
    assert res == []

def test_demo_mtg_PlantFrame1():
    """ Test dataflow demo PlantFrame1 """
    res = noserun(('demo.mtg','PlantFrame1'),{},pm=pm)
    assert res == []


