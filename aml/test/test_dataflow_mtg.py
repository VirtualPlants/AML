from openalea.core.alea import *

pm = PackageManager()
pm.init(verbose=False)

def test_demo_mtg_PlantFrame2():
    """ Test dataflow demo PlantFrame2 """
    res = run(('demo.mtg','PlantFrame2'),{},pm=pm)
    assert res == []

def test_demo_mtg_multiplemtg():
    """ Test dataflow demo PlantFrame1 """
    res = run(('demo.mtg','MTGActionOnMultipleMTGs'),{},pm=pm)
    assert res == []

def test_demo_mtg_selection():
    """ Test dataflow demo PlantFrame1 """
    res = run(('demo.mtg','MTGComponentSelection'),{},pm=pm)

    assert res == []

def test_demo_mtg_extraction():
    """ Test dataflow demo PlantFrame1 """
    res = run(('demo.mtg','MTGExtraction'),{},pm=pm)
    assert res == []

