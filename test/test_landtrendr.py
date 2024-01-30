import unittest
import ee
import os
from datetime import date
from ltgee import LandTrendr

# TODO: create generic assets path for testing
ee.Initialize(project=os.getenv('GEE_PROJECT'))
SAMPLES_BASE_PATH = os.getenv('GEE_LT_PY_SAMPLES_BASE_PATH')
TEST_AOI = ee.Geometry.Point([-122.11499066457009, 44.477400937191916])
LT_PARAMS = {
    "start_date": date(1985, 6, 1),
    "end_date": date(2017, 9, 1),
    "ftv_list": ['TCB', 'TCG', 'TCW',],
    "mask_labels": ['cloud', 'shadow', 'snow', 'water'],
    "aoi": TEST_AOI
}
BASE_INDEX = "NBR"
RECOVERY_THRESHOLDS = [0.25, 0.5, 0.75, 1]
MAX_SEGMENTS = [2, 4, 8, 10]
FTV_LIST = ['TCB', 'TCG', 'TCW']
CHANGE_PARAMS = {
    "delta":  'loss',
    "sort":   'greatest',
    "year":  {"start": LT_PARAMS["start_date"].year + 1, 'end': LT_PARAMS["end_date"].year},
    "mag":    {"value": 200,  "operator": '>'},
    "dur":    {"value": 4,    "operator": '<'},
    "preval": {"value": 300,  "operator": '>'},
    "mmu":    {"value": 5},
}
YEARS = [year for year in range(
    LT_PARAMS['start_date'].year, LT_PARAMS['end_date'].year + 1)]

# TODO: Refactor to better parameterize test cases to avoid redefining runTest() for each test case
# The excessive inheretence causes issues with unittest discovery
class LandTrendrBaseTest(unittest.TestCase):
    def __init__(self, debug, params, param_name, param_value, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if param_name != "index":
            run_params = {param_name: param_value}
            self.lt = LandTrendr(debug=debug, index=BASE_INDEX,
                                 run_params=run_params, **params)
        else:
            self.lt = LandTrendr(debug=debug, index=param_value, **params)
        self.param_name = param_name.upper()
        self.param_value = str(param_value).replace(".", "").upper()


class LandTrendrSRColTest(LandTrendrBaseTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sr_col_js_list = [
            ee.Image(f"{SAMPLES_BASE_PATH}/{year}_SRC") for year in YEARS]
        self.sr_col_py_list = self.lt.sr_collection.toList(
            self.lt.sr_collection.size().getInfo())

    def col_length(self):
        self.assertEqual(self.sr_col_py_list.size().getInfo(), len(YEARS))

    def col_element_diff(self):
        for i in range(len(YEARS)):
            self.assertTrue(ee.Image(self.sr_col_py_list.get(i)).eq(
                self.sr_col_js_list[i]).getInfo())

    def runTest(self):
        self.col_length()
        self.col_element_diff()


class LandTrendrLTColTest(LandTrendrBaseTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.lt_col_js_list = [ee.Image(
            f"{SAMPLES_BASE_PATH}/{year}_LTC_{self.param_value}") for year in YEARS]
        self.lt_col_py_list = self.lt.lt_collection.toList(
            self.lt.lt_collection.size().getInfo())

    def col_length(self):
        self.assertEqual(self.lt_col_py_list.size().getInfo(), len(YEARS))

    def col_element_diff(self):
        for i in range(len(YEARS)):
            self.assertTrue(ee.Image(self.lt_col_py_list.get(i)).eq(
                self.lt_col_js_list[i]).getInfo())

    def runTest(self):
        self.col_length()
        self.col_element_diff()


class LandTrendrLTDataTest(LandTrendrBaseTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.param_name == "INDEX":
            self.lt_data_js = ee.Image(
                f"{SAMPLES_BASE_PATH}/LTD_{self.param_value}")
        else:
            self.lt_data_js = ee.Image(
                f"{SAMPLES_BASE_PATH}/LTD_{BASE_INDEX}_{self.param_name}-{self.param_value}")

    def lt_band_names(self):
        self.assertEqual(self.lt.data.bandNames().getInfo(),
                         self.lt_data_js.bandNames().getInfo())

    def lt_data_diff(self):
        self.assertTrue(self.lt.data.eq(self.lt_data_js).getInfo())

    def runTest(self):
        self.lt_band_names()
        self.lt_data_diff()


if __name__ == '__main__':
    test_suite = unittest.TestSuite()
    test_suite.addTest(LandTrendrSRColTest(
        debug=True, params=LT_PARAMS, param_name="index", param_value="NBR"))
    test_suite.addTest(LandTrendrLTColTest(
        debug=True, params=LT_PARAMS, param_name="index", param_value="B1"))
    for index in LandTrendr._valid_indices:
        test_suite.addTest(LandTrendrLTDataTest(
            debug=False, params=LT_PARAMS, param_name="index", param_value=index))
    for recovery_threshold in RECOVERY_THRESHOLDS:
        test_suite.addTest(LandTrendrLTDataTest(
            debug=False, params=LT_PARAMS, param_name="recoveryThreshold", param_value=recovery_threshold))
    for max_segments in MAX_SEGMENTS:
        test_suite.addTest(LandTrendrLTDataTest(
            debug=False, params=LT_PARAMS, param_name="maxSegments", param_value=max_segments))
    test_suite.addTest(LandTrendrLTDataTest(
        debug=False, params=LT_PARAMS, param_name="preventOneYearRecovery", param_value=True))
    test_suite.addTest(LandTrendrLTDataTest(
        debug=False, params=LT_PARAMS, param_name="preventOneYearRecovery", param_value=False))
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(test_suite)
