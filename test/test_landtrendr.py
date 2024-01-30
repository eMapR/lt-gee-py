import unittest
import ee
import os
from datetime import date
from ltgee import LandTrendr

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
CHANGE_PARAMS = {
    "delta":  'loss',
    "sort":   'greatest',
    "year":  {"start": LT_PARAMS["start_date"].year + 1, 'end': LT_PARAMS["end_date"].year},
    "mag":    {"value": 200,  "operator": '>'},
    "dur":    {"value": 4,    "operator": '<'},
    "preval": {"value": 300,  "operator": '>'},
    "mmu":    {"value": 5},
}
YEARS = [year for year in
         range(LT_PARAMS['start_date'].year, LT_PARAMS['end_date'].year + 1)]

# Parameter options to test
BASE_INDEX = "NBR"
RECOVERY_THRESHOLDS = [0.25, 0.5, 0.75, 1]
MAX_SEGMENTS = [2, 4, 8, 10]
VALID_INDICES = LandTrendr._valid_indices


def test_coll_list_length(self):
    self.assertEqual(self.py_list.size().getInfo(), len(YEARS))


def test_coll_list_element_diff(self):
    for i in range(len(YEARS)):
        with self.subTest(i=i):
            self.assertTrue(
                ee.Image(self.py_list.get(i)).eq(
                    self.js_list[i]).getInfo())


def test_img_band_names(self):
    self.assertEqual(self.py_img.bandNames().getInfo(),
                     self.js_img.bandNames().getInfo())


def test_img_data_diff(self):
    self.assertTrue(self.py_img.eq(self.js_img).getInfo())


def build_test_class(test_attr_string, debug=False, params=LT_PARAMS, param_name="index", param_value=BASE_INDEX):
    """
    Create a test class for LandTrendr using the above test functions. The test class will be named according to the parameters used.

    Args:
        test_attr_string (str): The LandTrendr class attribute string for the test.
        debug (bool, optional): Debug flag (see LandTrendr.py). Defaults to False.
        params (dict, optional): Parameters for LandTrendr (see LandTrendr.py). Defaults to LT_PARAMS.
        param_name (str, optional): Name of the parameter to adjust. Defaults to "index".
        param_value (str, optional): Value of the parameter to adjust. Must be set with param_name. Defaults to BASE_INDEX.

    Returns:
        class: The test class for LandTrendr.
    """

    if param_name != "index":
        run_params = {param_name: param_value}
        lt = LandTrendr(debug=debug, index=BASE_INDEX,
                        run_params=run_params, **params)
    else:
        lt = LandTrendr(debug=debug, index=param_value, **params)
    param_name = param_name.upper()
    param_value = str(param_value).replace(".", "").upper()

    class LandTrendrTest(unittest.TestCase):
        @classmethod
        def setUpClass(cls):
            match test_attr_string:
                case "SRColl":
                    cls.js_list = [
                        ee.Image(f"{SAMPLES_BASE_PATH}/{year}_SRC") for year in YEARS]
                    cls.py_list = lt.sr_collection.toList(
                        lt.sr_collection.size().getInfo())
                case "LTColl":
                    cls.js_list = [
                        ee.Image(f"{SAMPLES_BASE_PATH}/{year}_LTC_{param_value}") for year in YEARS]
                    cls.py_list = lt.lt_collection.toList(
                        lt.lt_collection.size().getInfo())
                case "LTData":
                    if param_name == "INDEX":
                        cls.js_img = ee.Image(
                            f"{SAMPLES_BASE_PATH}/LTD_{param_value}")
                    else:
                        cls.js_img = ee.Image(
                            f"{SAMPLES_BASE_PATH}/LTD_{BASE_INDEX}_{param_name}-{param_value}")
                    cls.py_img = lt.data

    match test_attr_string:
        case "SRColl" | "LTColl":
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_length", test_coll_list_length)
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_diff", test_coll_list_element_diff)
        case "LTData":
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_band_names", test_img_band_names)
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_diff", test_img_data_diff)
    return LandTrendrTest


if __name__ == '__main__':
    test_classes = []
    test_classes.append(build_test_class(
        "SRColl", debug=True, param_name="index", param_value="NBR"))
    test_classes.append(build_test_class(
        "LTColl", debug=True, param_name="index", param_value="B1"))
    for index in LandTrendr._valid_indices:
        test_classes.append(build_test_class(
            "LTData", param_name="index", param_value=index))
    for recovery_threshold in RECOVERY_THRESHOLDS:
        test_classes.append(build_test_class(
            "LTData", param_name="recoveryThreshold", param_value=recovery_threshold))
    for max_segments in MAX_SEGMENTS:
        test_classes.append(build_test_class(
            "LTData", param_name="maxSegments", param_value=max_segments))
    test_classes.append(build_test_class(
        "LTData", param_name="preventOneYearRecovery", param_value=True))
    test_classes.append(build_test_class(
        "LTData", param_name="preventOneYearRecovery", param_value=False))

    test_loader = unittest.TestLoader()
    test_suite = unittest.TestSuite(
        [test_loader.loadTestsFromTestCase(cls) for cls in test_classes])

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(test_suite)
