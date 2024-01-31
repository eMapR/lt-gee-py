import unittest
import ee
import os
from datetime import date
from ltgee import LandTrendr

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


def compare_images(img1, img2):
    bands = img1.bandNames().getInfo()
    if "LandTrendr" in bands:
        # TODO: implement parsing and testing for landtrendr Data arrays
        return False
    res = img1.subtract(img2).divide(img2).abs().reduceRegion(
        reducer=ee.Reducer.mean(), geometry=TEST_AOI).getInfo()
    # TODO: There is potentially a floating point error to address however the difference is always small
    return [True if _ < 1e-7 else False for _ in res.values()]


def test_sr_list_length(self):
    self.assertEqual(self.py_sr_list.size().getInfo(), len(YEARS))


def test_lt_list_length(self):
    self.assertEqual(self.py_lt_list.size().getInfo(), len(YEARS))


def test_sr_list_element_diff(self):
    for i in range(len(YEARS)):
        with self.subTest(i=i):
            res = compare_images(
                ee.Image(self.py_sr_list.get(i)), self.js_sr_list[i])
            self.assertTrue(all(res))


def test_lt_list_element_diff(self):
    for i in range(len(YEARS)):
        with self.subTest(i=i):
            res = compare_images(
                ee.Image(self.py_lt_list.get(i)), self.js_lt_list[i])
            self.assertTrue(all(res))


def test_img_band_names(self):
    self.assertEqual(self.py_img.bandNames().getInfo(),
                     self.js_img.bandNames().getInfo())


def test_img_data_diff(self):
    res = compare_images(self.py_img, self.js_img)
    self.assertTrue(all(res))


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
                case "CollBuilders":
                    cls.js_sr_list = [
                        ee.Image(f"{SAMPLES_BASE_PATH}/{year}_SRC") for year in YEARS]
                    cls.py_sr_list = lt.sr_collection.toList(
                        lt.sr_collection.size().getInfo())
                    cls.js_lt_list = [
                        ee.Image(f"{SAMPLES_BASE_PATH}/{year}_LTC_{param_name}_{param_value}") for year in YEARS]
                    cls.py_lt_list = lt.lt_collection.toList(
                        lt.lt_collection.size().getInfo())
                case "LTData":
                    if param_name == "INDEX":
                        cls.js_img = ee.Image(
                            f"{SAMPLES_BASE_PATH}/LTD_{param_name}_{param_value}")
                    else:
                        cls.js_img = ee.Image(
                            f"{SAMPLES_BASE_PATH}/LTD_INDEX_{BASE_INDEX}_{param_name}_{param_value}")
                    cls.py_img = lt.data

    match test_attr_string:
        case "CollBuilders":
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_sr_length", test_sr_list_length)
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_lt_length", test_lt_list_length)
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_sr_diff", test_sr_list_element_diff)
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_lt_diff", test_lt_list_element_diff)
        case "LTData":
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_band_names", test_img_band_names)
            setattr(
                LandTrendrTest, f"test_{test_attr_string}_{param_name}_{param_value}_diff", test_img_data_diff)
    return LandTrendrTest


if __name__ == '__main__':
    ee.Initialize()
    test_classes = []
    test_classes.append(build_test_class("CollBuilders", debug=True))
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
