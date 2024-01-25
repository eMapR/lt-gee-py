# lt-gee-py
Python interface to the Google Earth Engine implementation of the LandTrendr spectral-temporal segmentation algorithm.

## Introduction

LandTrendr is set of spectral-temporal segmentation algorithms that are useful for change detection in a time series of moderate resolution satellite imagery (primarily Landsat) and for generating trajectory-based spectral time series data largely absent of inter-annual signal noise. LT was originally implemented in IDL (Interactive Data Language), but with the help of engineers at Google, it has been ported to the GEE platform.

The LandTrendr class is a light wrapper around the Google Earth Engine API to that includes convenience methods to generate images in the format required for the algorithm on GEE. 

## Getting Started

### Download and Install packages and dependencies

- Install the Python API for Google Earth Engine. This is the only dependency thus far.

```
conda install -c conda-forge earthengine-api
```

- Install package

```
pip install lt-gee-py
```

## Basic Usage

```python
import ee
from lt-gee-py import LandTrendr
# Initialize access to Google's EE servers
project_name = "my_project_name
ee.Initialize(project_name)

# Initialize variables for LandTrendr algorithm
lt_params = {
    "start_date": date(1985, 6,1),
    "end_date": date(2017, 9,1),
    "index": 'NBR',
    "ftv_list": ['TCB', 'TCG', 'TCW', 'NBR'],
    "mask_labels": ['cloud', 'shadow', 'snow', 'water'],

    "aoi": ee.Geometry.Point(-122.8848, 43.7929),
    "run_params": {
            "maxSegments": 6,
            "spikeThreshold": 0.9,
            "vertexCountOvershoot":  3,
            "preventOneYearRecovery":  True,
            "recoveryThreshold":  0.25,
            "pvalThreshold":  .05,
            "bestModelProportion":  0.75,
            "minObservationsNeeded": 6,
        }
}

# Instantiating LandTrendr object. Note: The object will immediately request to run the algorithm on Google's servers.
lt = LandTrendr(**lt_params)

# Access resulting image using the data attribute.
lt_data = lt.data
```

## [Manuscript](http://www.mdpi.com/2072-4292/10/5/691) 

## Citation

>Kennedy, R.E., Yang, Z., Gorelick, N., Braaten, J., Cavalcante, L., Cohen, W.B., Healey, S. (2018). Implementation of the LandTrendr Algorithm on Google Earth Engine. Remote Sensing. 10, 691.

Except as otherwise noted, the content of this repository and accompanying description site (https://emapr.github.io/LT-GEE/) are licensed under the [Creative Commons Attribution 4.0 License](https://creativecommons.org/licenses/by/4.0/), and code samples are licensed under the [Apache 2.0 License](https://www.apache.org/licenses/LICENSE-2.0).