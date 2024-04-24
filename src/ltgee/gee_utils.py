import ee
import math


def count_clear_view_pixels(collection: ee.ImageCollection) -> ee.Image:
    """
    Counts the number of clear view pixels for the given collection.

    Args:
        collection (ee.ImageCollection): The image collection.

    Returns:
        ee.Image: The image containing the number of clear view pixels.
    """
    binary = collection.map(lambda image: image.select(
        0).multiply(0).add(1).unmask(0))
    return binary.sum()


def calculate_median_diff(img: ee.Image, median: ee.Image) -> ee.Image:
    """
    Calculates the difference between an image and a given median value.

    Parameters:
        img (ee.Image): The input image.
        median (ee.Image): The median image.

    Returns:
        ee.Image: The difference image with an additional band representing the input image.
    """
    diff = img.subtract(median).pow(ee.Image.constant(2))
    return diff.reduce('sum').addBands(img)


def forest_mask(aoi: ee.Geometry) -> ee.Image:
    """
    Generate a forest mask for a given area of interest (AOI).

    Parameters:
    aoi (ee.Geometry): The area of interest for which the forest mask is generated.

    Returns:
    ee.Image: The forest mask image clipped to the AOI.
    """
    for_col = ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V/Global")
    img_for = for_col.toBands()
    forest_image = img_for.select('2015_forest_type')
    selected_forests = forest_image.expression(
        'Band >= 0 ? 1 : 0', {
            'Band': forest_image
        }
    )
    return selected_forests.clip(aoi)


def water_mask(aoi: ee.Geometry) -> ee.Image:
    """
    Generate a water mask to the given area of interest (aoi).

    Parameters:
    aoi (ee.Geometry): The area of interest to apply the water mask to.

    Returns:
    ee.Image: The water mask applied to the given area of interest.
    """
    mapped_water = ee.Image("JRC/GSW1_1/GlobalSurfaceWater")
    mapped_water_binary = mapped_water.expression(
        'band > 99 ? 0 : 1', {
            'band': mapped_water.select('recurrence')
        }
    )
    return mapped_water_binary.clip(aoi)


def standardize_collection(collection: ee.ImageCollection) -> ee.ImageCollection:
    """
    Standardizes the given collection by subtracting the mean and dividing by the standard deviation.

    Args:
        collection (ee.ImageCollection): The image collection to be standardized.

    Returns:
        ee.ImageCollection: The standardized image collection.
    """
    mean = collection.reduce(ee.Reducer.mean())
    std_dev = collection.reduce(ee.Reducer.stdDev())
    mean_adj = collection.map(lambda img: img.subtract(mean)
                              .set('system:time_start', img.get('system:time_start')))
    return mean_adj.map(lambda img: img.divide(std_dev).set('system:time_start', img.get('system:time_start')))


def tc_transform(img: ee.Image,
                 color_bands: list[str] = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']) -> ee.Image:
    b = ee.Image(img).select(color_bands)

    brt_coeffs = ee.Image.constant(
        [0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303])
    grn_coeffs = ee.Image.constant(
        [-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446])
    wet_coeffs = ee.Image.constant(
        [0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109])
    sum_reducer = ee.Reducer.sum()

    brightness = b.multiply(brt_coeffs).reduce(sum_reducer)
    greenness = b.multiply(grn_coeffs).reduce(sum_reducer)
    wetness = b.multiply(wet_coeffs).reduce(sum_reducer)
    angle = (greenness.divide(brightness)).atan().multiply(
        180/math.pi).multiply(100)

    return brightness.addBands(greenness) \
        .addBands(wetness) \
        .addBands(angle) \
        .select([0, 1, 2, 3], ['TCB', 'TCG', 'TCW', 'TCA'])\
        .set('system:time_start', img.get('system:time_start'))


def nbr_transform(img: ee.Image) -> ee.Image:
    return img.normalizedDifference(['B4', 'B7']) \
        .multiply(1000) \
        .select([0], ['NBR']) \
        .set('system:time_start', img.get('system:time_start'))


def ndfi_transform(img: ee.Image) -> ee.Image:
    params = {
        'cfThreshold': 0.01,
        'soil': [2000, 3000, 3400, 5800, 6000, 5800],
        'gv': [500, 900, 400, 6100, 3000, 1000],
        'npv': [1400, 1700, 2200, 3000, 5500, 3000],
        'shade': [0, 0, 0, 0, 0, 0],
        'cloud': [9000, 9600, 8000, 7800, 7200, 6500]
    }
    gv = params['gv']
    shade = params['shade']
    npv = params['npv']
    soil = params['soil']
    cloud = params['cloud']

    unmix_image = ee.Image(img).unmix([gv, shade, npv, soil, cloud], True, True).rename(
        ['band_0', 'band_1', 'band_2', 'band_3', 'band_4'])
    new_image = ee.Image(img).addBands(unmix_image)

    ndfi = unmix_image.expression(
        '((GV / (1 - SHADE)) - (NPV + SOIL)) / ((GV / (1 - SHADE)) + NPV + SOIL)', {
            'GV': unmix_image.select('band_0'),
            'SHADE': unmix_image.select('band_1'),
            'NPV': unmix_image.select('band_2'),
            'SOIL': unmix_image.select('band_3')
        })

    ndvi = ee.Image(img).normalizedDifference(['B4', 'B3']).rename('NDVI')

    evi = ee.Image(img).expression(
        'float(2.5*(((B4/10000) - (B3/10000)) / ((B4/10000) + (6 * (B3/10000)) - (7.5 * (B1/10000)) + 1)))',
        {
            'B4': img.select(['B4']),
            'B3': img.select(['B3']),
            'B1': img.select(['B1'])
        }).rename('EVI')

    to_exp = new_image.addBands([ndfi.rename(['NDFI']), ndvi, evi]) \
        .select(['band_0', 'band_1', 'band_2', 'band_3', 'NDFI', 'NDVI', 'EVI', 'B1', 'B2', 'B3', 'B4', 'B5']) \
        .rename(['GV', 'Shade', 'NPV', 'Soil', 'NDFI', 'NDVI', 'EVI', 'Blue', 'Green', 'Red', 'NIR', 'SWIR1'])

    return to_exp.select(['NDFI']).multiply(1000).set(
        'system:time_start', img.get('system:time_start'))


def ndvi_transform(img: ee.Image,
                   nir_band_name: str = 'B4',
                   red_band_name: str = 'B5') -> ee.Image:
    ndvi = img.normalizedDifference([nir_band_name, red_band_name]) \
        .multiply(1000) \
        .select([0], ['NDVI']) \
        .set('system:time_start', img.get('system:time_start'))
    return ndvi


def ndsi_transform(img: ee.Image,
                   swir_band_name: str = 'B5',
                   green_band_name: str = 'B2') -> ee.Image:
    ndsi = img.normalizedDifference([green_band_name, swir_band_name]) \
        .multiply(1000) \
        .select([0], ['NDSI']) \
        .set('system:time_start', img.get('system:time_start'))
    return ndsi


def ndmi_transform(img: ee.Image,
                   nir_band_name: str = 'B4',
                   swir_band_name: str = 'B5') -> ee.Image:
    ndmi = img.normalizedDifference([nir_band_name, swir_band_name]) \
        .multiply(1000) \
        .select([0], ['NDMI']) \
        .set('system:time_start', img.get('system:time_start'))
    return ndmi


def evi_transform(img: ee.Image,
                  nir_band_name: str = 'B4',
                  red_band_name: str = 'B3',
                  blue_band_name: str = 'B1') -> ee.Image:
    evi = img.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
            'NIR': img.select(nir_band_name),
            'RED': img.select(red_band_name),
            'BLUE': img.select(blue_band_name)
        }) \
        .multiply(1000) \
        .select([0], ['EVI']) \
        .set('system:time_start', img.get('system:time_start'))
    return evi
