import ee

from datetime import datetime
from typing import Optional

from .gee_utils import (count_clear_view_pixels,
                        water_mask,
                        forest_mask,
                        calculate_median_diff,
                        tc_transform,
                        ndvi_transform,
                        ndmi_transform,
                        ndsi_transform,
                        nbr_transform,
                        evi_transform,
                        ndfi_transform,
                        standardize_collection)


class Sentinel2Composite(ee.ImageCollection):
    """
    SentinelComposite class for building a cloud-masked Sentinel-2 surface reflectance composite.

    Args:
        area_of_interest (ee.Geometry): Area of interest.
        start_date (datetime): Start date for the image collection.
        end_date (datetime): End date for the image collection.
        cloud_filter (int): Max cloud cover percentage allowed in images.
        cloud_prb_thresh (int): Cloud probability threshold for masking.
        nir_drk_thresh (float): NIR reflectance threshold for cloud shadows.
        cloud_prj_dist (float): Max distance (in km) to project cloud shadows.
        buffer (int): Buffer distance (in meters) for cloud edge.
        exclude (dict): Dictionary for exluding images in the compositing
    """
    _sr_band_scale = 1e4

    def __init__(self,
                 area_of_interest: ee.Geometry,
                 start_date: datetime,
                 end_date: datetime,
                 cloud_filter: int = 20,
                 cloud_prb_thresh: int = 65,
                 nir_drk_thresh: float = 0.2,
                 cloud_prj_dist: float = 1.5,
                 buffer: int = 50,
                 exclude: dict = {}):
        self.area_of_interest = area_of_interest
        self.start_date = start_date
        self.end_date = end_date
        self.cloud_filter = cloud_filter
        self.cloud_prb_thresh = cloud_prb_thresh
        self.nir_drk_thresh = nir_drk_thresh
        self.cloud_prj_dist = cloud_prj_dist
        self.buffer = buffer
        self.exclude = exclude
        super().__init__(self._build_sr_collection())
        
    # TODO: assertion checks on thresholds
    def _build_sr_collection(self, debug: Optional[bool] = False):
        """
        Builds a medoid composite of Landsat surface reflectance TM-equivalent bands 1,2,3,4,5,7 from Sentinel 2 Harmonized
        This collection can be useful outside of use by LandTrendr, but is also the base for creating the input collection for LandTrendr.

        Returns:
            ee.ImageCollection: A collection where each image represents the medoid of observations per TM-equivalent surface reflectance bands 1-5 and 7, for a given year. There will be as many images as there are years in the range inclusive of the start year and end year.
        """
        dummy_collection = ee.ImageCollection(
            [ee.Image([0, 0, 0, 0, 0, 0]).mask(ee.Image(0))])
        return ee.ImageCollection(
            [self._build_medoid_mosaic(year, dummy_collection, debug) for year in range(self.start_date.year, self.end_date.year + 1)])

    def _build_medoid_mosaic(self,
                             year: int,
                             dummy_collection: ee.ImageCollection,
                             debug: Optional[bool] = False):
        collection = self._get_combined_sr_collection(year)
        collection = collection.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12'])
        
        image_count = collection.size()
        final_collection = ee.ImageCollection(ee.Algorithms.If(image_count.gt(0), collection, dummy_collection))
        
        median = final_collection.median()
        med_diff_collection = final_collection.map(
            lambda image: calculate_median_diff(image, median))
        
        # TODO: selection based on actual names
        medoid_image = med_diff_collection.reduce(ee.Reducer.min(7))
        medoid_image = medoid_image.select([1, 2, 3, 4, 5, 6],
                                           LandsatComposite._band_names)
        
        return medoid_image\
            .set('system:time_start', ee.Date.fromYMD(year, self.start_date.month, self.start_date.day).millis())\
            .toUint16()

    def _get_combined_sr_collection(self, year: int):
        return self._get_sr_collection(year)

    def _get_sr_collection(self, year: int):
        if self.start_date.month > self.end_date.month:
            start_date = ee.Date.fromYMD(
                year - 1, self.start_date.month, self.start_date.day)
            end_date = ee.Date.fromYMD(
                year, self.end_date.month, self.end_date.day)
        else:
            start_date = ee.Date.fromYMD(
                year, self.start_date.month, self.start_date.day)
            end_date = ee.Date.fromYMD(
                year, self.end_date.month, self.end_date.day)
            
        s2_sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')\
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', self.cloud_filter))\
                    .filterBounds(self.area_of_interest)\
                    .filterDate(start_date, end_date)

        s2_cloud_prob = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')\
                        .filterBounds(self.area_of_interest)\
                        .filterDate(start_date, end_date)
                         

        # Join collections to attach cloud probability to each SR image
        collection = self._index_join(s2_sr, s2_cloud_prob, 's2cloudless')\
                    .map(self._preprocess_image)\
                    .set("system:time_start", start_date.millis())
        
        return self._remove_images(collection)

    def _preprocess_image(self, image: int):
        image = self._add_cloud_shdw_mask(image)
        image = self._apply_cloud_shdw_mask(image)
        return image

    def _remove_images(self, collection: ee.ImageCollection):
        """
        Removes images from a collection based on the given exclude criteria.

        Args:
            collection (ee.ImageCollection): The image collection to remove images from.

        Returns:
            ee.ImageCollection: The updated image collection with images removed.
        """
        if 'imageIds' in self.exclude:
            exclude_list = self.exclude['imageIds']
            for image_id in exclude_list:
                collection = collection.filter(ee.Filter.neq(
                    'system:index', image_id.split('/')[-1]))
        return collection

    def _index_join(self, collectionA, collectionB, property_name):
        """
        Joins two collections on the system:index property and attaches the second collection to the first as a band.
        """
        joined = ee.ImageCollection(ee.Join.saveFirst(property_name).apply(
            primary=collectionA,
            secondary=collectionB,
            condition=ee.Filter.equals(
                leftField='system:index', rightField='system:index'
            )
        ))
        return joined.map(lambda image: image.addBands(ee.Image(image.get(property_name))))

    def _add_cloud_bands(self, image):
        """
        Adds cloud bands to the image based on cloud probability and cloud mask.
        """
        cloud_prb = ee.Image(image.get('s2cloudless')).select('probability')
        is_cloud = cloud_prb.gt(self.cloud_prb_thresh).rename('clouds')
        return image.addBands(ee.Image([cloud_prb, is_cloud]))

    def _add_shadow_bands(self, image):
        """
        Adds shadow bands to the image to mask cloud shadows.
        """
        not_water = image.select('SCL').neq(6)
        dark_pixels = image.select('B8').lt(self.nir_drk_thresh * self._sr_band_scale).multiply(not_water).rename('dark_pixels')

        shadow_azimuth = ee.Number(90).subtract(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')))
        cloud_proj = (image.select('clouds').directionalDistanceTransform(shadow_azimuth, self.cloud_prj_dist * 10)
                    .reproject(crs=image.select(0).projection(), scale=100)
                    .select('distance')
                    .mask()
                    .rename('cloud_transform'))

        shadows = cloud_proj.multiply(dark_pixels).rename('shadows')
        return image.addBands(ee.Image([dark_pixels, cloud_proj, shadows]))

    def _add_cloud_shdw_mask(self, image):
        """
        Adds a combined cloud and shadow mask to the image.
        """
        image = self._add_cloud_bands(image)
        image = self._add_shadow_bands(image)

        is_cloud_shdw = image.select('clouds').add(image.select('shadows')).gt(0)
        is_cloud_shdw = (is_cloud_shdw.focalMin(2).focalMax(self.buffer * 2 / 20)
                       .reproject(crs=image.select([0]).projection(), scale=20)
                       .rename('cloudmask'))

        return image.addBands(is_cloud_shdw)

    def _apply_cloud_shdw_mask(self, image):
        """
        Applies the cloud and shadow mask to the image, leaving only clear areas.
        """
        not_cloud_shdw = image.select('cloudmask').Not()
        return image.select('B.*').updateMask(not_cloud_shdw)


class LandsatComposite(ee.ImageCollection):
    """
    ### LandsatComposite class for building a Landsat surface reflectance composite for use in LandTrendr.

    #### Args:
        collection (ee.ImageCollection, optional):  An ee.ImageCollection of Landsat surface reflectance images. Defaults to None.
        start_date (datetime, optional): The start date of the annual composites. Defaults to None.
        end_date (datetime, optional): The end date of the annual composites. Defaults to None.
        area_of_interest (ee.Geometry, optional): The area of interest for the composite. Defaults to None.
        mask_labels (list, optional): A list of mask labels to apply to the composite. Defaults to ['cloud', 'shadow', 'snow'].
        exclude (dict, optional): A dictionary of exclude criteria. Defaults to {}.
        debug (bool, optional): Whether to calculate the number of clear view pixels. Defaults to False.
    """
    _band_names = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']
    _mask_options = ['cloud', 'shadow', 'snow',
                     'water', 'waterplus', 'nonforest']

    def __init__(self,
                 collection: Optional[ee.ImageCollection] = None,
                 start_date: Optional[datetime] = None,
                 end_date: Optional[datetime] = None,
                 area_of_interest: Optional[ee.Geometry] = None,
                 mask_labels: list[str] = ['cloud', 'shadow', 'snow'],
                 exclude: dict = {},
                 debug: bool = False):
        if collection:
            super().__init__(collection)
        else:
            assert start_date, "start_date is required"
            assert end_date, "end_date is required"
            assert area_of_interest, "area_of_interest is required"
            self.start_date = start_date
            self.end_date = end_date
            self.area_of_interest = area_of_interest
            self.mask_labels = mask_labels
            self.exclude = exclude
            super().__init__(self._build_sr_collection(debug))

    @property
    def mask_labels(self):
        return self._mask_labels

    @mask_labels.setter
    def mask_labels(self, mask_labels: list):
        assert all([_ in self._mask_options for _ in mask_labels]
                   ), f"mask_labels must be a subset of {self._mask_options}"
        self._mask_labels = mask_labels

    def _build_sr_collection(self, debug: Optional[bool] = False):
        """
        Builds a medoid composite of Landsat surface reflectance TM-equivalent bands 1,2,3,4,5,7. 
        This collection can be useful outside of use by LandTrendr, but is also the base for creating the input collection for LandTrendr.

        Returns:
            ee.ImageCollection: A collection where each image represents the medoid of observations per TM-equivalent surface reflectance bands 1-5 and 7, for a given year. There will be as many images as there are years in the range inclusive of the start year and end year.
        """
        dummy_collection = ee.ImageCollection(
            [ee.Image([0, 0, 0, 0, 0, 0]).mask(ee.Image(0))])
        return ee.ImageCollection(
            [self._build_medoid_mosaic(year, dummy_collection, debug) for year in range(self.start_date.year, self.end_date.year + 1)])

    def _build_medoid_mosaic(self,
                             year: int,
                             dummy_collection: ee.ImageCollection,
                             debug: Optional[bool] = False):
        collection = self._get_combined_sr_collection(year)
        image_count = collection.size()
        final_collection = ee.ImageCollection(ee.Algorithms.If(
            image_count.gt(0), collection, dummy_collection))
        if debug:
            not_mask_count = ee.ImageCollection(
                [count_clear_view_pixels(final_collection)])
            if getattr(self, 'clear_pixel_count_collection', None) is not None:
                self.clear_pixel_count_collection = self.clear_pixel_count_collection.merge(
                    not_mask_count)
            else:
                self.clear_pixel_count_collection = not_mask_count
        median = final_collection.median()
        med_diff_collection = final_collection.map(
            lambda image: calculate_median_diff(image, median))
        return med_diff_collection\
            .reduce(ee.Reducer.min(7))\
            .select([1, 2, 3, 4, 5, 6], self._band_names)\
            .set('system:time_start', ee.Date.fromYMD(year, self.start_date.month, self.start_date.day).millis())\
            .toUint16()

    def _get_combined_sr_collection(self, year: int):
        lt5 = self._get_sr_collection(year, 'LT05')
        le7 = self._get_sr_collection(year, 'LE07')
        lc8 = self._get_sr_collection(year, 'LC08')
        lc9 = self._get_sr_collection(year, 'LC09')
        return lt5.merge(le7).merge(lc8).merge(lc9)

    def _get_sr_collection(self, year: int, sensor: str):
        if self.start_date.month > self.end_date.month:
            start_date = ee.Date.fromYMD(
                year - 1, self.start_date.month, self.start_date.day)
            end_date = ee.Date.fromYMD(
                year, self.end_date.month, self.end_date.day)
        else:
            start_date = ee.Date.fromYMD(
                year, self.start_date.month, self.start_date.day)
            end_date = ee.Date.fromYMD(
                year, self.end_date.month, self.end_date.day)
        sr_collection = ee.ImageCollection('LANDSAT/' + sensor + '/C02/T1_L2')\
            .filterBounds(self.area_of_interest)\
            .filterDate(start_date, end_date)\
            .map(lambda image: self._preprocess_image(image, sensor))\
            .set("system:time_start", start_date.millis())
        return self._remove_images(sr_collection)

    def _preprocess_image(self, image: int, sensor: str):
        # Accounting for band shift between landsat difference landsat images
        if sensor == 'LC08' or sensor == 'LC09':
            dat = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
                               self._band_names)
        else:
            dat = image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
                               self._band_names)
        dat = self._scale_unmask_image(dat)
        if len(self.mask_labels) > 0:
            dat = self._apply_masks(image.select('QA_PIXEL'), dat)
        return dat

    def _scale_unmask_image(self, image: ee.Image):
        return image.multiply(0.0000275).add(-0.2).multiply(10000).toUint16().unmask()

    def _apply_masks(self, qa: ee.Image, dat: ee.Image):
        mask = ee.Image(1)
        # TODO: Refactor to allow dynamically allow new masks
        for mask_label in self.mask_labels:
            match mask_label:
                case 'water':
                    mask = qa.bitwiseAnd(1 << 7).eq(0).multiply(mask)
                case 'shadow':
                    mask = qa.bitwiseAnd(1 << 4).eq(0).multiply(mask)
                case 'snow':
                    mask = qa.bitwiseAnd(1 << 5).eq(0).multiply(mask)
                case 'cloud':
                    mask = qa.bitwiseAnd(1 << 3).eq(0).multiply(mask)
                case 'waterplus':
                    mask = mask.mask(water_mask(self.area_of_interest))
                case 'nonforest':
                    mask = mask.mask(forest_mask(self.area_of_interest))
        return dat.mask(mask)

    def _remove_images(self, collection: ee.ImageCollection):
        """
        Removes images from a collection based on the given exclude criteria.

        Args:
            collection (ee.ImageCollection): The image collection to remove images from.

        Returns:
            ee.ImageCollection: The updated image collection with images removed.
        """
        if 'imageIds' in self.exclude:
            exclude_list = self.exclude['imageIds']
            for image_id in exclude_list:
                collection = collection.filter(ee.Filter.neq(
                    'system:index', image_id.split('/')[-1]))

        if 'slcOff' in self.exclude:
            if self.exclude['slcOff'] is True:
                collection = collection.filter(ee.Filter.And(
                    ee.Filter.eq('SPACECRAFT_ID', 'LANDSAT_7'),
                    ee.Filter.gt('system:time_start', 1054425600000)
                ).Not())
        return collection


class LtCollection(ee.ImageCollection):
    """

    ### LtCollection class for building a LandTrendr collection for use in LandTrendr.

    #### Args:
    sr_collection (ee.ImageCollection, optional):  An ee.ImageCollection of Landsat surface reflectance images. Defaults to None.
    composite_params (dict, optional):  A dictionary of parameters to build a new LandsatComposite. Defaults to None.
    index (str, optional): The spectral index to base temporal segmentation on. Defaults to 'NBR'.
    ftv_list (list, optional): A list of spectral indices to fit to the segmentation structure of the segmentation index. Defaults to [].

    """
    _index_flip = ['NBR', 'NDVI', 'NDSI', 'NDMI',
                   'EVI', 'TCG', 'TCW', 'TCA', 'B4', 'NDFI']
    _valid_indices = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'NBR',
                      'NDMI', 'NDVI', 'NDSI', 'EVI', 'NDFI', 'TCB', 'TCG', 'TCW', 'TCA']
    _valid_indices_alt = ['TCC', 'TCM', 'TCS', 'ENC', 'ENM', 'ENS', 'ENC1',
                          'ENM1', 'ENS1', 'B5z', 'B7z', 'TCWz', 'TCAz', 'NDMIz', 'NBRz']

    def __init__(self,
                 sr_collection: LandsatComposite | ee.ImageCollection | dict,
                 index: str = 'NBR',
                 ftv_list: Optional[list[str]] = []):
        self.index = index
        self.ftv_list = ftv_list
        if isinstance(sr_collection, dict):
            self.sr_collection = LandsatComposite(**sr_collection)
        else:
            self.sr_collection = sr_collection
        assert self.sr_collection, "sr_collection must be provided or composite_params must be provided to build a new collection."
        self.needs_index_flip = self.index in self._index_flip
        super().__init__(self._build_lt_collection(self.sr_collection))

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, index: str):
        assert index in self._valid_indices or index in self._valid_indices_alt, \
            f"Index must be one of {self._valid_indices} or {self._valid_indices_alt}"
        self._index = index

    @property
    def ftv_list(self):
        return self._ftv_list

    @ftv_list.setter
    def ftv_list(self, ftv_list: list):
        assert all([_ in self._valid_indices for _ in ftv_list]
                   ), f"ftv_list must be a subset of {self._valid_indices}"
        self._ftv_list = ftv_list

    def _build_lt_collection(self,
                             collection: ee.ImageCollection) -> ee.ImageCollection:
        """
        Builds a collection as input to LandTrendr. It will prepare a collection where the first band is the spectral index to base temporal segmentation on, and the subsequent bands will be fitted to segmentation structure of the segmentation index.

        Args:
            collection (ee.ImageCollection): The input annual image collection.

        Returns:
            ee.ImageCollection
        """

        match self.index:
            case 'TCC':
                return self._make_tc_composite(collection, 'mean')
            case 'TCM':
                return self._make_tc_composite(collection, 'max')
            case 'TCS':
                return self._make_tc_composite(collection, 'sum')
            case 'ENC':
                return self._make_ensemble_composite(collection, 'mean')
            case 'ENM':
                return self._make_ensemble_composite(collection, 'max')
            case 'ENS':
                return self._make_ensemble_composite(collection, 'sum')
            case 'ENC1':
                return self._make_ensemble_composite_alt(collection, 'mean')
            case 'ENM1':
                return self._make_ensemble_composite_alt(collection, 'max')
            case 'ENS1':
                return self._make_ensemble_composite_alt(collection, 'sum')
            case 'B5z':
                return self._standardize_index(collection, 'B5')
            case 'B7z':
                return self._standardize_index(collection, 'B7')
            case 'TCWz':
                return self._standardize_index(collection, 'TCW')
            case 'TCAz':
                return self._standardize_index(collection, 'TCA')
            case 'NDMIz':
                return self._standardize_index(collection, 'NDMI')
            case 'NBRz':
                return self._standardize_index(collection, 'NBR')
            case _:
                return self._transform_sr_collection(collection, self.ftv_list, "ftv")

    def _transform_sr_collection(self,
                                 collection: ee.ImageCollection,
                                 band_list: list[str],
                                 prefix: str = None) -> ee.ImageCollection:
        """
        Transforms the images within an annual surface reflectance collection built by buildSRcollection to a list of provided indices or bands.

        Args:
            collection (ee.ImageCollection): The surface reflectance collection to transform.
            band_list (list): The list of bands to transform the collection to.

        Returns:
            ee.ImageCollection: The transformed surface reflectance collection that includes one image per year based on an image collection built by buildSRcollection function transformed to the indices provided in the bandList parameter..
        """
        return collection.map(lambda image: self._make_default_composite(image, band_list, prefix))

    def _make_default_composite(self,
                                image: ee.Image,
                                band_list: list[str],
                                prefix: str = None) -> ee.Image:
        """
        Generates a default feature composite for the given image.

        Args:
            image (ee.Image): The input image.
            band_list (list[str]): The list of bands to include in the composite.
            prefix (str): The prefix to add to the band names.

        Returns:
            ee.Image: The feature composite image.
        """
        all_stack = self._calculate_index(image, self.index)
        for band in band_list:
            band_image = self._calculate_index(image, band, False)
            if prefix:
                band_image = band_image.select(
                    [band], [prefix + '_' + band.lower()])
            all_stack = all_stack.addBands(band_image).set(
                'system:time_start', image.get('system:time_start'))
        return all_stack

    def _standardize_index(self,
                           collection: ee.ImageCollection,
                           index: str) -> ee.ImageCollection:
        """
        Standardizes the images in a collection based on a given index.

        Args:
            collection (ee.ImageCollection): The image collection.
            index (str): The index to standardize.

        Returns:
            ee.ImageCollection: The standardized image collection.
        """
        z_collection = collection.map(
            lambda image: self._calculate_index(image, index))
        return standardize_collection(z_collection).map(lambda image: image.multiply(1000).set('system:time_start', image.get('system:time_start')))

    def _calculate_index(self,
                         image: ee.Image,
                         index: str,
                         flip: bool = True) -> ee.Image:
        index = index.upper()
        match index:
            case 'B1' | 'B2' | 'B3' | 'B5' | 'B4' | 'B7':
                index_image = image.select(index).toFloat()
            case 'NBR':
                index_image = nbr_transform(image)
            case 'NDMI':
                index_image = ndmi_transform(image)
            case 'NDVI':
                index_image = ndvi_transform(image)
            case 'NDSI':
                index_image = ndsi_transform(image)
            case 'EVI':
                index_image = evi_transform(image)
            case 'TCB' | 'TCG' | 'TCW' | 'TCA':
                index_image = tc_transform(image).select([index])
            case 'NDFI':
                index_image = ndfi_transform(image)
            case _:
                raise ValueError('The index you provided is not supported')
        if flip and self.needs_index_flip:
            index_image = index_image.multiply(-1)
        return index_image.set('system:time_start', image.get('system:time_start'))

    def _reducer(image_collection: ee.ImageCollection,
                 reducer: str) -> ee.Image:
        match reducer:
            case 'mean':
                return image_collection.mean()
            case 'max':
                return image_collection.max()
            case 'sum':
                return image_collection.sum()
            case _:
                raise ValueError('The reducer you provided is not supported')

    def _make_tc_composite(self,
                           collection: ee.ImageCollection,
                           reducer: str) -> ee.ImageCollection:
        """
        Creates a tasseled-cap (brightness, greenness, and wetness) composite for the given collection using the specified reducer.

        Args:
            collection (ee.ImageCollection): The image collection to create the composite from.
            reducer (str): The reducer to use for compositing.

        Returns:
            ee.ImageCollection: The tasseled-cap composite.
        """
        tc_composite = collection.map(self._tc_composite)
        tcb = tc_composite.select('TCB')
        tcg = tc_composite.select('TCG')
        tcw = tc_composite.select('TCW')

        tcb_standard = standardize_collection(tcb)
        tcg_standard = standardize_collection(tcg)
        tcw_standard = standardize_collection(tcw)
        tc_standard = tcb_standard.combine(tcg_standard).combine(tcw_standard)

        return tc_standard.map(lambda image: self._tc_reducer(image, reducer))

    def _tc_composite(self,
                      image: ee.Image) -> ee.Image:
        tcb = self._calculate_index(image, 'TCB')
        tcg = self._calculate_index(image, 'TCG')
        tcw = self._calculate_index(image, 'TCW')
        return tcb.addBands(tcg).addBands(tcw).set('system:time_start', image.get('system:time_start'))

    def _tc_reducer(self,
                    image: ee.Image,
                    reducer: str) -> ee.Image:
        image_collection = ee.ImageCollection.fromImages(
            [image.select(['TCB'], ['Z']),
             image.select(['TCG'], ['Z']),
             image.select(['TCW'], ['Z'])]
        )
        reduced_image = self._reducer(image_collection, reducer)
        return reduced_image.multiply(1000).set('system:time_start', image.get('system:time_start'))

    def _make_ensemble_composite(self,
                                 collection: ee.ImageCollection,
                                 reducer: str) -> ee.Image:
        """
        Creates an ensemble (see bands below) composite for the given collection using the specified reducer.

        Args:
            collection (ee.ImageCollection): The input image collection.
            reducer (str): The reducer to apply when combining the bands.

        Returns:
            ee.Image: The ensemble composite image.
        """

        stack = collection.map(self._ensemble_composite)

        b5 = stack.select('B5')
        b7 = stack.select('B7')
        tcw = stack.select('TCW')
        tca = stack.select('TCA')
        ndmi = stack.select('NDMI')
        nbr = stack.select('NBR')

        b5_standard = standardize_collection(b5)
        b7_standard = standardize_collection(b7)
        tcw_standard = standardize_collection(tcw)
        tca_standard = standardize_collection(tca)
        ndmi_standard = standardize_collection(ndmi)
        nbr_standard = standardize_collection(nbr)
        ensemble = b5_standard.combine(b7_standard)\
            .combine(tcw_standard)\
            .combine(tca_standard)\
            .combine(ndmi_standard)\
            .combine(nbr_standard)

        return ensemble.map(lambda image: self._ensemble_reducer(image, reducer))

    def _ensemble_composite(self,
                            image: ee.Image) -> ee.Image:
        b5 = self._calculate_index(image, 'B5')
        b7 = self._calculate_index(image, 'B7')
        tcw = self._calculate_index(image, 'TCW')
        tca = self._calculate_index(image, 'TCA')
        ndmi = self._calculate_index(image, 'NDMI')
        nbr = self._calculate_index(image, 'NBR')

        return b5.addBands(b7)\
            .addBands(tcw)\
            .addBands(tca)\
            .addBands(ndmi)\
            .addBands(nbr)\
            .set('system:time_start', image.get('system:time_start'))

    def _ensemble_reducer(self,
                          image: ee.Image,
                          reducer: str) -> ee.Image:
        image_collection = ee.ImageCollection.fromImages(
            [image.select(['B5'], ['Z']),
             image.select(['B7'], ['Z']),
             image.select(['TCW'], ['Z']),
             image.select(['TCA'], ['Z']),
             image.select(['NDMI'], ['Z']),
             image.select(['NBR'], ['Z'])]
        )
        reduced_image = self._reducer(image_collection, reducer)
        return reduced_image.multiply(1000).set('system:time_start', image.get('system:time_start'))

    def _make_ensemble_composite_alt(self,
                                     collection: ee.ImageCollection,
                                     reducer: str) -> ee.ImageCollection:
        """
        Creates an ensemble (see bands below) composite for the given collection using the specified reducer.

        Args:
            collection (ee.ImageCollection): The input image collection.
            reducer (str): The reducer to apply on the ensemble composite.

        Returns:
            ee.Image: The ensemble composite image.
        """
        e_composite = collection.map(self._ensemble_composite_alt)
        b5 = e_composite.select('B5')
        tcb = e_composite.select('TCB')
        tcg = e_composite.select('TCG')
        nbr = e_composite.select('NBR')

        b5_standard = standardize_collection(b5)
        tcb_standard = standardize_collection(tcb)
        tcg_standard = standardize_collection(tcg)
        nbr_standard = standardize_collection(nbr)
        ensemble_alt = b5_standard.combine(tcb_standard)\
            .combine(tcg_standard)\
            .combine(nbr_standard)

        return ensemble_alt.map(lambda image: self._ensemble_reducer_alt(image, reducer))

    def _ensemble_composite_alt(self,
                                image: ee.Image) -> ee.Image:
        b5 = self._calculate_index(image, 'B5')
        tcb = self._calculate_index(image, 'TCB')
        tcg = self._calculate_index(image, 'TCG')
        nbr = self._calculate_index(image, 'NBR')

        return b5.addBands(tcb)\
            .addBands(tcg)\
            .addBands(nbr)\
            .set('system:time_start', image.get('system:time_start'))

    def _ensemble_reducer_alt(self,
                              image: ee.Image,
                              reducer: str) -> ee.Image:
        image_collection = ee.ImageCollection.fromImages(
            [image.select(['B5'], ['Z']),
             image.select(['TCB'], ['Z']),
             image.select(['TCG'], ['Z']),
             image.select(['NBR'], ['Z'])]
        )
        reduced_image = self._reducer(image_collection, reducer)
        return reduced_image.multiply(1000).set('system:time_start', image.get('system:time_start'))


class LandTrendr(ee.Image):
    """
    ### Landtrendr class for performing landtrendr analysis using Google Earth Engine.

    #### Args:
    lt_collection (:obj:, optional):         LandTrendr collection. Defaults to None.
    run_params (:obj:, optional):            Dictionary of run parameters (See below).
        {
            maxSegments (int, optional):                Maximum number of segments. Defaults to 6.
            spikeThreshold (float, optional):           Spike threshold. Defaults to 0.9.
            vertexCountOvershoot (int, optional):       Vertex count overshoot. Defaults to 3.
            preventOneYearRecovery (bool, optional):    Prevent one year recovery. Defaults to True.
            recoveryThreshold (float, optional):        Recovery threshold. Defaults to 0.25.
            pvalThreshold (float, optional):            P-value threshold. Defaults to 0.25.
            bestModelProportion (float, optional):      Best model proportion. Defaults to 0.75.
            minObservationsNeeded (int, optional):      Minimum number of observations needed. Defaults to 6.
        }
    run (bool, optional):                               Whether to run the run method immediately. Defaults to True.
    """

    _default_run_params = {
        "maxSegments": 6,
        "spikeThreshold": 0.9,
        "vertexCountOvershoot":  3,
        "preventOneYearRecovery":  False,
        "recoveryThreshold":  0.25,
        "pvalThreshold":  0.1,
        "bestModelProportion":  1.25,
        "minObservationsNeeded": 6,
    }

    def __init__(self,
                 lt_collection: LtCollection | ee.ImageCollection | dict,
                 run_params: dict = _default_run_params,
                 run: bool = True):
        self.run_params = run_params
        self.initialized = False
        if isinstance(lt_collection, dict):
            self.lt_collection = LtCollection(**lt_collection)
        else:
            self.lt_collection = lt_collection

        if run:
            self.run()

    def __getattribute__(self, name):
        attr = super().__getattribute__(name)

        # 'initilize' is the run method of ee.Image
        if callable(attr) and name not in ["run", "initialize"] and not self.initialized:
            raise AttributeError(
                "LandTrendr object has not been initialized. Please run the 'run' method.")
        else:
            return attr

    @property
    def run_params(self):
        return self._run_params

    @run_params.setter
    def run_params(self, run_params: dict):
        assert all([_ in self._default_run_params.keys()
                    for _ in run_params.keys()]), f"run_params must be a subset of {self._default_run_params.keys()}"
        if hasattr(self, '_data'):
            self._run_params |= run_params
        else:
            self._run_params = self._default_run_params | run_params
            self._data = None

    def run(self):
        """
        Initiates the LandTrendr algorithm on Google's servers using the specified run_params and generates an image.
        """
        super().__init__(ee.Algorithms.TemporalSegmentation.LandTrendr(
            timeSeries=self.lt_collection, **self.run_params))
        self.initialized = True

    def get_change_map(self,
                       change_params: dict) -> ee.Image:
        """
        Generates a set of map layers describing either vegetation loss or gain events with attributes including: year of change detection, spectral delta, duration of change event, pre-change event spectral value, and the rate of spectral change. Each attribute is a band of an ee.Image.

        Args:
            change_params (dict): A dictionary containing the parameters for change detection.
                {
                    'delta':                'loss' | 'gain' | 'all',
                    'sort' (optional):      'greatest' | 'least' | 'newest' | 'oldest' | 'fastest' | 'slowest',
                    'years' (optional):     {'start': int, 'end': int},
                    'mag' (optional):       {'value': int, 'operator': '>' | '<' , 'dsnr': bool},
                    'dur' (optional):       {'value': int, 'operator': '>' | '<'},
                    'preval' (optional):    {'value': int, 'operator': '>' | '<'},
                    'mmu' (optional):       {'value': int > 1}
                }

        Returns:
            ee.Image: An image with bands for attributes of change events meeting filtering criteria including:
                Year of change event detection: 'yod' (year)
                Magnitude of change event: 'mag' (absolute value of change event spectral delta)
                Duration of change event: 'dur' (years)
                Pre-change event spectral value: 'preval' (spectral value)
                Rate of spectral change for event 'rate' (mag/dur)
                DSNR 'dsnr' (mag/fit rmse) multiplied by 100 to retain two decimal precision with Int16 data.
        """
        # Backward compatibility for dsnr
        if 'dsnr' not in change_params['mag']:
            change_params['mag']['dsnr'] = False

        # Get the segment info
        assert change_params['delta'] in \
            ['loss', 'gain', 'all'], \
            "delta must be one of 'loss', 'gain', or 'all'"
        seg_info = self.get_segment_data(
            change_params['delta'], self.lt_collection.needs_index_flip)
        change_mask = seg_info.arraySlice(0, 4, 5).gt(0)
        seg_info = seg_info.arrayMask(change_mask)

        # Filter by year
        if 'years' in change_params:
            yod_arr = seg_info.arraySlice(0, 0, 1).add(1)
            year_mask = yod_arr.gte(ee.Number(change_params['years']['start'])).And(
                yod_arr.lte(change_params['years']['end']))
            seg_info = seg_info.arrayMask(year_mask)

        # Filter by mag
        mag_band = {'axis': 0, 'start': 4, 'end': 5}
        if 'mag' in change_params:
            assert change_params['mag']['operator'] in \
                ['>', '<'], "mag operator must be one of '>' or '<'"
            if change_params['mag']['dsnr']:
                mag_band = {'axis': 0, 'start': 7, 'end': None}
            match change_params['mag']['operator']:
                case '<':
                    mag_mask = seg_info.arraySlice(**mag_band).lt(
                        change_params['mag']['value'])
                case '>':
                    mag_mask = seg_info.arraySlice(**mag_band).gt(
                        change_params['mag']['value'])
                case _:
                    raise ValueError(
                        "provided mag operator does not match either '>' or '<'")
            seg_info = seg_info.arrayMask(mag_mask)

        # Filter by dur
        if 'dur' in change_params:
            assert change_params['dur']['operator'] in \
                ['>', '<'], "dur operator must be one of '>' or '<'"
            dur_band = {'axis': 0, 'start': 5, 'end': 6}
            match change_params['dur']['operator']:
                case '<':
                    dur_mask = seg_info.arraySlice(**dur_band).lt(
                        change_params['dur']['value'])
                case '>':
                    dur_mask = seg_info.arraySlice(**dur_band).gt(
                        change_params['dur']['value'])
                case _:
                    raise ValueError(
                        "The dur operator does not match either '>' or '<'")
            seg_info = seg_info.arrayMask(dur_mask)

        # Filter by preval
        if 'preval' in change_params:
            assert change_params['preval']['operator'] in \
                ['>', '<'], "preval operator must be one of '>' or '<'"
            preval_band = {'axis': 0, 'start': 2, 'end': 3}
            match change_params['preval']['operator']:
                case '<':
                    preval_mask = seg_info.arraySlice(**preval_band).lt(
                        change_params['preval']['value'])
                case '>':
                    preval_mask = seg_info.arraySlice(**preval_band).gt(
                        change_params['preval']['value'])
                case _:
                    raise ValueError(
                        "The preval operator does not match either '>' or '<'")
            seg_info = seg_info.arrayMask(preval_mask)

        # Sort by dist type
        if 'sort' in change_params:
            assert change_params['sort'] in \
                ['greatest', 'least', 'newest', 'oldest', 'fastest', 'slowest'], \
                "sort must be one of 'greatest', 'least', 'newest', 'oldest', 'fastest', or 'slowest'"
            match change_params['sort']:
                case 'greatest':
                    sort_by_this = seg_info.arraySlice(0, 4, 5).multiply(-1)
                case 'least':
                    sort_by_this = seg_info.arraySlice(0, 4, 5)
                case 'newest':
                    sort_by_this = seg_info.arraySlice(0, 0, 1).multiply(-1)
                case 'oldest':
                    sort_by_this = seg_info.arraySlice(0, 0, 1)
                case 'fastest':
                    sort_by_this = seg_info.arraySlice(0, 5, 6)
                case 'slowest':
                    sort_by_this = seg_info.arraySlice(0, 5, 6).multiply(-1)
            seg_info = seg_info.arraySort(sort_by_this)

        change_array = seg_info.arraySlice(1, 0, 1)

        # Make an image from the array of attributes for the change of interest
        arr_row_names = [['startYear', 'end year', 'preval',
                          'postval', 'mag', 'dur', 'rate', 'dsnr']]
        change_image = change_array\
            .arrayProject([0]).arrayFlatten(arr_row_names)
        yod = change_image.select('startYear').add(1).toInt16().rename('yod')
        change_image = change_image.addBands(yod).select(
            ['yod', 'mag', 'dur', 'preval', 'rate', 'dsnr'])

        # Mask for change/no change
        change_image = change_image\
            .updateMask(change_image.select('mag').gt(0))

        # Filter by MMU on year of change detection
        if 'mmu' in change_params:
            assert change_params['mmu']['value'] >= 1, "mmu value must be greater than 1"
            mmu_lyr = change_image.select('yod')
            mmu_mask = self._apply_mmu(mmu_lyr, change_params['mmu']['value'])
            change_image = change_image.updateMask(mmu_mask)

        return change_image

    def get_segment_data(self,
                         delta: int,
                         index_flip: bool,
                         options: Optional[dict] = None) -> ee.Image:
        """
        Generates an array of information about spectral-temporal segments from the breakpoint vertices identified by LandTrendr. Returns either all spectral-temporal segments, or just vegetation loss segments, or just vegetation growth segments.

        Args:
            index (int): The index of the segment.
            delta (str): The type of delta to calculate. Can be 'all', 'gain', or 'loss'.
            options (bool, optional): Additional options for the segment data retrieval. Defaults to None.

        Returns:
            ee.Image: An image array with dimensions: 8 (rows) x nSegments (cols). Each row describes an attribute of the segments identified by LandTrendr per pixel time series. Each column represents a segment in the time series per pixel ordered from earliest to latest in the series.
                Row 1: segment start year
                Row 2: segment end year
                Row 3: segment start value
                Row 4: segment end value
                Row 5: segment spectral delta
                Row 6: segment duration
                Row 7: segment rate of spectral change
                Row 8: segment DSNR*
        """

        # Deal with options
        _options = {'right': False}  # Defaults
        if options is not None and isinstance(options, bool):
            _options['right'] = options

        # TODO: Refactor to avoid excessive array slicing
        lt_band = self.select('LandTrendr')  # select the LandTrendr band
        rmse = self.select('rmse')  # select the rmse band
        # slice out the 'Is Vertex' row - yes(1)/no(0)
        vertex_mask = lt_band.arraySlice(0, 3, 4)
        # use the 'Is Vertex' row as a mask for all rows
        vertices = lt_band.arrayMask(vertex_mask)
        # slice out the vertices as the start of segments
        left_list = vertices.arraySlice(1, 0, -1)
        # slice out the vertices as the end of segments
        right_list = vertices.arraySlice(1, 1, None)
        # get year dimension of LT data from the segment start vertices
        start_year = left_list.arraySlice(0, 0, 1)
        # get spectral index dimension of LT data from the segment start vertices
        start_val = left_list.arraySlice(0, 2, 3)
        # get year dimension of LT data from the segment end vertices
        end_year = right_list.arraySlice(0, 0, 1)
        # get spectral index dimension of LT data from the segment end vertices
        end_val = right_list.arraySlice(0, 2, 3)
        # subtract the segment start year from the segment end year to calculate the duration of segments
        dur = end_year.subtract(start_year)
        # subtract the segment start index value from the segment end index value to calculate the delta of segments
        mag = end_val.subtract(start_val)
        rate = mag.divide(dur)  # calculate the rate of spectral change
        dsnr = mag.divide(rmse)  # make mag relative to fit rmse

        # Whether to return all segments or either dist or grow
        match delta:
            case 'all':
                # If the data should be set to the correct orientation, adjust it
                if _options['right']:
                    if index_flip:
                        start_val = start_val.multiply(-1)
                        end_val = end_val.multiply(-1)
                        mag = mag.multiply(-1)
                        rate = rate.multiply(-1)
                        dsnr = dsnr.multiply(-1)

                # Now just get out - return the result
                return ee.Image.cat([start_year, end_year, start_val, end_val, mag, dur, rate, dsnr]) \
                    .unmask(ee.Image(ee.Array([[-9999]]))) \
                    .toArray(0)

            case 'gain' | 'loss':
                match delta:
                    case 'gain':
                        change_type_mask = mag.lt(0)
                    case 'loss':
                        change_type_mask = mag.gt(0)

                flip = -1 if index_flip else 1
                return ee.Image.cat([
                    start_year.arrayMask(change_type_mask),
                    end_year.arrayMask(change_type_mask),
                    start_val.arrayMask(change_type_mask).multiply(flip),
                    end_val.arrayMask(change_type_mask).multiply(flip),
                    mag.arrayMask(change_type_mask).abs(),
                    dur.arrayMask(change_type_mask),
                    rate.arrayMask(change_type_mask).abs(),
                    dsnr.arrayMask(change_type_mask).abs(),
                ]) \
                    .unmask(ee.Image(ee.Array([[-9999]]))) \
                    .toArray(0)

    @staticmethod
    def get_segment_count(segment_data: ee.Image) -> ee.Image:
        """
        Given a segment data array produced by the getSegmentData function, this function returns the number of segments identified by LandTrendr as an ee.Image.

        Args:
            segment_data (ee.Image): an image array returned from the get_segment_data function

        Returns:
            ee.Image: A single-band ee.Image describing the number of segments per pixel time series identified by LandTrendr.
        """
        return segment_data.arrayLength(1).select([0], ['segCount']).toByte()

    def get_fitted_data(self,
                        index: str,
                        start_date: datetime,
                        end_date: datetime) -> ee.Image:
        """
        Generates an annual band stack for a given index provided as ftvList indices to either buildLTcollection or runLT. It flattens the FTV array format to a band per year for a given FTV index.

        Args:
            index (str): The index for which to retrieve the fitted data.
            start_date (datetime): The start date of the time series.
            end_date (datetime): The end date of the time series.

        Returns:
            ee.Image: An image representing fitted-to-vertex annual spectral data for whatever index was provided as the index parameter. There will be as many bands as there are years in the range inclusive of the start year and end year.
        """
        search = 'ftv_' + index.lower() + '_fit'
        return self.select(search).arrayFlatten([[str(_) for _ in range(start_date.year, end_date.year + 1)]])

    def collection_to_band_stack(self,
                                 collection: ee.ImageCollection,
                                 start_date: datetime,
                                 end_date: datetime,
                                 mask_fill: int = 0) -> ee.Image:
        """
        Transforms an image collection into an image stack where each band of each image in the collection is concatenated as a band into a single image. Useful for mapping a function over a collection, like transforming surface reflectance to NDVI, and then transforming the resulting collection into a band sequential time series image stack.

        Args:
            collection (ee.ImageCollection): The Earth Engine image collection to convert.
            start_date (datetime): The start date of the time series.
            end_date (datetime): The end date of the time series.
            mask_fill (int, optional): The value to fill masked pixels with. Default is 0.

        Returns:
            ee.Image: The band stack image representing a band sequential time series of image bands from each image in the given collection between the start year and end year. Note that masked values in the image collection will be filled with 0

        """
        unmasked_collection = collection.map(
            lambda image: image.unmask(mask_fill))
        collection_array = unmasked_collection.toArrayPerBand()
        bands = unmasked_collection.first().bandNames().getInfo()
        all_stack = ee.Image()

        for band in bands:
            band_ts = collection_array.select(band).arrayFlatten(
                [[str(_) for _ in range(start_date.year, end_date.year + 1)]])
            all_stack = ee.Image.cat([all_stack, band_ts])

        return all_stack.slice(1, None).toUint16()

    def get_fitted_rgb_col(self,
                           bands: list[str],
                           vis_params: dict,
                           start_date: datetime,
                           end_date: datetime) -> ee.ImageCollection:
        """
        Creates a collection of RGB visualization images from three FTV bands resulting from a call to LandTrendr segmentation. This is useful for creating thumbnails, filmstrips, and GIFs.

        Args:
            lt (LandTrendr): The LandTrendr object.
            bands (list): A list of band names.
            vis_params (dict): Visualization parameters for the RGB image.

        Returns:
            ee.ImageCollection: An image collection with an RGB image for each year between and including the start year and end year.
        """

        r = self.get_fitted_data(
            bands[0], start_date=start_date, end_date=end_date)
        g = self.get_fitted_data(
            bands[1], start_date=start_date, end_date=end_date)
        b = self.get_fitted_data(
            bands[2], start_date=start_date, end_date=end_date)
        rgb_list = []
        for year in range(start_date.year, end_date.year + 1):
            year_str = str(year)
            rgb_list.append(r.select(year_str)
                            .addBands(g.select(year_str))
                            .addBands(b.select(year_str))
                            .rename(['R', 'G', 'B']))
        rgb_col = ee.ImageCollection(rgb_list)\
            .map(lambda image: image.visualize(**vis_params))\
            .map(lambda image:
                 image.set({'system:time_start': ee.Date.fromYMD(start_date.year, start_date.month, start_date.day).millis(), 'composite_year': start_date.year}))

        return rgb_col

    def _apply_mmu(self, image: ee.Image, mmu_value: int) -> ee.Image:
        mmu_image = image.select([0])\
            .gte(ee.Number(1))\
            .selfMask()\
            .connectedPixelCount()
        min_area = mmu_image.gte(ee.Number(mmu_value)).selfMask()
        return min_area.reproject(image.projection().atScale(30)).unmask()

    # Function to extract vertices from LandTrendr results and stack bands
    def getLTvertStack(lt, runParams):
        lt = lt.select('LandTrendr')
    
        emptyArray = []  # Empty array to hold another array whose length will vary depending on maxSegments parameter
        vertLabels = []  # Empty array to hold band names whose length will vary depending on maxSegments parameter
    
        # Loop through the maximum number of vertices in segmentation and fill empty arrays
        for i in range(1, runParams['maxSegments'] + 2):
            vertLabels.append("vert_" + str(i))  # Make a band name for the given vertex
            emptyArray.append(0)  # Fill in emptyArray
    
        # Make an image to fill holes in result 'LandTrendr' array where vertices found are not equal to maxSegments + 1
        zeros = ee.Image(ee.Array([emptyArray, emptyArray, emptyArray]))
    
        # Labels for 2 dimensions of the array that will be cast to each other in the final step of creating the vertex output
        lbls = [['yrs_', 'src_', 'fit_'], vertLabels]
    
        # Slice out the 4th row, which identifies vertices
        vmask = lt.arraySlice(0, 3, 4)
    
        # Create the vertex stack
        ltVertStack = (lt.arrayMask(vmask)  # Uses the sliced-out isVert row as a mask to include only vertices
                         .arraySlice(0, 0, 3)  # Slices out the vert year row, raw spectral row, and fitted spectral row
                         .addBands(zeros)  # Adds the zeros matrix to ensure 7 vertex slots are represented
                         .toArray(1)  # Concatenates the zeros matrix to the vertex data for 7 vertex slots
                         .arraySlice(1, 0, runParams['maxSegments'] + 1)  # Truncates the columns at the max segments allowed
                         .arrayFlatten(lbls, '')  # Flattens the 2D array into a 1D array with bands for vertices
                       )
    
        return ltVertStack  # Return the stack
