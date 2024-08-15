import fiona
from shapely.geometry import shape, Point
import numpy as np
from pathlib import Path
from shapely.ops import transform

def StatePolygon(state_name):
    """
    List of States:
        Maryland, Iowa, Delaware, Ohio, Pennsylvania, Nebraska, Washington, Puerto Rico, Alabama,
        Arkansas, New Mexico, Texas, California, Kentucky, Georgia, Wisconsin, Oregon,
        Missouri, Virginia, Tennessee, Louisiana, New York, Michigan, Idaho, Florida,
        Alaska, Illinois, Montana, Minnesota, Indiana, Massachusetts, Kansas, Nevada,
        Vermont, Connecticut, New Jersey, District of Columbia, North Carolina, Utah,
        North Dakota, South Carolina, Mississippi, Colorado, South Dakota, Oklahoma,
        Wyoming, West Virginia, Maine, Hawaii, New Hampshire, Arizona, Rhode Island
    :param state_name: input the state name in the list
    :return: polygon of the state
    """
    us_shape_name = Path(__file__).parent.parent / 'data/geo/US/cb_2018_us_state_20m.shp'
    us_states = fiona.open(us_shape_name)
    state_geo = None
    for us_state in us_states:
        cur_name = us_state['properties']['NAME']
        if state_name == cur_name:
            state_geo = shape(us_state['geometry'])
            break
    return state_geo


def utc_zone(lon, lat):
    """

    :param lon: Longitude of a location
    :param lat: Latitude of a location
    :return: UTC offset (UTC + return_value = LST)
    """
    utc_file = Path(__file__).parent.parent / 'data/geo/timezone/timezone.shp'
    utc_zones = fiona.open(utc_file)
    fire_point = Point(lon, lat)
    zone_value_res = None
    for utc_zone in utc_zones:
        zone_value = utc_zone['properties']['ZONE']
        zone_geo = shape(utc_zone['geometry'])
        if zone_geo.contains(fire_point):
            zone_value_res = zone_value
            break
    return zone_value_res


def UTC_offset(lon, lat):
    """_summary_

    Args:
        lon (_type_): _description_
        lat (_type_): _description_

    Returns:
        _type_: _description_
    """
    m, n = lon.shape
    lon_flatten = lon.flatten()
    lat_flatten = lat.flatten()
    utc_offset_matrix = np.zeros(m * n)
    for i in range(0, len(lon_flatten)):
        current_lon = lon_flatten[i]
        current_lat = lat_flatten[i]
        current_utc_offset = utc_zone(current_lon, current_lat)
        utc_offset_matrix[i] = current_utc_offset
    utc_offset_matrix = np.reshape(utc_offset_matrix, (m, n))
    return utc_offset_matrix


def interest_domain_bound(polygon_list, tol=0.01, crs=None):
    if crs is None:
        x_bound = []
        y_bound = []
        for cur_poly in polygon_list:
            cur_bound = cur_poly.bounds
            x_bound.append(cur_bound[0])
            x_bound.append(cur_bound[2])
            y_bound.append(cur_bound[1])
            y_bound.append(cur_bound[3])
        x_bound = np.array(x_bound)
        y_bound = np.array(y_bound)
        x_min, x_max, y_min, y_max = np.min(x_bound), np.max(x_bound), np.min(y_bound), np.max(y_bound)
        x_min = min([(1 - tol) * x_min, (1 + tol) * x_min])
        x_max = max([(1 - tol) * x_max, (1 + tol) * x_max])
        y_min = min([(1 - tol) * y_min, (1 + tol) * y_min])
        y_max = max([(1 - tol) * y_max, (1 + tol) * y_max])
    else:
        x_bound = []
        y_bound = []
        for cur_poly in polygon_list:
            cur_poly = transform(crs, cur_poly)
            cur_bound = cur_poly.bounds
            x_bound.append(cur_bound[0])
            x_bound.append(cur_bound[2])
            y_bound.append(cur_bound[1])
            y_bound.append(cur_bound[3])
        x_bound = np.array(x_bound)
        y_bound = np.array(y_bound)
        x_min, x_max, y_min, y_max = np.min(x_bound), np.max(x_bound), np.min(y_bound), np.max(y_bound)
        x_min = min([(1 - tol) * x_min, (1 + tol) * x_min])
        x_max = max([(1 - tol) * x_max, (1 + tol) * x_max])
        y_min = min([(1 - tol) * y_min, (1 + tol) * y_min])
        y_max = max([(1 - tol) * y_max, (1 + tol) * y_max])
    return x_min, x_max, y_min, y_max