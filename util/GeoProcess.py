import fiona
from shapely.geometry import shape


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
    us_shape_name = "./data/geo/US/cb_2018_us_state_20m.shp"
    us_states = fiona.open(us_shape_name)
    state_geo = None
    for us_state in us_states:
        cur_name = us_state['properties']['NAME']
        if state_name == cur_name:
            state_geo = shape(us_state['geometry'])
            break
    return state_geo
