# This file contains the functions used in the Paleogeography_Analysis notebook within this repository.

import numpy as np
from numpy.core.umath_tests import inner1d
import pandas as pd
import matplotlib.pyplot as plt
import math
import pygplates
import cartopy
import cartopy.crs as ccrs
from shapely.geometry.polygon import Polygon





def lat_from_pole(ref_loc_lon, ref_loc_lat, pole_plon, pole_plat):
    """
    Calculate paleolatitude for a reference location based on a paleomagnetic pole

    Required Parameters
    ----------
    ref_loc_lon: longitude of reference location in degrees
    ref_loc_lat: latitude of reference location
    pole_plon: paleopole longitude in degrees
    pole_plat: paleopole latitude in degrees
    """

    ref_loc = (ref_loc_lon, ref_loc_lat)
    pole = (pole_plon, pole_plat)
    paleo_lat = 90 - angle(pole, ref_loc)
    return float(paleo_lat)





def angle(D1, D2):
    """
    Calculate the angle between two directions.

    Parameters
    ----------
    D1 : Direction 1 as (declination, inclination)
    D2 : Direction 2 as (declination, inclination)

    Returns
    -------
    angle : angle between the directions as a single-element array

    Examples
    --------
    >>> pmag.angle([350.0,10.0],[320.0,20.0])
    array([ 30.59060998])
    """
    D1 = np.array(D1)
    if len(D1.shape) > 1:
        D1 = D1[:, 0:2]  # strip off intensity
    else:
        D1 = D1[:2]
    D2 = np.array(D2)
    if len(D2.shape) > 1:
        D2 = D2[:, 0:2]  # strip off intensity
    else:
        D2 = D2[:2]
    X1 = dir2cart(D1)  # convert to cartesian from polar
    X2 = dir2cart(D2)
    angles = []  # set up a list for angles
    for k in range(X1.shape[0]):  # single vector
        angle = np.arccos(np.dot(X1[k], X2[k])) * \
            180. / np.pi  # take the dot product
        angle = angle % 360.
        angles.append(angle)
    return np.array(angles)





def dir2cart(d):
    """
    Converts a list or array of vector directions in degrees (declination,
    inclination) to an array of the direction in cartesian coordinates (x,y,z)

    Parameters
    ----------
    d : list or array of [dec,inc] or [dec,inc,intensity]

    Returns
    -------
    cart : array of [x,y,z]

    Examples
    --------
    >>> pmag.dir2cart([200,40,1])
    array([-0.71984631, -0.26200263,  0.64278761])
    """
    ints = np.ones(len(d)).transpose(
    )  # get an array of ones to plug into dec,inc pairs
    d = np.array(d)
    rad = np.pi/180.
    if len(d.shape) > 1:  # array of vectors
        decs, incs = d[:, 0] * rad, d[:, 1] * rad
        if d.shape[1] == 3:
            ints = d[:, 2]  # take the given lengths
    else:  # single vector
        decs, incs = np.array(float(d[0])) * rad, np.array(float(d[1])) * rad
        if len(d) == 3:
            ints = np.array(d[2])
        else:
            ints = np.array([1.])
    cart = np.array([ints * np.cos(decs) * np.cos(incs), ints *
                     np.sin(decs) * np.cos(incs), ints * np.sin(incs)]).transpose()
    return cart





def lat_lon_2_cart(lat, lon):
    """
    Convert lat/lon coordinates to cartesian.

    inputs:
    - lat = latitude (-90 to 90)
    - lon = longitude (-180 to 180)

    outputs:
    - cart = (x, y, z) cartesian coordinates on the unit sphere
    """
    # convert to radians
    lat = math.radians(lat)
    lon = math.radians(lon)
    # calculations
    x = math.cos(lon) * math.cos(lat)
    y = math.sin(lon) * math.cos(lat)
    z = math.sin(lat)
    cart = (x, y, z)

    return cart





def cart_2_lat_lon(cart):
    """
    Convert cartesian coordinates to lat/lon.

    inputs:
    - cart = (x, y, z) cartesian coordinates on the unit sphere

    outputs:
    - lat = latitude (-90 to 90)
    - lon = longitude (-180 to 180)
    """
    # calculations
    lon = math.atan2(cart[1], cart[0])
    lat = math.atan2(cart[2], math.sqrt(cart[0]**2 + cart[1]**2))
    # convert to degrees
    lat = math.degrees(lat)
    lon = math.degrees(lon)

    return lat, lon





def fast_cross(a, b):
    """
    3D matrix cross multiplication.
    source: http://ssb.stsci.edu/doc/stsci_python_x/stsci.sphere.doc/html/_modules/stsci/sphere/great_circle_arc.html

    inputs:
    - a, b = matrices to be cross multiplied

    outputs:
    - cp = cross product
    """
    cp = np.empty(np.broadcast(a, b).shape)
    aT = a.T
    bT = b.T
    cpT = cp.T
    cpT[0] = aT[1]*bT[2] - aT[2]*bT[1]
    cpT[1] = aT[2]*bT[0] - aT[0]*bT[2]
    cpT[2] = aT[0]*bT[1] - aT[1]*bT[0]

    return cp





def cross_and_normalize(A, B):
    """
    3D matrix cross multiplication and normalized.
    source: http://ssb.stsci.edu/doc/stsci_python_x/stsci.sphere.doc/html/_modules/stsci/sphere/great_circle_arc.html

    inputs:
    - A, B = matrices to be cross multiplied

    outputs:
    - TN = normalized cross product
    """
    T = fast_cross(A, B)
    # normalization
    l = np.sqrt(np.sum(T ** 2, axis=-1))
    l = np.expand_dims(l, 2)
    # might get some divide-by-zeros, but we don't care
    with np.errstate(invalid='ignore'):
        TN = T / l

    return TN





def intersection(A, B, C, D):
    """
    Point of intersection between two great circle arcs.
    source: http://ssb.stsci.edu/doc/stsci_python_x/stsci.sphere.doc/html/_modules/stsci/sphere/great_circle_arc.html

    inputs:
    - A, B = (*x*, *y*, *z*) triples or Nx3 arrays of triples. Endpoints of the first great circle arc.
    - C, D = (*x*, *y*, *z*) triples or Nx3 arrays of triples. Endpoints of the second great circle arc.

    outputs:
    - T = (*x*, *y*, *z*) triples or Nx3 arrays of triples, If the given arcs intersect,
          the intersection is returned.  If the arcs do not intersect, the triple is set to all NaNs.
    """
    A = np.asanyarray(A)
    B = np.asanyarray(B)
    C = np.asanyarray(C)
    D = np.asanyarray(D)

    A, B = np.broadcast_arrays(A, B)
    C, D = np.broadcast_arrays(C, D)

    ABX = fast_cross(A, B)
    CDX = fast_cross(C, D)
    T = cross_and_normalize(ABX, CDX)
    T_ndim = len(T.shape)

    if T_ndim > 1:
        s = np.zeros(T.shape[0])
    else:
        s = np.zeros(1)
    s += np.sign(inner1d(fast_cross(ABX, A), T))
    s += np.sign(inner1d(fast_cross(B, ABX), T))
    s += np.sign(inner1d(fast_cross(CDX, C), T))
    s += np.sign(inner1d(fast_cross(D, CDX), T))
    if T_ndim > 1:
        s = np.expand_dims(s, 2)

    cross = np.where(s == -4, -T, np.where(s == 4, T, np.nan))

    # If they share a common point, it's not an intersection.  This
    # gets around some rounding-error/numerical problems with the
    # above.
    equals = (np.all(A == C, axis=-1) |
              np.all(A == D, axis=-1) |
              np.all(B == C, axis=-1) |
              np.all(B == D, axis=-1))

    equals = np.expand_dims(equals, 2)

    T = np.where(equals, np.nan, cross)

    return T





def plot_reconstruction(reconstructed_feature_geometries_list, color_list, lon_0):
    """
    Plot a global reconstruction from pygplates.

    Parameters
    ----------
    reconstructed_feature_geometries_list : list of lists
        list of lists of reconstructed features
    color_list : list of colors
        list of matplotlib color for geometries
    lon_0 : float
        the central longitude for viewing

    Returns
    -------
    None.
    """
    # initialize map
    fig = plt.figure(figsize=(12,10))

    ax = plt.subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=lon_0))

    ax.set_title(str(reconstructed_feature_geometries_list[0][0].get_reconstruction_time()) + ' Ma')
    ax.gridlines(xlocs=np.arange(-180,181,60),ylocs=np.arange(-90,91,30),linestyle='--')

    # loop over each reconstructed geometry list
    for i in range(len(reconstructed_feature_geometries_list)):

        # loop over each reconstructed geometry
        for j in range(len(reconstructed_feature_geometries_list[i])):

            # pull out lat/lon vertices
            lat_lon_array = reconstructed_feature_geometries_list[i][j].get_reconstructed_geometry().to_lat_lon_array()
            lats = lat_lon_array[:,0]
            lons = lat_lon_array[:,1]

            # wrapping to fix incorrect plotting across dateline
            for k in range(len(lons)):
                if lons[k]<0:
                    lons[k] = 180 + 180 + lons[k]

            # zip the result
            poly = Polygon(zip(lons, lats))

            # add the polygon to the map
            ax.add_geometries([poly], ccrs.PlateCarree(), facecolor=color_list[i], edgecolor='k', alpha=0.5)
            #ax.plot(lons,lats,transform=ccrs.PlateCarree(),color=color_list[i],linewidth=2)

    plt.show()





def check_polygon_in_band(polygon, lat_min, lat_max):
    """
    Check to see whether any part of a given polygon is inside a given latitude band.

    inputs:
    - polygon = pygpplates polygon
    - lat_min = the minimum latitude of the latitude band
    - lat_max = the maximum latitude of the latitude band

    ouputs:
    - in_band = boolean - True if inside, False if outside
    """
    # pull out lat/lon vertices
    lat_lon_array = polygon.to_lat_lon_array()
    lats = lat_lon_array[:,0]

    # check to see if any part of the polygon falls into our latitude band
    in_band = False
    for j in range(len(lats)):
        if lats[j]>lat_min and lats[j]<lat_max:
            in_band = True
            break

    return in_band





def get_area_in_band(polygon, lat_min, lat_max):
    """
    Calculate the area of a given polygon inside a given latitude band.

    inputs:
    - polygon = pygpplates polygon
    - lat_min = the minimum latitude of the latitude band
    - lat_max = the maximum latitude of the latitude band

    outputs:
    - area = the area of the polygon within the latitude band (in km^2)
    - band_polygon = pygplates polygon, with the parts outside of the latitude band removed
    """
    # pull out lat/lon vertices
    lat_lon_array = polygon.to_lat_lon_array()
    lats = lat_lon_array[:,0]
    lons = lat_lon_array[:,1]

    # storage lists
    bookmarks = []
    lat_add_list = []
    lon_add_list = []

    # iterate through the points
    for i in range(1,len(lats)):
        top_cross = False
        bot_cross = False

        # case where we move into the band from above
        if lats[i-1]>lat_max and lats[i]<lat_max:
            top_cross = True
        # case where we move out of the band from below
        if lats[i-1]<lat_max and lats[i]>lat_max:
            top_cross = True
        # case where we move out of the band from above
        if lats[i-1]>lat_min and lats[i]<lat_min:
            bot_cross = True
        # case where we move into the band from below
        if lats[i-1]<lat_min and lats[i]>lat_min:
            bot_cross = True

        # do calculations if we cross
        if top_cross or bot_cross:

            # convert the endpoints of the polygon segment into cartesian
            A = lat_lon_2_cart(lats[i-1], lons[i-1])
            B = lat_lon_2_cart(lats[i]  , lons[i])

            # get the intersection point (for the top and bottom cases), and convert back to lat/lon
            if top_cross:
                C_top = lat_lon_2_cart(lat_max, min([lons[i],lons[i-1]]))
                D_top = lat_lon_2_cart(lat_max, max([lons[i],lons[i-1]]))
                T = intersection(A, B, C_top, D_top)
            else:
                C_bot = lat_lon_2_cart(lat_min, min([lons[i],lons[i-1]]))
                D_bot = lat_lon_2_cart(lat_min, max([lons[i],lons[i-1]]))
                T = intersection(A, B, C_bot, D_bot)
            lat_add, lon_add = cart_2_lat_lon(T)

            # add to the storage lists
            bookmarks.append(i)
            lat_add_list.append(lat_add)
            lon_add_list.append(lon_add)

    # now insert the stored values into the original arrays
    new_lats = np.insert(lats, bookmarks, lat_add_list)
    new_lons = np.insert(lons, bookmarks, lon_add_list)

    # only keep values below the maximum latitude (with small buffer)
    mask = np.less(new_lats, lat_max+0.1)
    new_lats = new_lats[mask]
    new_lons = new_lons[mask]

    # only keep values above the minimum latitude
    mask = np.greater(new_lats, lat_min-0.1)
    new_lats = new_lats[mask]
    new_lons = new_lons[mask]

    # create a Polygon, if we are left with enough points
    if len(new_lats) >= 3:
        band_polygon = pygplates.PolygonOnSphere(zip(new_lats,new_lons))

        # get the area in km2
        area = band_polygon.get_area() * 6371.009**2

    # if we don't...
    else:
        area = 0
        band_polygon = None

    return area, band_polygon





def plot_polygons(polygon_list, facecolor, lon_0):
    """
    Plot pygplates polygons.

    inputs:
    - polygon_list = list of pygplates polygons
    - facecolor = string of matplotlib colour for geometries
    - lon_0 = the central longitude for viewing

    outputs:
    - fig = figure handle
    - ax_map = axis handle
    """
    # initialize map
    fig = plt.figure(figsize=(12,10))

    ax = plt.subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=lon_0))
    ax.gridlines(xlocs=np.arange(-180,181,60),ylocs=np.arange(-90,91,30),linestyle='--')

    # loop over each polygon
    for i in range(len(polygon_list)):

        if polygon_list[i] != None:
            # pull out lat/lon vertices
            lat_lon_array = polygon_list[i].to_lat_lon_array()
            lats = lat_lon_array[:,0]
            lons = lat_lon_array[:,1]

            # wrapping to fix incorrect plotting across dateline
            for j in range(len(lons)):
                if lons[j]<0:
                    lons[j] = 180 + 180 + lons[j]

            # zip the result
            poly = Polygon(zip(lons, lats))

            # add the polygon to the map
            ax.add_geometries([poly], ccrs.PlateCarree(), facecolor=facecolor, edgecolor='k', alpha=0.5)
            #ax.plot(lons,lats,transform=ccrs.PlateCarree(),color=color,linewidth=2)

    return fig, ax
    plt.show()





def get_areas_in_bands(reconstructed_feature_geometries, lat_mins, lat_maxs):
    """
    Get the area of all features in each latitude band.

    inputs:
    - reconstructed_feature_geometries = list of reconstructed features
    - lat_mins = list/np array of latitude minimums
    - lat_maxs = list/np array of latitude maximums

    outputs:
    - areas = list of total area in each latitude band
    - area_polygons = all polygons for which areas were calculated
    """
    # storage vectors
    areas = []
    area_polygons = []

    # iterate over each latitude band
    for i in range(len(lat_mins)):

        accumulated_area = 0

        # iterate over each polygon
        for j in range(len(reconstructed_feature_geometries)):

            current_polygon = reconstructed_feature_geometries[j].get_reconstructed_geometry()

            # check if the polygon is in the band
            in_band = check_polygon_in_band(current_polygon, lat_mins[i], lat_maxs[i])

            if in_band:
                # do the calculation
                area, band_polygon = get_area_in_band(current_polygon, lat_mins[i], lat_maxs[i])

                # store results
                accumulated_area = accumulated_area + area
                area_polygons.append(band_polygon)

        # store total area for the band
        areas.append(accumulated_area)

    return areas, area_polygons





def initialize_LIP_dict(LIP_feature_collection):
    """
    Initialize the dictionary which contains the LIP fraction remaining for all LIPs.

    inputs:
    - LIP_feature_collection = feature collection of LIPs

    outputs:
    - LIP_fracs = dictionary with keys = LIP names, values = LIP fraction remaining
    """
    # get the unique ID associated with each LIP geometry
    LIP_Ids = []
    for feature_collection in LIP_feature_collection:
        for feature in feature_collection:
            LIP_Id = feature.get_feature_id().get_string()
            LIP_Ids.append(LIP_Id)

    # create a dictionary: key = LIP Id, value = LIP fraction remaining
    ones = [1]*len(LIP_Ids)
    LIP_fracs = dict(zip(LIP_Ids, ones))

    return LIP_fracs





def weather_LIPs(reconstructed_feature_geometries, LIP_fracs, Ts, Rs, t_step, thickness, density):
    """
    A model for weathering LIPs.

    inputs:
    - reconstructed_feature_geometries = list of reconstructed features
    - LIP_fracs = dictionary with keys = LIP Ids, values = LIP fraction remaining
    - Ts = dataframe of zonal temperature
    - Rs = dataframe of zonal runoff
    - t_step = time since last weathering calculation (million years)
    - thickness = LIP thickness (km)
    - density = LIP density (kg/km^3)

    outputs:
    - modifies LIP_fracs in place
    """
    # only calculate if there are geometries
    if len(reconstructed_feature_geometries) > 0:
        # the time
        t = int(reconstructed_feature_geometries[0].get_reconstruction_time())

        # only calculate if the t is in our climate model
        T_string = 'T_C_'+str(t)
        if T_string in Ts.columns:

            # iterate over each polygon
            for i in range(len(reconstructed_feature_geometries)):

                # get the Id
                current_Id = reconstructed_feature_geometries[i].get_feature().get_feature_id().get_string()

                # only calculate if there is some LIP left
                if LIP_fracs[current_Id] > 0:

                    # get the polygon, ID, area, mass, and centre latitude
                    current_polygon = reconstructed_feature_geometries[i].get_reconstructed_geometry()

                    og_mass = current_polygon.get_area() * (6371.009**2) * thickness * density

                    current_area = (current_polygon.get_area() * 6371.009**2) * LIP_fracs[current_Id]
                    current_mass = current_area * thickness * density

                    centre_point = current_polygon.get_interior_centroid()
                    centre_lat = centre_point.to_lat_lon_array()[0][0]

                    # get the T and R at this time and latitude
                    for j in range(len(Ts.index)):
                        if centre_lat<Ts['lat_maxs'][j] and centre_lat>Ts['lat_mins'][j]:
                            T = Ts['T_C_'+str(t)][j]
                    for j in range(len(Rs.index)):
                        if centre_lat<Rs['lat_maxs'][j] and centre_lat>Rs['lat_mins'][j]:
                            R = Rs['R_m/yr_'+str(t)][j]

                    # get the weathering rate (t/km^2/yr)
                    fw = (R * 18.41 * np.exp(0.0533 * T)) / 0.649

                    # weather (kg)
                    weather_mass = fw * current_area * abs(t_step)*1e6 * 1000

                    # update the mass
                    current_mass = current_mass - weather_mass

                    # update the remaining fraction in place
                    if current_mass > 0:
                        LIP_fracs[current_Id] = current_mass / og_mass
                    else:
                        LIP_fracs[current_Id] = 0

        # if the t is not in our climate model, fill the dictionary with nans
        else:
            # if the dictionary is already filled with nans, we can move on
            if math.isnan(LIP_fracs[reconstructed_feature_geometries[0].get_feature().get_feature_id().get_string()]):
                pass
            else:
                for key, value in LIP_fracs.iteritems():
                    LIP_fracs[key] = float('NaN')





def get_LIP_areas_in_bands(reconstructed_feature_geometries, lat_mins, lat_maxs, LIP_fracs):
    """
    Get the area of all LIP features in each latitude band, accounting for erosion.

    inputs:
    - reconstructed_feature_geometries = list of reconstructed features
    - lat_mins = list/np array of latitude minimums
    - lat_maxs = list/np array of latitude maximums
    - LIP_fracs = dictionary with keys = LIP Ids, values = LIP fraction remaining

    outputs:
    - areas = list of total area in each latitude band (accounting for weathering)
    - area_polygons = all polygons for which areas were calculated (not accounting for weathering)
    """
    # only calculate if there are geometries
    if len(reconstructed_feature_geometries) > 0:
        # if our LIP_fracs has nan's
        if math.isnan(LIP_fracs[reconstructed_feature_geometries[0].get_feature().get_feature_id().get_string()]):
            areas = [np.nan]*len(lat_mins)
            area_polygons = None

        # if our LIP_fracs has values
        else:
            # storage vectors
            areas = []
            area_polygons = []

            # iterate over each latitude band
            for i in range(len(lat_mins)):

                accumulated_area = 0

                # iterate over each polygon
                for j in range(len(reconstructed_feature_geometries)):

                    current_polygon = reconstructed_feature_geometries[j].get_reconstructed_geometry()
                    current_Id = reconstructed_feature_geometries[j].get_feature().get_feature_id().get_string()

                    # check if the polygon is in the band
                    in_band = check_polygon_in_band(current_polygon, lat_mins[i], lat_maxs[i])

                    if in_band:
                        # do the calculation
                        area, band_polygon = get_area_in_band(current_polygon, lat_mins[i], lat_maxs[i])

                        # scale the area appropriately
                        area = area * LIP_fracs[current_Id]

                        # store results
                        accumulated_area = accumulated_area + area
                        area_polygons.append(band_polygon)

                # store total area for the band
                areas.append(accumulated_area)

    else:
        areas = [0]*len(lat_mins)
        area_polygons = []

    return areas, area_polygons





def get_LIP_areas_in_bands2(reconstructed_feature_geometries, lat_mins, lat_maxs, thresh=None, halflife=None, cover_thresh=None, covered_LIP_names=None):
    """
    Get the area of all features in each latitude band, with optional calculations for:
    - features instantly disappearing after `thresh` years
    - features decaying exponentially
    - covered features instantly disappearing after `cover_thresh` years, and other features
      decaying exponentially

    Parameters
    ----------
    reconstructed_feature_geometries : list
        list of reconstructed features
    lat_mins : array-like
        array-like of latitude minimums
    lat_maxs : array_like
        array_like of latitude maximums
    thresh : array-like (optional)
        after this many Myrs, the feature will instantly disappear
    halflife : array-like (optional)
        half life of exponential decay
    cover_thresh : array-like
        after this many Myrs, the covered feature will instantly disappear

    Returns
    -------
    areas : array
        list of total area in each latitude band
    area_polygons : list
        list of all polygons for which areas were calculated

    Optional Returns
    ----------------
    areas_thresh : list of arrays
        list of total area in each latitude band, accounting for thresholding method
    areas_decay : list of arrays
        list of total area in each latitude band, accounting for exponential decay method
    areas_cover : list of arrays
        list of total area in each latitude band, accounting for cover method
    """
    # storage vectors
    areas = np.array([])
    area_polygons = []

    # optional storage vectors
    if thresh != None:
        areas_thresh_temp = []
    if halflife != None:
        areas_decay_temp = []
        # convert halflife to decay constant
        lamb = np.log(2)/halflife
    if cover_thresh != None:
        areas_cover_temp = []

    # iterate over each latitude band
    for i in range(len(lat_mins)):

        accumulated_area = 0

        if thresh != None:
            accumulated_area_thresh = np.zeros(len(thresh))
        if halflife != None:
            accumulated_area_decay = np.zeros(len(halflife))
        if cover_thresh != None:
            accumulated_area_cover = np.zeros(len(cover_thresh))

        # iterate over each polygon
        for j in range(len(reconstructed_feature_geometries)):

            if thresh != None or halflife != None or cover_thresh != None:
                # get the begin date, reconstruction date, and feature age
                begin_date, end_date = reconstructed_feature_geometries[j].get_feature().get_valid_time()
                now_date = reconstructed_feature_geometries[j].get_reconstruction_time()
                feature_age = begin_date - now_date
                feature_name = reconstructed_feature_geometries[j].get_feature().get_name()

            # get the actual polygon
            current_polygon = reconstructed_feature_geometries[j].get_reconstructed_geometry()

            # check if the polygon is in the band
            in_band = check_polygon_in_band(current_polygon, lat_mins[i], lat_maxs[i])

            if in_band:
                # do the calculation
                area, band_polygon = get_area_in_band(current_polygon, lat_mins[i], lat_maxs[i])

                # store results
                accumulated_area = accumulated_area + area
                area_polygons.append(band_polygon)

                if thresh != None:
                    for k in range(len(thresh)):
                        # check if the polygon is past its due by date
                        if begin_date < (now_date+thresh[k]):
                            still_alive = True
                        else:
                            still_alive = False

                        if still_alive:
                            accumulated_area_thresh[k] = accumulated_area_thresh[k] + area

                if halflife != None:
                    for k in range(len(halflife)):
                        # scale the area based on the exponential decay equation
                        decay_area = area * np.exp(-lamb[k]*feature_age)

                        # store results
                        accumulated_area_decay[k] = accumulated_area_decay[k] + decay_area

                if cover_thresh != None:
                    if feature_name in covered_LIP_names:
                        for k in range(len(cover_thresh)):
                            # check if the polygon is past its due by date
                            if begin_date < (now_date+cover_thresh[k]):
                                still_alive = True
                            else:
                                still_alive = False

                            if still_alive:
                                accumulated_area_cover[k] = accumulated_area_cover[k] + area
                    else:
                        for k in range(len(cover_thresh)):
                            # scale the area based on the exponential decay equation
                            decay_area = area * np.exp(-lamb[k]*feature_age)

                            # store results
                            accumulated_area_cover[k] = accumulated_area_cover[k] + decay_area

        # store total area for the band
        areas = np.append(areas, accumulated_area)

        if thresh != None:
            areas_thresh_temp.append(accumulated_area_thresh)
        if halflife != None:
            areas_decay_temp.append(accumulated_area_decay)
        if cover_thresh != None:
            areas_cover_temp.append(accumulated_area_cover)

    # basically, flip our outputs so that each array in our list is for a given input thresh/halflife
    if thresh != None:
        areas_thresh = []
        for i in range(len(thresh)):
            this_array = np.array([])
            for j in range(len(areas_thresh_temp)):
                this_array = np.append(this_array, areas_thresh_temp[j][i])
            areas_thresh.append(this_array)
    if halflife != None:
        areas_decay = []
        for i in range(len(halflife)):
            this_array = np.array([])
            for j in range(len(areas_decay_temp)):
                this_array = np.append(this_array, areas_decay_temp[j][i])
            areas_decay.append(this_array)
    if cover_thresh != None:
        areas_cover = []
        for i in range(len(cover_thresh)):
            this_array = np.array([])
            for j in range(len(areas_cover_temp)):
                this_array = np.append(this_array, areas_cover_temp[j][i])
            areas_cover.append(this_array)

    # returns
    if thresh == None and halflife == None and cover_thresh == None:
        return areas, area_polygons
    elif thresh != None and halflife == None and cover_thresh == None:
        return areas, area_polygons, areas_thresh
    elif thresh != None and halflife != None and cover_thresh == None:
        return areas, area_polygons, areas_thresh, areas_decay
    elif thresh != None and halflife == None and cover_thresh != None:
        return areas, area_polygons, areas_thresh, areas_cover
    elif thresh != None and halflife != None and cover_thresh != None:
        return areas, area_polygons, areas_thresh, areas_decay, areas_cover
    elif thresh == None and halflife != None and cover_thresh == None:
        return areas, area_polygons, areas_decay
    elif thresh == None and halflife != None and cover_thresh != None:
        return areas, area_polygons, areas_decay, areas_cover
    elif thresh == None and halflife == None and cover_thresh != None:
        return areas, area_polygons, areas_cover