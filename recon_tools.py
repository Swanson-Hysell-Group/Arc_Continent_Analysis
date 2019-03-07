import math
from numpy.core.umath_tests import inner1d
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import numpy as np
from shapely.geometry.polygon import Polygon
import pygplates
import pandas as pd


########## BASIC CALCULATIONS ##########


def lat_lon_2_cart(lat, lon):
    """
    Convert lat/lon coordinates to cartesian.

    Parameters
    ----------
    lat : array-like
        latitude (-90 to 90)
    lon : array-like
        longitude (-180 to 180)

    Returns
    -------
    cart : tuple
        (x, y, z) cartesian coordinates on the unit sphere
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

    Parameters
    ----------
    cart : tuple
        (x, y, z) cartesian coordinates on the unit sphere

    Returns
    -------
    lat : array-like
        latitude (-90 to 90)
    lon : array-like
        longitude (-180 to 180)
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

    Parameters
    ----------
    a, b : numpy arrays
        matrices to be cross multiplied

    Returns
    -------
    cp : numpy array
        cross product
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

    Parameters
    ----------
    A, B : numpy arrays
        matrices to be cross multiplied

    Returns
    -------
    TN : numpy array
        normalized cross product
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

    Parameters
    ----------
    A, B : (*x*, *y*, *z*) triples or Nx3 arrays of triples
        Endpoints of the first great circle arc.
    C, D : (*x*, *y*, *z*) triples or Nx3 arrays of triples
        Endpoints of the second great circle arc.

    Returns
    -------
    T : (*x*, *y*, *z*) triples or Nx3 arrays of triples
        If the given arcs intersect, the intersection is returned. If the arcs do not intersect,
        the triple is set to all NaNs.
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


########## ZONAL CALCULATIONS ##########


def check_polygon_in_band(polygon, lat_min, lat_max):
    """
    Check to see whether any part of a given polygon is inside a given latitude band.

    Parameters
    ----------
    polygon : pygpplates polygon
    lat_min : float
        the minimum latitude of the latitude band
    lat_max : float
        the maximum latitude of the latitude band

    Returns
    -------
    in_band : boolean
        True if inside, False if outside
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

    Parameters
    ----------
    polygon : pygpplates polygon
    lat_min : float
        the minimum latitude of the latitude band
    lat_max : float
        the maximum latitude of the latitude band

    Returns
    -------
    area : float
        the area of the polygon within the latitude band (in km^2)
    band_polygon : pygplates polygon
        with the parts outside of the latitude band removed
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


def get_areas_in_bands(reconstructed_feature_geometries, lat_mins, lat_maxs):
    """
    Get the area of all features in each latitude band.

    Parameters
    ----------
    reconstructed_feature_geometries : list
        list of reconstructed features
    lat_mins : array-like
        array-like of latitude minimums
    lat_maxs : array_like
        array_like of latitude maximums

    Returns
    -------
    areas : array
        list of total area in each latitude band
    area_polygons : list
        list of all polygons for which areas were calculated
    """
    # storage vectors
    areas = np.array([])
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
        areas = np.append(areas, accumulated_area)

    return areas, area_polygons


def get_length_in_band(polyline, lat_min, lat_max):
    """
    Calculate the length of a given polyline inside a given latitude band.

    Parameters
    ----------
    polyline : pygpplates polyline
    lat_min : float
        the minimum latitude of the latitude band
    lat_max : float
        the maximum latitude of the latitude band

    Returns
    -------
    length : float
        the length of the polyline within the latitude band (in km)
    band_polylines : list
        pygplates polylines with the parts outside of the latitude band removed
    """
    # pull out lat/lon vertices
    lat_lon_array = polyline.to_lat_lon_array()
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

    # mask points above and below the latitude band (with small buffer)
    top_mask = np.less(new_lats, lat_max+0.1)
    bot_mask = np.greater(new_lats, lat_min-0.1)
    mask = top_mask & bot_mask

    # initiate final output variabiles
    length = 0
    band_polylines = []

    # initiate lat lon arrays for the current polyline that is being split out
    split_lats = np.array([])
    split_lons = np.array([])
    hold = False

    # iterate through the mask
    for i in range(len(mask)):
        if mask[i]:
            # begin holding, and store points
            hold = True
            split_lats = np.append(split_lats, new_lats[i])
            split_lons = np.append(split_lons, new_lons[i])

        # i.e. the first time we reach a point outside of the band
        elif hold:
            # create a polyline, if we are left with enough points
            if len(split_lats) >= 2:
                band_polyline = pygplates.PolylineOnSphere(zip(split_lats,split_lons))

                # get the length in km
                length_split = band_polyline.get_arc_length() * pygplates.Earth.mean_radius_in_kms

                # store
                length = length + length_split
                band_polylines.append(band_polyline)

            # reset stuff
            split_lats = np.array([])
            split_lons = np.array([])
            hold = False

    # if the last point was within the band
    if hold:
        # create a polyline, if we are left with enough points
        if len(split_lats) >= 2:
            band_polyline = pygplates.PolylineOnSphere(zip(split_lats,split_lons))

            # get the length in km
            length_split = band_polyline.get_arc_length() * pygplates.Earth.mean_radius_in_kms

            # store
            length = length + length_split
            band_polylines.append(band_polyline)

    return length, band_polylines


def get_lengths_in_bands(reconstructed_feature_geometries, lat_mins, lat_maxs):
    """
    Get the area of all features in each latitude band.

    Parameters
    ----------
    reconstructed_feature_geometries : list
        list of reconstructed features
    lat_mins : array-like
        array-like of latitude minimums
    lat_maxs : array_like
        array_like of latitude maximums

    Returns
    -------
    accumulated_lengths : array
        array of total area in each latitude band
    length_polylines : polylines
        list of all polylines for which areas were calculated
    polyline_attributes : Pandas Dataframe
        Dataframe with names, lengths, latitudes and recon_times for each polyline
    """
    # storage vectors
    lengths = np.array([])
    accumulated_lengths = np.array([])
    length_polylines = []
    names = []
    polyline_lat_min = []
    polyline_lat_max = []
    recon_times = []

    # iterate over each latitude band
    for i in range(len(lat_mins)):

        accumulated_length = 0

        # iterate over each polygon
        for j in range(len(reconstructed_feature_geometries)):

            current_polyline = reconstructed_feature_geometries[j].get_reconstructed_geometry()


            # check if the polygon is in the band
            in_band = check_polygon_in_band(current_polyline, lat_mins[i], lat_maxs[i])

            if in_band:
                # do the calculation
                length, band_polyline = get_length_in_band(current_polyline, lat_mins[i], lat_maxs[i])
                name = reconstructed_feature_geometries[j].get_feature().get_name()
                recon_time = reconstructed_feature_geometries[j].get_reconstruction_time()

                # store results
                accumulated_length = accumulated_length + length
                length_polylines.append(band_polyline)
                names.append(name)
                lengths = np.append(lengths, length)
                polyline_lat_min.append(lat_mins[i])
                polyline_lat_max.append(lat_maxs[i])
                recon_times.append(recon_time)

        # store total area for the band
        accumulated_lengths = np.append(accumulated_lengths, accumulated_length)

        polyline_attributes = pd.DataFrame({'name':names,'length':lengths,'polyline_lat_min':polyline_lat_min,'polyline_lat_max':polyline_lat_max,'recon_time':recon_times})

    return accumulated_lengths, length_polylines, polyline_attributes

########## PLOTTING FUNCTIONS ##########

def plot_feature(ax,feature,color='grey',linewidth=1):
    for n in range(len(feature)):

        # pull out lat/lon vertices
        lat_lon_array = feature[n].to_lat_lon_array()
        lats = lat_lon_array[:,0]
        lons = lat_lon_array[:,1]

        ax.plot(lons,lats, transform=ccrs.Geodetic(), color=color, linewidth=linewidth)


def plot_reconstructed_feature(ax,reconstructed_feature,color='grey',linewidth=1):
    for n in range(len(reconstructed_feature)):

        # pull out lat/lon vertices
        lat_lon_array = reconstructed_feature[n].get_reconstructed_geometry().to_lat_lon_array()
        lats = lat_lon_array[:,0]
        lons = lat_lon_array[:,1]

        ax.plot(lons,lats, transform=ccrs.Geodetic(), color=color, linewidth=linewidth)


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

            # zip the result
            poly = Polygon(zip(lons, lats))

            # add the polygon to the map
            ax.add_geometries([poly], ccrs.PlateCarree(), facecolor=color_list[i], edgecolor='k', alpha=0.5)
            #ax.plot(lons,lats,transform=ccrs.Geodetic(),color=color_list[i],linewidth=1)

    plt.show()


def plot_polygons(polygon_list, color, lon_0):
    """
    Plot pygplates polygons.

    Parameters
    ----------
    polygon_list : list
        list of pygplates polygons
    color : color
        matplotlib color for geometries
    lon_0 : float
        the central longitude for viewing

    Returns
    -------
    None.
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

            # zip the result
            poly = Polygon(zip(lons, lats))

            # add the polygon to the map
            ax.add_geometries([poly], ccrs.PlateCarree(), facecolor=color, edgecolor='k', alpha=0.5)
            #ax.plot(lons,lats,transform=ccrs.Geodetic(),color=color,linewidth=1)

    plt.show()
