#! /usr/bin/env python
# ----------------------------------------------------------
# Author: 				Jesse Thompson, Environment and Climate Change Canada
# Date:   				2018-01-17 16:07:59
# Email: 				jessethompsonj@gmail.com
# Modified by:   	Jesse Thompson, Environment and Climate Change Canada
# Modified time:	2018-04-27 15:43:32
# Mark Shephard (ECCC) : 2018-06-11
#  - fixed bug going from negative to postive lon range
# Enrico Dammers(ECCC) : 2018-08-28
# improved multi thread functionality for multi month/year big regions

from __future__ import division, print_function
import os
import sys
import argparse
import matplotlib as mpl
import datetime
import numpy as np


if "DISPLAY" not in os.environ:
    mpl.use('Agg')
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# CONSTANTS AND OTHER GLOBAL VARS/FLAGS
# ----------------------------------------------------------------------------

# =============================================================================
# Defining constants
# =============================================================================

G = 6.674e-11  # Gravitational Constant
c = 2.998e8  # Speed of light
h = 6.626e-34  # Planck's constant
hbar = h / (2 * np.pi)
k = 1.381e-23  # Boltzmann Constant
sigma = 5.670e-8  # Stefan-Boltzmann Constant
m_e = 9.109e-31  # Electron mass
M_s = 1.989e30  # Mass of Sun
R_s = 6.963e8  # Radius of Sun
L_s = 3.828e26  # Luminosity of Sun
m_p = 1.6726219e-27  # Mass of proton
a = 4 * sigma / c
pi = np.pi
ep_pp_coeff = 1.07e-7 * 1e-5 * (1e-6) ** 4
ep_cno_coeff = 8.24e-26 * 1e-5 * (1e-6) ** 19.9
nonrelgenpress = (3 * pi ** 2) ** (2 / 3) / 5 * hbar ** 2 / m_e * m_p ** (-5 / 3)
mach_ep = np.finfo(np.float64).eps
tiny_float = 1e-20
gamma = 5 / 3  # ideal gas constant

# =================Defining mu, XYZ=====================================
X = 0.7381
Y = 0.2485
Z = 1 - (X + Y)
mu = 1 / (2 * X + 0.75 * Y + 0.5 * Z)
# mu = mean molecular weight for fully ionized gas


# =============================================================

# DEFAULTS
binsize = 0.1
oversample_binsize = 0.5

# use_pickle = True
average_flag = False
world_flag = False

n_bins = 14

# SATELLITE CONSTANTS
cris_fovDia = 0.016808
fill_value = -999.5

# TEMPORARY FOLDER FILES
KMLDIR = "/home/mas001/projects/retv/cris/global_v1_5_beta/plots/KML_files_VIIRS-CrIS/"

currentdir = prt.formatDir(
    os.getcwd()
) if "cdf" not in os.getcwd() else prt.formatDir(os.getcwd()[:-3])

lower_lims = [  # MUST HAVE SAME LENGTH AS CBAR_COLORS
    0.,
#    .25,
    .5,
#    .75,
    1.0,
    1.25,
    1.5,
    1.75,
    2.0,
    2.5,
    3.,
    4.,
    5.,
    7.5,
    10.
#    15.
]

# (r,g,b) -> red(1,0,0), green(0,1,0), blue(0,0,1)
# 256 from 0 to 255 scaled between 0 and 1

cbar_colors = [
    (.9, .9, .9),  # LIGHT GREY (0-0.25)
    (.8, .8, .8),  # LIGHT GREY (0-0.25)
    (.5, .5, .5),  # GREY (0.25-0.5)
#    (0, .6, 1),  # BLUE (0.5-0.75)
#    (0.2, 1, 1),  # LIGHT BLUE (0.75 - 1)
#    (0.4, 1., 0.6),  # TEAL (1.0 - 1.5)
    (0, .5, 0),  # DARKGREEN (1.5 - 2.0)
    (0, .9, 0),  # GREEN (1.5 - 2.0)
    (0.6, 1, 0.),  # YELLOW GREEN (2.0 - 2.5)
    (1, 1, 0),  # YELLOW (2.5 - 3.0)
    (1., 0.8, 0.2),  # YELLOW ORANGE (3.0 - 4)
    (1, 0.6, 0),  # ORANGE (4.0 - 5.0)
    (0.8, 0.2, 0),  # ORANGE RED (5.0 - 7.5)
    (1, 0, 0),  # RED (7.5 - 10.0)
    (0.8, 0.2, 0.4),  # RED/FUSCHIA (10.0 - 15)
    (0.8, 0, 1)  # FUSCHIA (15 +)
]

#cbar_colors = [
#    (.8, .8, .8),  # LIGHT GREY (0-0.25)
#    (.5, .5, .5),  # GREY (0.25-0.5)
#    (0, .6, 1),  # BLUE (0.5-0.75)
#    (0.2, 1, 1),  # LIGHT BLUE (0.75 - 1)
#    (0.4, 1., 0.6),  # TEAL (1.0 - 1.5)
#    (0, 1, 0),  # GREEN (1.5 - 2.0)
#    (0.6, 1, 0.),  # YELLOW GREEN (2.0 - 2.5)
#    (1, 1, 0),  # YELLOW (2.5 - 3.0)
#    (1., 0.8, 0.2),  # YELLOW ORANGE (3.0 - 4)
#    (1, 0.6, 0),  # ORANGE (4.0 - 5.0)
#    (0.8, 0.2, 0),  # ORANGE RED (5.0 - 7.5)
#    (1, 0, 0),  # RED (7.5 - 10.0)
#    (0.8, 0.2, 0.4),  # RED/FUSCHIA (10.0 - 15)
#    (0.8, 0, 1)  # FUSCHIA (15 +)
#]
#
#
# plt.register_cmap(
#     'xretv',
#     mpl.colors.LinearSegmentedColormap.from_list(
#         'xretv',
#         cbar_colors,
#         N=n_bins
#     )
# )
#
# cmap = mpl.colors.LinearSegmentedColormap.from_list('xretv', [(0. / 125, 'grey'),
#                                                   (5. / 125, 'blue'),
#                                                   (10. / 125, 'lightblue'),
#                                                   (15. / 125, 'green'),
#                                                   (25. / 125, 'lightgreen'),
#                                                   (50. / 125, 'gold'),
#                                                   (75. / 125., 'orange'),
#                                                   (100 / 125., 'red'),
#                                                   (1., 'purple')])
plt.register_cmap(
    'xretv',
    mpl.colors.LinearSegmentedColormap.from_list(
        'xretv',
        cbar_colors
    )
)
# plt.register_cmap('xretv', cbar_colors)
cm = mpl.cm.get_cmap('xretv')


# ----------------------------------------------------------------------------
# FUNCTIONS FOR PLOTTING
# ----------------------------------------------------------------------------

def nonlinear_custom_color_xretv(val):
    '''
    In order to have a non-linear colormap, We must force the
    values into bins, and return their colour value instead
    '''
    import numpy as np
    norm = mpl.colors.Normalize(vmin=0, vmax=len(cbar_colors) - 1)
    vmax = lower_lims[-1]

    for ind, lowerbound in enumerate(lower_lims):
        if val <= 0 or np.isnan(val):
            return cm(norm(0))
        elif val >= vmax:
            return cm(norm(14))
        elif lowerbound < val < lower_lims[ind + 1]:
            return cm(norm(ind + 0.01))
        else:
            continue


def write_to_nc_file(currentdir, date, xretv_data, datakey, lat_range, lon_range, filetail=False, Day_Night_Flag="All"):
    '''In order to be able to create larger yearly plots,
     the averaging has to be done per set of monthly'''
    import netCDF4
    # Generate directory string
    write_dir = "{}{}".format(currentdir, "plots/")
    if not os.path.exists(write_dir):
        os.makedirs(write_dir)
    write_dir = "{}{}".format(currentdir, "plots/nc_dir/")
    # Make directory if it doesn't exist
    if not os.path.exists(write_dir):
        os.makedirs(write_dir)
    # create netCDF file
    # Get bounds from the polygons
    if not lat_range[0] and not lon_range[0]:
        west, east = [
            min(longitudes) if min(longitudes) >= -180 else -180,
            max(longitudes) if max(longitudes) <= 180 else 180
        ]
        south, north = [
            min(latitudes) if min(latitudes) >= -90 else -90,
            max(latitudes) if max(latitudes) <= 90 else 90
        ]
    else:
        west, east = [min(lon_range), max(lon_range)]
        south, north = [min(lat_range), max(lat_range)]
#    print('test2',datakey,len(datakey))
    if type(datakey) != list or len(datakey)==1:
        if len(datakey) == 1:
            datakey2 = datakey[0]
        savename = "{}_{}_{}_{}_{}_DNF_{}_{}_level3.nc".format(
            date,
            ot.regnum2strnum(round(west if not lon_range[0] else lon_range[0], 1)),
            ot.regnum2strnum(round(east if not lon_range[0] else lon_range[1], 1)),
            ot.regnum2strnum(round(south if not lat_range[0] else lat_range[0], 1)),
            ot.regnum2strnum(round(north if not lat_range[0] else lat_range[1], 1)),
            Day_Night_Flag.lower(),
            filetail if filetail else datakey2
        )
    else:
        savename = "{}_{}_{}_{}_{}_DNF_{}_{}_level3.nc".format(
            date,
            ot.regnum2strnum(round(west if not lon_range[0] else lon_range[0], 1)),
            ot.regnum2strnum(round(east if not lon_range[0] else lon_range[1], 1)),
            ot.regnum2strnum(round(south if not lat_range[0] else lat_range[0], 1)),
            ot.regnum2strnum(round(north if not lat_range[0] else lat_range[1], 1)),
            Day_Night_Flag.lower(),
            filetail if filetail else 'all'
        )
    # savename_txt = "{}_{}_{}_{}_{}_DNF_{}_{}_level3.txt".format(
    #     date,
    #     ot.regnum2strnum(round(west if not lon_range[0] else lon_range[0], 1)),
    #     ot.regnum2strnum(round(east if not lon_range[0] else lon_range[1], 1)),
    #     ot.regnum2strnum(round(south if not lat_range[0] else lat_range[0], 1)),
    #     ot.regnum2strnum(round(north if not lat_range[0] else lat_range[1], 1)),
    #     Day_Night_Flag.lower(),
    #     filetail if filetail else datakey
    # )
    list_of_vari = ['binsize', 'Longitude', 'Latitude', 'polyLon', 'polyLat', 'retv', 'retv_sdev', 'Filter_Mask', 'oversample_binsize',
                    'xretv_median', 'N_Oversampled_Points', 'N_Bin_Points', 'N_Dates', 'day_thresh_perc']

    naming_convention = ['description','units']
    expl_list = [
        ['Grid spacing in degrees latitude and longitude (e.g. 0.1 is a 0.1x0.1 degree grid) of each gridded variable','degrees'],
        ['Longitude of the centre of each grid box','degrees (Max range: -180W to +180E)'],
        ['Latitude of the centre of each grid box','degrees (Max range: +90N to -90 S)'],
        ['Longitude of the 4 corners of the each grid box', 'degrees (+ for East, and - for West)'],
        ['Latitude of the 4 corners of the each grid box','degrees (+ for North, and - for South)'],
        ['Gridded mean surface retrieved values over the time period specified in the filename','volume mixing ratio (ppmv)'],
        ['Standard deviation of the gridded surface retrieved values','volume mixing ratio (ppmv)'],
        ['Binary filter mask based on recommended spatial (wgt_thresh_perc) and temporal (day_thresh_perc) thresholds','0 : did not passed recommended filter test (these grid boxes should be filtered out) 1 : passed recommended filtering'],
        ['Specified oversampling distance used around each grid box in degrees (e.g. 0.5)','degrees'],
        ['Gridded median surface ammonia values over the time period specified in the filename', 'volume mixing ratio (ppmv)'],
        ['Number of observations within the grid including the oversampling region around it', 'unitless'],
        ['Number of observations within the grid not including the oversampling region around the grid.', 'untiless'],
        ['Number of unique days with an observation within the grid including the surrounding oversample region, For example, a monthly gridded product in June would have a maximum of 30 unique days.', 'unitless'],
        ['Fraction of the number of unique days compared to the maximum potential days with observations. For example, a specified threshold value of 0.2 would require for a monthly gridded product in June', 'unitless'],
        ['(with a maximum of 30 unique days) that there be at least  6 unique days within June that there were observations contributing to the gridded mean variable for the Mask_Filter=1.', 'unitless']]

    print(currentdir + "plots/nc_dir/" + savename)
    nc = netCDF4.Dataset(currentdir + "plots/nc_dir/" + savename, 'w')
    nc.createDimension('N', len(xretv_data['Longitude']))
    nc.createDimension('poly', 4)

    for key in xretv_data.keys() + []:
        if 'poly' in key:
            vari = nc.createVariable(key, float, ('N', 'poly'))
        else:
            vari = nc.createVariable(key, float, ('N'))
        vari[:] = xretv_data[key]
        if key in list_of_vari:
            key_index = list_of_vari.index(key)
            for idx, attr in enumerate(expl_list[key_index]):
                # print(key,key_index,idx,attr,naming_convention[idx])
                # print(expl_list[key_index])
                vari.setncattr(naming_convention[idx], attr)

    # list_of_lists = []
    # for i in list_of_vari:
    #     list_of_lists.append(((list_of_vari[i], expl_list[i]), ('Units', units_list[i])
    #                           print 'lala'
    #
    #                           if key in keylist:
    #                           key_id = keylist.index[key)
    #                          for attr in list_of_lists[key_id]:
    #     vari.setncattr(attr[0], attr[1])

    #obtain the filter values from the filetail
    print('filetail', filetail)
    split_filetail=filetail.split('_')
    #obtain the index before the value and then get the value
    ind_bin = split_filetail.index('bin')
    ind_os = split_filetail.index('os')
    ind_wgt = split_filetail.index('wgt')
    ind_dthres = split_filetail.index('dthres')
    binsize = split_filetail[ind_bin+1]
    oversample_binsize = split_filetail[ind_os+1]
    wgt_thresh_perc = split_filetail[ind_wgt+1]
    day_thresh_perc = split_filetail[ind_dthres+1]
    nc.setncattr('binsize', binsize)
    nc.setncattr('oversample_binsize', oversample_binsize)
    nc.setncattr('wgt_thresh_perc', wgt_thresh_perc)
    nc.setncattr('day_thresh_perc', day_thresh_perc)
    print("Filters:", binsize, oversample_binsize, wgt_thresh_perc, day_thresh_perc)


    # gbin = nc.createVariable('binsize', float)
    # os_bin = nc.createVariable('oversample_binsize', float)
    # wgt = nc.createVariable('wgt_thresh_perc', float)
    # day = nc.createVariable('day_thresh_perc', float)
    #
    # gbin[:] = binsize
    # os_bin[:] = oversample_binsize
    # wgt[:] = wgt_thresh_perc
    # day[:] = day_thresh_perc

    nc.close()
    # tmp = open(currentdir + "plots/nc_dir/" + savename_txt, 'w')
    # tmp.write('tmp\n')
    # tmp.close()
    return


def read_from_nc(currentdir, dates, datakey, lat_range, lon_range, filetail, Day_Night_Flag='All', wgt_thresh_perc=.05,
                 day_thresh_perc=.05):
    '''In order to be able to create larger yearly plots,
     the averaging has to be done per set of monthly
     here we read the monthly files'''
    import netCDF4
    import pandas as pd
    import numpy as np
    import glob
    # create netCDF file
    # Get bounds from the polygons
    if not lat_range[0] and not lon_range[0]:
        west, east = [
            min(longitudes) if min(longitudes) >= -180 else -180,
            max(longitudes) if max(longitudes) <= 180 else 180
        ]
        south, north = [
            min(latitudes) if min(latitudes) >= -90 else -90,
            max(latitudes) if max(latitudes) <= 90 else 90
        ]
    else:
        west, east = [min(lon_range), max(lon_range)]
        south, north = [min(lat_range), max(lat_range)]
    txtname = currentdir + "plots/nc_dir/" + "{}_{}_{}_{}_{}_{}_{}_final.txt".format(
        '_'.join(dates),
        ot.regnum2strnum(round(west if not lon_range[0] else lon_range[0], 1)),
        ot.regnum2strnum(round(east if not lon_range[0] else lon_range[1], 1)),
        ot.regnum2strnum(round(south if not lat_range[0] else lat_range[0], 1)),
        ot.regnum2strnum(round(north if not lat_range[0] else lat_range[1], 1)), Day_Night_Flag.lower(), filetail)
    # find list of files
    # files = glob.glob(currentdir + "plots/nc_dir/" + savename)
    read_txt = open(txtname)
    lines = read_txt.readlines()
    print(lines)
    # raise ValueError
    files = lines
    df = pd.DataFrame()
    for fil in files:
        # print(fil)
        df_tmp = pd.DataFrame()
        nc = netCDF4.Dataset(fil.strip('\n'))
        for key in nc.variables.keys():
            # print(key,nc.variables[key].shape)
            df_tmp[key] = [uu for uu in nc.variables[key][:]]
        print(fil,'\n',np.unique(df_tmp.Longitude))
        nc.close()
        df = df.append(df_tmp)
    # print(np.array([np.array(uu) for uu in df['polyLon'].values]).shape)
    df['lon_lat'] = ['%s_%s' % (u1, u2) for u1, u2 in zip(df['Longitude'], df['Latitude'])]
    if type(datakey) != list:
        df['%s_wgt' % datakey] = df[datakey] * df['weight']
    else:
        for datak in datakey:
            df['%s_wgt' % datak] = df[datak] * df['weight']
    df['polyLon1'] = [uu[0] for uu in df['polyLon']]
    df['polyLon2'] = [uu[1] for uu in df['polyLon']]
    df['polyLon3'] = [uu[2] for uu in df['polyLon']]
    df['polyLon4'] = [uu[3] for uu in df['polyLon']]
    df['polyLat1'] = [uu[0] for uu in df['polyLat']]
    df['polyLat2'] = [uu[1] for uu in df['polyLat']]
    df['polyLat3'] = [uu[2] for uu in df['polyLat']]
    df['polyLat4'] = [uu[3] for uu in df['polyLat']]
    grouped = df.groupby('lon_lat')
    grp2 = grouped.sum() / grouped.count()
    if type(datakey) != list:
        grp2[datakey] = grouped['%s_wgt' % datakey].sum() / grouped['weight'].sum()
    else:
        for datak in datakey:
            grp2[datak] = grouped['%s_wgt' % datak].sum() / grouped['weight'].sum()
    grp2['N'] = grouped.sum()['N']
    grp2['N_Dates'] = grouped.sum()['N_Dates']
    grp2['N_Bin_Points'] = grouped.sum()['N_Bin_Points']
    grp2['N_Oversampled_Points'] = grouped.sum()['N_Oversampled_Points']

    # back to the old dict
    xretv_out = {}
    for key in [key2 for key2 in grp2.keys() if key2 not in ['%s_wgt' % uu for uu in datakey]]:
        if 'poly' not in key: xretv_out[key] = np.array(grp2[key].values)
    xretv_out['polyLon'] = np.array([np.array([u1, u2, u3, u4]) for u1, u2, u3, u4 in
                                     zip(grp2['polyLon1'], grp2['polyLon2'], grp2['polyLon3'], grp2['polyLon4'])])
    xretv_out['polyLat'] = np.array([np.array([u1, u2, u3, u4]) for u1, u2, u3, u4 in
                                     zip(grp2['polyLat1'], grp2['polyLat2'], grp2['polyLat3'], grp2['polyLat4'])])
    for key in xretv_out.keys():
        print(xretv_out[key].shape, key)
    # add in wgt threshold / day threshold here
    # add in wgt threshold / day threshold here
    # for now the N_Dates etc are averaged, so you get a 12 month average in a yearset...
    # In the future maybe take it out of the standard summation/divide
    # also no filter is not a choice yet, if we do want that we need to include it in launch by parts
    print('keys', xretv_out.keys())
    print('N_dates, max min', max(xretv_out['N_Dates']), min(xretv_out['N_Dates']))
    print('Weight_Frac, max min', max(xretv_out['Weight_Frac']), min(xretv_out['Weight_Frac']))
    # if ((bindict['N_Dates'] < day_threshold) or
    #         (Weight_Frac < wgt_thresh_perc)):  # and no_filter:
    #     bindict[datakey] = -999.5

    # get filter parameters
    #max_weight = grouped['Max_weight'].sum()
    #weight = grouped['weight'].sum()
    Weight_Frac = grouped['Weight_Frac'].sum()
    Filter_Mask = grouped['Filter_Mask'].sum()
    N_dates = grouped['N_Dates'].sum()
    start = datetime.datetime.strptime(dates[0], "%Y%m%d").date()
    end = datetime.datetime.strptime(dates[1], "%Y%m%d").date() if len(dates) == 2 else (
            start + datetime.timedelta(days=0))
    num_days = (end - start).days
    day_threshold = round((day_thresh_perc * num_days), 0)
    # apply filtering
    print('N_dates', np.max(N_dates), np.min(N_dates))
    print('Weight_Frac', np.max(Weight_Frac), np.min(Weight_Frac))
    #for datak in datakey:
        #xretv_out[datak][((N_dates < day_threshold) |
        #                    (Weight_Frac < wgt_thresh_perc))] = -999.5
        #xretv_out[datak][(Filter_Mask < 1)] = -999.5
    # add in wgt threshold / day threshold here
    # add in wgt threshold / day threshold here
    # add in wgt threshold / day threshold here

    return xretv_out


def unpack_helper(args):
    return rpt.unpack_plotdata_monthly_yearly(*args)


def unpack_helper_single(args):
    return rpt.unpack_plotdata(*args)


def avg_bin_helper(args):
    return avg_by_bin(*args)


def launch_parts_jobs(cmd_list, dates, currentdir, lon_range, lat_range, filetail=None, Day_Night_Flag='All'):
    '''Split averaging into small parts and merge after'''
    import glob
    import time

    start = datetime.datetime.strptime(dates[0], "%Y%m%d").date()
    end = datetime.datetime.strptime(dates[1], "%Y%m%d").date() if len(dates) == 2 else (
            start + datetime.timedelta(days=0))
    nmonths = (end.year - start.year) * 12 + (end.month - start.month)

    month_list = [
        datetime.datetime(int(start.year + (start.month + nm - 1) / 12), int(np.mod(start.month - 1 + nm, 12) + 1),
                          start.day) for
        nm in range(int(nmonths + 1))]
    month_list_end = [datetime.datetime(int(start.year + (start.month + nm - 1) / 12),
                                        int(np.mod(start.month - 1 + nm, 12) + 1),
                                        start.day) - datetime.timedelta(1) for nm in range(nmonths + 1)][1:] + [end]
    base_cmd_tmp_b = [arg for arg in cmd_list if arg not in ['--parts']]

    # remove f-p filters from set
    idx_to_remove_tmp = [idx for idx, uu in enumerate(base_cmd_tmp_b) if 'thres' in uu]
    idx_to_remove = list(np.array(idx_to_remove_tmp)) + list(np.array(idx_to_remove_tmp) + 1)
    base_cmd_tmp = [arg for idx, arg in enumerate(base_cmd_tmp_b) if idx not in idx_to_remove]
    end_cmd_tmp = [arg for arg in cmd_list if arg not in ['--parts', '--to_nc']]

    # remove date from commands
    arg_index = base_cmd_tmp.index('-g')
    base_cmd = base_cmd_tmp[:arg_index] + base_cmd_tmp[arg_index + 3:]
    new_cmd_list = [" ".join(base_cmd) + " -j -g %2.4i%2.2i%2.2i %2.4i%2.2i%2.2i --no_filter --parts_out_flag" % (
        d1.year, d1.month, d1.day, d2.year, d2.month, d2.day)
                    for d1, d2 in zip(month_list, month_list_end)]
    print('new_cmd_list', new_cmd_list)
    print('end_cmd_tmp', end_cmd_tmp)
    # raise
    quit_flag = 0
    for cmd in new_cmd_list:
        os.system(cmd)

    west, east = [min(lon_range), max(lon_range)]
    south, north = [min(lat_range), max(lat_range)]

    dir_west = ot.regnum2strnum(round(west if not lon_range[0] else lon_range[0], 1))
    dir_east = ot.regnum2strnum(round(east if not lon_range[1] else lon_range[1], 1))
    dir_south = ot.regnum2strnum(round(south if not lat_range[0] else lat_range[0], 1))
    dir_north = ot.regnum2strnum(round(north if not lat_range[1] else lat_range[1], 1))

    glob_tot = glob.glob(currentdir + "plots/nc_dir/" + "*{}_{}_{}_{}_DNF_{}_{}*helper.txt".format(
        dir_west,
        dir_east,
        dir_south,
        dir_north, Day_Night_Flag.lower(), filetail.split('dthres')[0])
                         )
    print('glob_tot', glob_tot)
    print(currentdir + "plots/nc_dir/" + "*{}_{}_{}_{}_DNF_{}_{}*helper.txt".format(
        dir_west,
        dir_east,
        dir_south,
        dir_north, Day_Night_Flag.lower(),
        filetail.split('dthres')[0]))

    while len(glob_tot) != len(month_list) and quit_flag != 1:
        time.sleep(10)
        glob_tot = glob.glob(currentdir + "plots/nc_dir/" + "*{}_{}_{}_{}_DNF_{}_{}*helper.txt".format(
            dir_west,
            dir_east,
            dir_south,
            dir_north, Day_Night_Flag.lower(),
            filetail.split('dthres')[0])
                             )
        # why o why gaat dit fout voor de "all" run?
        out = os.popen('jobst | grep "%s" |grep "plot_avg_q"'%os.environ.get('USER')).read()
        print('len(out)', len(out), out)
        print('len(glob_tot)', len(glob_tot),len(month_list))
        if len(out) == 1:  # 1 as the run itself is also run with plot_avg_q when submitted
            quit_flag = 1

    # create a filelist for the next loop
    to_file_list = open(currentdir +
                        "plots/nc_dir/" +
                        "{}_{}_{}_{}_{}_{}_{}_final.txt".format('_'.join(dates),
                                                                dir_west,
                                                                dir_east,
                                                                dir_south,
                                                                dir_north,
                                                                Day_Night_Flag.lower(),
                                                                filetail), 'w')

    file_list = glob.glob(currentdir + "plots/nc_dir/" + "*{}_{}_{}_{}_DNF_{}_{}*helper.txt".format(
        dir_west,
        dir_east,
        dir_south,
        dir_north, Day_Night_Flag.lower(), filetail.split('dthres')[0])
                          )
    print('file_list')
    print(file_list)
    file_list.sort()
    [to_file_list.write('%s\n' % (uu[:-10] + 'level3.nc')) for uu in file_list]
    to_file_list.close()
    # remove temp txt files
    os.system('rm ' + currentdir + "plots/nc_dir/" + "*{}_{}_{}_{}_DNF_{}_{}*helper.txt".format(
        dir_west,
        dir_east,
        dir_south,
        dir_north, Day_Night_Flag.lower(), filetail.split('dthres')[0])
              )

    # restart script and read all monthly files?
    new_cmd = " ".join(end_cmd_tmp) + " --from_nc"
    os.system(new_cmd)
    return


def get_vrange(datakey):
    '''
    Function that return max and min vals based off data being plotted
    Returns None, None if not recognized
    '''
    if datakey in ['xretv','xa']:
        return [0, len(cbar_colors)]
    elif datakey in ['tot_col']:
        return [0, 3e16]
    elif datakey in ['DOF']:
        return [0, 1.3]
    else:
        return [None, None]


def get_datatitle(datakey):
    '''
    Function that return Title of the data based off data being plotted
    Will simply return string of input if not recognized
    '''
    if datakey in ["xretv"]:
        return "Surface NH3"
    elif datakey in ["tot_col"]:
        return "Total Column NH3"
    else:
        return str(datakey)


def unpack_and_plot(
        plotting_data,
        directory,
        datakey,
        Day_Night_Flag="All",
        lat_range=[None],
        lon_range=[None],
        save_flag=True,
        show_flag=False,
        kml_flag=False,
        filetail=None,
        title=None,
        maskoceans_flag=True,
        use_kml_file=True,
        date="",
        average_flag=False
):
    '''
    Function that takes a plotting data object (see retv_plot_tools, as
    well as the directory, and plots.
    Capable of generating saved files and kml files
    '''

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    import numpy as np
    import copy

    # We just want to plot filtered data so we will filter it here
    #print("plotting_data datakey", plotting_data[datakey])
    #print("min/max", [np.min(plotting_data[datakey]), np.max(plotting_data[datakey])])
    #print("plotting_data Filter Mask", plotting_data['Filter_Mask'])

    # This puts in fill values
    #plotting_data[datakey][(plotting_data['Filter_Mask']==0)]=-999.5

    # This just resizes the list to only include the valid points
    if average_flag:
        plotting_data[datakey][(plotting_data['Filter_Mask']!=0)]

    #print("min/max", [np.min(plotting_data[datakey]), np.max(plotting_data[datakey])])

    #exit



    if maskoceans_flag and not kml_flag:
        prt.printMessage("Removing data points over oceans", '--')
        plotting_data.remove_oceans(key=datakey)

    if kml_flag:
        try:
            import simplekml as sk
        except ImportError:
            prt.printError("Error loading simplekml module, must not be installed", exit_flag=False)
            kml_flag = False
            pass

    # Glitch in numpy makes annoying FutureWarning
    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)

    prt.printMessage("Unpacking of Values complete, now plotting Polygons for {}".format(datakey), "--")

    if not plotting_data[datakey].size:
        prt.printError("No plotting data for this day", exit_flag=False)
        return

    # Determine max/min values from datakey
    vmin, vmax = get_vrange(datakey)
    # For point density plots it will return None, vmin, vmax must be generated according to count values
    if vmin is None or vmax is None:
        vmin, vmax = [np.min(plotting_data[datakey]), np.max(plotting_data[datakey])]

    # Determine units and set up norm from datakey
    if datakey in ['xretv', 'xretv_bin','xa']:
        # norm = mpl.colors.Normalize(vmin=0.0, vmax=15.0)
        label = "ppbv"
    elif datakey in ['tot_col', 'tot_col_bin']:
        label = r'molecules cm$^{-2}$'
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        label = ""
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    # Setup KML, figure, and subplot
    fig, ax = plt.subplots(1, figsize=(10, 10))
    kml = sk.Kml() if kml_flag else None

    # List of lats & lons
    longitudes = plotting_data['Longitude']
    latitudes = plotting_data['Latitude']

    # Get bounds from the polygons
    if not lat_range[0] and not lon_range[0]:
        west, east = [
            min(longitudes) if min(longitudes) >= -180 else -180,
            max(longitudes) if max(longitudes) <= 180 else 180
        ]
        south, north = [
            min(latitudes) if min(latitudes) >= -90 else -90,
            max(latitudes) if max(latitudes) <= 90 else 90
        ]
    else:
        west, east = [min(lon_range), max(lon_range)]
        south, north = [min(lat_range), max(lat_range)]

    # overlay borders/coastlines to plot
    prt.printMessage("Initialising Basemap...", "--")

    rpt.drawBasemap(s=south, n=north, e=east, w=west, res='c')

    plt.scatter(west - 1, south - 1, c=0, vmin=vmin, vmax=vmax, s=10, cmap=cm)

    # Now plot the polygon on top of that
    # GENERATE POLYGONS AND COLORS, PLOT
    if 'polyColl' not in plotting_data.keys():
        plotting_data.polyGen()

    polyColl = copy.copy(plotting_data['polyColl'])
    white = (1., 1., 1.)
    if datakey == 'xretv' or datakey == 'xa':
        polycolors = [
            nonlinear_custom_color_xretv(
                plotting_data[datakey][i] * 1e3
            ) if plotting_data[datakey][i] >= 0 else white for i in xrange(plotting_data[datakey].shape[0])
        ]
    # if datakey == 'xretv':
    #     polycolors = [
    #         cm(
    #             norm(
    #                 plotting_data[datakey][i]*1e3
    #             )
    #         ) if plotting_data[datakey][i] >= 0 else white for i in xrange(plotting_data[datakey].shape[0])
    #     ]
    else:
        polycolors = [
            cm(
                norm(
                    plotting_data[datakey][i]
                )
            ) if plotting_data[datakey][i] >= 0 else white for i in xrange(plotting_data[datakey].shape[0])
        ]
    polyColl.set_color(polycolors)
    ax.add_collection(polyColl)

    # If the user wants a KML file,
    # an entirely different plotting routine is used,
    # to avoid doubling up on the regular plot
    if kml:

        # use of xrange(len()) for indices. Ugly, but faster than enumerate
        for ind in xrange(plotting_data[datakey].shape[0]):

            # Determine colour based off datakey
            if datakey == 'xretv' or datakey == 'xa':
                color = nonlinear_custom_color_xretv(1e3*(plotting_data[datakey][ind]))
            # if datakey == 'xretv':
            #     color = cm(norm(plotting_data[datakey][ind]))
            else:
                color = cm(norm(plotting_data[datakey][ind]))

            # Polygon for KML
            # If color is none, then point was masked
            # and should be avoided anyways
            if color is not None:
                pol = kml.newpolygon(altitudemode=sk.AltitudeMode.relativetoground)
                pol.outerboundaryis = zip(
                    plotting_data['polyLon'][ind],
                    plotting_data['polyLat'][ind],
                    [3000 for i in plotting_data['polyLat'][ind]]
                )
                pol.style.polystyle.color = rpt.convert_to_hex(color)
                pol.style.polystyle.fill = 1
                pol.style.polystyle.outline = 0

    # Format plot --------------------------------------------------------
    title = title if title is not None else "Date:  {}  Data:  {}".format(
        date,
        datakey,
    )
    plt.title(title, fontsize=22)
    plt.xlabel("\nLongitude", fontsize=18)
    plt.ylabel("Latitude\n\n", fontsize=18)
    plt.xlim(west, east)
    plt.ylim(south, north)

    # Formatting for colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    if datakey == "xretv" or datakey == 'xa':
        cbar.set_ticks(range(len(cbar_colors)))  # this range needs to be changed for Py3
        cbar.ax.set_yticklabels([str(i) for i in lower_lims] + ["^"])
        cbar.set_label(label, fontsize=16)
    else:
        cbar.set_ticks(np.linspace(vmin, vmax, len(cbar_colors) + 1))
        cbar.set_label(label, fontsize=16)

    # Generate filename
    savename = "{}_{}_{}_{}_{}_DNF_{}_{}_level3".format(
        date,
        ot.regnum2strnum(round(west if not lon_range[0] else lon_range[0], 1)),
        ot.regnum2strnum(round(east if not lon_range[0] else lon_range[1], 1)),
        ot.regnum2strnum(round(south if not lat_range[0] else lat_range[0], 1)),
        ot.regnum2strnum(round(north if not lat_range[0] else lat_range[1], 1)),
        Day_Night_Flag.lower(),
        filetail if filetail else datakey
    )

    # For saving a PNG image
    if save_flag:
        # Generate directory string
        plotdir = "{}{}".format(directory, "plots/")
        # Make directory if it doesn't exist
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        # Save final KML
        if kml:
            kml.save("{}{}.kml".format(plotdir, savename))
        # Save the plot
        plt.savefig("{}{}.png".format(plotdir, savename))
        print("Saving file to: {}{}.png".format(plotdir, savename))
        # for rico pdf
        # plt.savefig("{}{}.pdf".format(plotdir, savename))

    # If show flag, show
    if show_flag:
        plt.show()

    # clear plot for memory
    plt.close("all")
    return


# ----------------------------------------------------------------------------
# FUNCTIONS FOR BINNED AND AVERAGED PLOTS
# ----------------------------------------------------------------------------


def plot_average_of_range(dates,
                          directories,
                          lat_range,
                          lon_range,
                          plots,
                          show_flag=False,
                          save_flag=True,
                          use_combined_files=False,
                          kml_flag=False,
                          binsize=binsize,
                          oversample_binsize=oversample_binsize,
                          dict_flag=True,
                          wgt_thresh_perc=.05,
                          day_thresh_perc=.05,
                          no_filter=False,
                          quality_flags=[4],
                          Day_Night_Flag="All",
                          nthreads=1,
                          to_nc=False,
                          from_nc=False,
                          remove_outliers=False
                          ):
    '''
    Plots average over many directories, given binning and oversample size
    '''
    from multiprocessing.pool import Pool
    import retv_plot_tools as rpt
    import gc
    from functools import partial
    for datak in plots:
        prt.printMessage("Plotting {}".format(datak), "==")

    prt.printMessage(
        "Importing {} directories in {} threads...".format(
            len(directories),
            nthreads
        ),
        "--"
    )
    print("to_nc: ", to_nc)
    print("from_nc: ", from_nc)
    if not from_nc:
        if nthreads == 1:
            comb_plt_dat = []
            comb_append = comb_plt_dat.append
            if use_combined_files:
                [comb_append(rpt.unpack_plotdata_monthly_yearly(directory, datakey=plots,
                                                                use_combined_files=use_combined_files)) for
                 directory in directories]
            else:
                [comb_append(rpt.unpack_plotdata(directory, datakey=plots)) for directory in
                 directories]
        else:
            io_pool = Pool(processes=nthreads)
            print(use_combined_files, 'multi_thread')
            if use_combined_files:
                comb_plt_dat = io_pool.map(unpack_helper, [(direc, plots) for direc in directories])
            else:
                comb_plt_dat = io_pool.map(unpack_helper_single, [(direc, plots) for direc in directories])
            print("done pool, closing")
            io_pool.close()
            print("closed pool, joining")
            io_pool.join()
            print("pool joined")

        prt.printMessage("Combining data into one object...", "--")
        global top_plotting_data
        top_plotting_data = rpt.retvplotdat(comb_plt_dat)
        if not top_plotting_data.keys():
            prt.printError("No results for this range.", exit_flag=False)
            # continue
        prt.printMessage("Filtering Data...", "--")
        print(top_plotting_data.keys())
        top_plotting_data.filterData(
            lat_range=lat_range,
            lon_range=lon_range,
            quality_flags=quality_flags,
            Day_Night_Flag=Day_Night_Flag
        )
        if any([top_plotting_data[key].size == 0 for key in top_plotting_data.keys()]):
            print("No data returned from filter, skipping directory.")
            # continue
        # Sort points:
        top_plotting_data.sortAscending()

        # Done with these, delete it for space
        del comb_plt_dat, directories, quality_flags
        gc.collect()

        if any(
                [
                    top_plotting_data[key].size == 0 for key in top_plotting_data.keys()
                ]
        ):
            print("\nData returned empty after filtering, continuing to next datakey\n")
            # continue

        prt.printMessage("Binning and averaging data", "=-")
    if not from_nc:
        xretv_data = bin_and_avg(
            top_plotting_data,
            plots,
            lat_range=lat_range,
            lon_range=lon_range,
            binsize=binsize,
            oversample_binsize=oversample_binsize,
            wgt_thresh_perc=wgt_thresh_perc,
            day_thresh_perc=day_thresh_perc,
            no_filter=no_filter,
            Day_Night_Flag=Day_Night_Flag,
            dict_flag=dict_flag,
            nthreads=nthreads,
            remove_outliers=remove_outliers
        )
        # Send to plotting routine. Regular datakey, then point density
        #print('test',plots, type(plots),len(plots))
        if type(plots)!=list or len(plots)==1:
            title = make_title_format(
                get_datatitle(plots[0] if len(plots)==1 else plots),
                dates,
                Day_Night_Flag)
        else:
            title = make_title_format(
                get_datatitle("all"),
                dates,
                Day_Night_Flag)
        #title = "{}  |  {} to {:}\nDay/Night: {}".format(
        #    get_datatitle(datakey),
        #    min(top_plotting_data['Date']),
        #    max(top_plotting_data['Date']),
        #    Day_Night_Flag
        #)

        date = "{}_{}".format(
            min(dates), max(dates)
        )
    if from_nc:
        xretv_data = rpt.retvplotdat(
            read_from_nc(currentdir, dates, plots, lat_range=lat_range, lon_range=lon_range,
                         filetail="{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
                             plots[0] if len(plots)==1 else 'all',
                             binsize,
                             oversample_binsize,
                             wgt_thresh_perc,
                             day_thresh_perc), Day_Night_Flag=Day_Night_Flag, wgt_thresh_perc=wgt_thresh_perc,
                         day_thresh_perc=day_thresh_perc))
        # Send to plotting routine. Regular datakey, then point density
        if type(plots)!=list or len(plots)==1:
            title = make_title_format(
                get_datatitle(plots[0] if len(plots)==1 else plots),
                dates,
                Day_Night_Flag)
        else:
            title = make_title_format(
                get_datatitle("all"),
                dates,
                Day_Night_Flag)
        #title = "{}  |  {} to {:}\nDay/Night: {}".format(
        #    get_datatitle(datakey),
        #    min(dates),
        #    max(dates),
        #    Day_Night_Flag
        #)
        date = "{}_{}".format(
            min(dates), max(dates)
        )
    if not xretv_data:
        prt.printError("No Plotting Data")

    if to_nc:
        if type(plots)!=list or len(plots)==1:
            filetail = "{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
                                 plots[0] if len(plots) == 1 else plots,
                                 binsize,
                                 oversample_binsize,
                                 wgt_thresh_perc,
                                 day_thresh_perc)
        else:
            filetail = "{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
                                 'all',
                                 binsize,
                                 oversample_binsize,
                                 wgt_thresh_perc,
                                 day_thresh_perc)
        write_to_nc_file(currentdir, date, xretv_data, plots, lat_range=lat_range, lon_range=lon_range,
                         filetail=filetail, Day_Night_Flag=Day_Night_Flag)
        # now quit and return
        # exit()  # remove when all is working as planned
        # return
    # xretv_data.printStats()

    # Generate polygons once:
    xretv_data.polyGen()
    print('keys', xretv_data.keys())
    for key in xretv_data.keys():
        print(key, np.shape(xretv_data[key]))
    # Regular binned average plot
    for datak in plots:
        title = make_title_format(
                get_datatitle(datak),
                dates,
                Day_Night_Flag)
        unpack_and_plot(
            xretv_data,
            currentdir,
            datak,
            Day_Night_Flag=Day_Night_Flag,
            lat_range=lat_range,
            lon_range=lon_range,
            save_flag=save_flag,
            show_flag=show_flag,
            date=date,
            title=title,
            use_kml_file=True,
            kml_flag=kml_flag,
            filetail="{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
                datak,
                binsize,
                oversample_binsize,
                wgt_thresh_perc,
                day_thresh_perc
            ),
            average_flag=average_flag
        )

    # Total points in bin + oversampled region
    if type(plots)!=list or len(plots)==1:
        title = title.replace(get_datatitle(plots[0] if len(plots) == 1 else plots), "Total Point Density")
    else:
        title = title.replace(get_datatitle("all"), "Total Point Density")
    unpack_and_plot(
        xretv_data,
        currentdir,
        "N",
        Day_Night_Flag=Day_Night_Flag,
        lat_range=lat_range,
        lon_range=lon_range,
        save_flag=save_flag,
        show_flag=show_flag,
        date=date,
        title=title,
        kml_flag=kml_flag,
        filetail="{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
            "TotalPointDensity",
            binsize,
            oversample_binsize,
            wgt_thresh_perc,
            day_thresh_perc
        ),
            average_flag=average_flag
    )
    if type(plots)!=list or len(plots)==1:
        title = title.replace(get_datatitle(plots[0] if len(plots) == 1 else plots), "Bin Point Density")
    else:
        title =  title.replace(get_datatitle("all"), "Bin Point Density")
    # Binned total only
    unpack_and_plot(
        xretv_data,
        currentdir,
        "N_Bin_Points",
        Day_Night_Flag=Day_Night_Flag,
        lat_range=lat_range,
        lon_range=lon_range,
        save_flag=save_flag,
        show_flag=show_flag,
        date=date,
        title=title,
        kml_flag=kml_flag,
        filetail="{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
            "BinPointDensity",
            binsize,
            oversample_binsize,
            wgt_thresh_perc,
            day_thresh_perc
        ),
            average_flag=average_flag
    )
    if type(plots)!=list or len(plots)==1:
        title = title.replace(get_datatitle(plots[0] if len(plots) == 1 else plots), "Date Density")
    else:
        title = title.replace(get_datatitle("all"), "Date Density")
    # Daily density
    unpack_and_plot(
        xretv_data,
        currentdir,
        "N_Dates",
        Day_Night_Flag=Day_Night_Flag,
        lat_range=lat_range,
        lon_range=lon_range,
        date=date,
        title=title,
        kml_flag=kml_flag,
        filetail="{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
            "DateDensity",
            binsize,
            oversample_binsize,
            wgt_thresh_perc,
            day_thresh_perc
        ),
            average_flag=average_flag
    )

    print("")

    prt.printMessage("Plot complete", "-=")

    return


def bin_and_avg(
        comb_plot_dat,
        datakey,
        lat_range,
        lon_range,
        binsize=binsize,
        oversample_binsize=oversample_binsize,
        wgt_thresh_perc=.05,
        day_thresh_perc=.05,
        no_filter=False,
        Day_Night_Flag="All",
        dict_flag=True,
        nthreads=1.,
        remove_outliers=False
):
    '''
    Function that bins and averages the points when given in the form
    of a dictionary with Lats, Lons, and a datakey
    '''

    # ------------------------------------------------------------------------
    # SETUP
    # ------------------------------------------------------------------------
    import numpy as np
    from multiprocessing.pool import Pool
    import retv_plot_tools as rpt
    import gc

    # Set up grid
    # minlat, maxlat = [min(comb_plot_dat['Latitude']), max(comb_plot_dat['Latitude'])]
    # minlon, maxlon = [min(comb_plot_dat['Longitude']), max(comb_plot_dat['Longitude'])]
    minlat, maxlat = [min(lat_range), max(lat_range)]
    minlon, maxlon = [min(lon_range), max(lon_range)]

    lonbins, latbins = [np.arange(minlon, maxlon + binsize, binsize), np.arange(minlat, maxlat + binsize, binsize)]
    lon_arr, lat_arr = np.meshgrid(lonbins, latbins)

    # Determine number of days, and thus the oversampling threshold
    total_dates = [int(tmpdate) for tmpdate in list(set(comb_plot_dat['Date']))]
    num_days = len(total_dates)

    day_threshold = round((day_thresh_perc * num_days), 0)

    # Initialize list to store the polygons and average ammonia vals
    xretv_data, pointcount_data = [[] for i in xrange(2)]
    # Output dictionary (to be written to file)
    outdict = {}
    if type(datakey)==list:
        out_keys = [
            "Latitude",
            "Longitude",
            "N_Bin_Points",
            "N_Oversampled_Points",
            "N",
            "N_Dates",
            "Lat_BL",
            "Lat_TL",
            "Lat_TR",
            "Lat_BR",
            "Lon_BL",
            "Lon_TL",
            "Lon_TR",
            "Lon_BR",
        ]
        for tmp_key in datakey:
            out_keys += [
                tmp_key + "_sdev",
                tmp_key + "_median",
                tmp_key]
    else:
        out_keys = [
                       "Latitude",
                       "Longitude",
                       "sdev",
                       "median",
                       "N_Bin_Points",
                       "N_Oversampled_Points",
                       "N",
                       "N_Dates",
                       "Lat_BL",
                       "Lat_TL",
                       "Lat_TR",
                       "Lat_BR",
                       "Lon_BL",
                       "Lon_TL",
                       "Lon_TR",
                       "Lon_BR",datakey]

    #        out_keys_2 = [
    #                str(len(lonbins)),
    #                str(len(latbins)),
    #        ]

    # Determine denominator for weight threshold:
    print('get_total_weights')
    max_weight = rpt.get_total_weight(binsize, oversample_binsize)

    # ----------------------------------------------------------------------------
    # BEGIN BINNING & AVERAGING
    # ----------------------------------------------------------------------------

    # loop over all possible coordinates
    gc.collect()  # free up some space
    grid_vector = zip(lon_arr.ravel(), lat_arr.ravel())
    print('grid_vector', len(grid_vector))
    print(datakey)
    if type(datakey) is list:
        print("Data key was not a string and its type was changed")
    #print(datakey)
    pool_list = [
        (
            coord if type(coord) is tuple else tuple(coord),
            binsize,
            oversample_binsize,
            '#'.join(datakey),
            wgt_thresh_perc,
            max_weight,
            day_threshold, no_filter
        ) for coord in grid_vector
    ]
    print(grid_vector[0])
    print(pool_list[0])

    print(grid_vector[9])
    print(pool_list[9])

    if nthreads == 1:
        bin_data = []
        bin_append = bin_data.append
        print('pool items')
        [bin_append(avg_by_bin(pool_item)) for pool_item in pool_list]
        fin_data = [i[0] for i in bin_data]
        out_data = [i[1] for i in bin_data]
    else:
        print("Opening pool of {} threads".format(nthreads))
        try:
            p = Pool(processes=nthreads)
            # p = Pool(processes=int(10))  # for now max on 40, see what happens
        except OSError:
            print("ERROR: Encountered error in Initialising pool.")
            print("Retrying with half as many threads")
            nthreads = np.floor(
                nthreads / 2.
            ) if np.floor(nthreads / 2.) >= 1 else 1
            try:
                p = Pool(processes=int(nthreads))
            except OSError:
                print("ERROR: Pooling failed again")
                print("Reducing threads to 2")
                p = Pool(processes=2)

        if not remove_outliers:
            bin_data = p.map(avg_by_bin, pool_list)
        else:
            bin_data = p.map(avg_bin_helper, [(pool_item, remove_outliers) for pool_item in pool_list])
        print("Pool done, now closing")
        p.close()
        print("Pool closed, now joining")
        p.join()
        print("Pool joined")
        fin_data = [i[0] for i in bin_data]
        out_data = [i[1] for i in bin_data]
        for n in range(20):
            print(bin_data[n][1])
    print('reached_final2')
    final_data = rpt.retvplotdat(fin_data)
    if dict_flag:
        prt.printMessage("Writing output gridded file...", "--")

        outdict_filename = "{}_{}_{}_{}_{}_{}_DNF_{}_{}_bin_{}_os_{}_wgt_{}_dthres_{}_level3".format(
            min(comb_plot_dat['Date']),
            max(comb_plot_dat['Date']),
            ot.regnum2strnum(round(minlon, 1)),
            ot.regnum2strnum(round(maxlon, 1)),
            ot.regnum2strnum(round(minlat, 1)),
            ot.regnum2strnum(round(maxlat, 1)),
            Day_Night_Flag,
            datakey,
            binsize,
            oversample_binsize,
            wgt_thresh_perc,
            day_thresh_perc
        )
        print('test1',out_keys)
        # for n in range(100):
        #     print(out_data[n])
        out_data = rpt.combine_outdata(out_data)

        for outkey in out_keys:
            if outkey in ["Lat_BL", "Lon_BL","Lat_TL", "Lon_TL","Lat_TR", "Lon_TR","Lat_BR", "Lon_BR"]:
                continue
            outdict[outkey] = out_data[outkey][:]
            print('out',outkey, outdict.keys())
        for ind, [latkey, lonkey] in enumerate([
            ("Lat_BL", "Lon_BL"),
            ("Lat_TL", "Lon_TL"),
            ("Lat_TR", "Lon_TR"),
            ("Lat_BR", "Lon_BR")
        ]):
            outdict[latkey] = out_data['polyLat'][:, ind]
            outdict[lonkey] = out_data['polyLon'][:, ind]
        print('out',outdict.keys())
        write_dict(outdict, outdict_filename, header_order=out_keys)

    del bin_data, fin_data

    return final_data


def avg_by_bin(gridinfo, remove_outliers=False):
    import numpy as np
    '''
    Takes in grid data, and computes the average value within the bin/grid

    STRUCTURE:
    ((BLLON, BLLAT), bin, os_bin, datakey)
    '''
    global top_plotting_data

    BL_lon, BL_lat = gridinfo[0]
    binsize = gridinfo[1]
    os_binsize = gridinfo[2]
    datakey = gridinfo[3].split('#')
    #print(datakey)
    wgt_thresh_perc = gridinfo[4]
    max_weight = gridinfo[5]
    day_threshold = gridinfo[6]
    no_filter = gridinfo[7]

    bindict = {}

    tmp_minlon = BL_lon
    tmp_maxlon = tmp_minlon + binsize
    tmp_minlat = BL_lat
    tmp_maxlat = tmp_minlat + binsize

    tmp_lon = tmp_minlon + binsize / 2.
    tmp_lat = tmp_minlat + binsize / 2.

    os_minlon = tmp_minlon - os_binsize
    os_maxlon = tmp_maxlon + os_binsize
    os_minlat = tmp_minlat - os_binsize
    os_maxlat = tmp_maxlat + os_binsize

    os_indices = (
            (top_plotting_data['Longitude'] >= os_minlon) &
            (top_plotting_data['Longitude'] <= os_maxlon) &
            (top_plotting_data['Latitude'] >= os_minlat) &
            (top_plotting_data['Latitude'] <= os_maxlat)
    )

    num_points = np.sum(os_indices)
    if type(datakey) == list:
        for key in ['Latitude', 'Longitude', 'Date'] + list(datakey):
            bindict[key] = top_plotting_data[key][os_indices]
    else:
        for key in ['Latitude', 'Longitude', 'Date',datakey]:
            bindict[key] = top_plotting_data[key][os_indices]

    # print(os_indices)
    # print(bindict[datakey])
    for key_tmp in datakey:
        if not bindict[key_tmp].size:
            return None, None

    bin_indices = (
            (bindict['Longitude'] >= tmp_minlon) &
            (bindict['Longitude'] <= tmp_maxlon) &
            (bindict['Latitude'] >= tmp_minlat) &
            (bindict['Latitude'] <= tmp_maxlat)
    )

    num_binned = np.sum(bin_indices)
    num_oversampled = num_points - num_binned
    # Find distances for each point
    bindict['distance'] = np.sqrt((bindict['Longitude'] - tmp_lon) ** 2 + (bindict['Latitude'] - tmp_lat) ** 2)
    # Find weights corresponding to distances
    bindict['weight'] = stat.get_weights(bindict['distance'], wgt_flag='gaus')
    # if type(datakey)!=list:
    #     statdict = stat.stat_mv(bindict, datakey, remove_outliers=remove_outliers)
    #     bindict[datakey] = statdict['mean_wgt']
    #     bindict['sdev'] = statdict['sdev']
    #     bindict['median'] = statdict['median']
    # else:
    for key_tmp in datakey:
        statdict = stat.stat_mv(bindict, key_tmp, remove_outliers=remove_outliers)
        bindict[key_tmp] = statdict['mean_wgt']
        bindict[key_tmp+'_sdev'] = statdict['sdev']
        bindict[key_tmp+'_median'] = statdict['median']
    bindict['polyLon'] = np.array([
        tmp_minlon,
        tmp_minlon,
        tmp_minlon + binsize,
        tmp_minlon + binsize
    ])
    bindict['polyLat'] = np.array([
        tmp_minlat,
        tmp_minlat + binsize,
        tmp_minlat + binsize,
        tmp_minlat
    ])
    # Lat lon coords of bin
    bindict['Longitude'] = tmp_lon
    bindict['Latitude'] = tmp_lat
    # Density Data
    bindict['N_Bin_Points'] = num_binned
    bindict['N_Oversampled_Points'] = num_oversampled
    bindict['N'] = num_points
    bindict['N_Dates'] = len(list(set(bindict['Date'])))
    # add in wgt threshold / day threshold here
    # add in wgt threshold / day threshold here
    # add in wgt threshold / day threshold here
    # add in wgt threshold / day threshold here
    # add in wgt threshold / day threshold here
    # add in wgt threshold / day threshold here
    bindict['Filter_Mask'] = 1
    # remove other one if split flag included AND add a statement to skip this
    if ((bindict['N_Dates'] < day_threshold) or
        (np.sum(bindict['weight']) / max_weight < wgt_thresh_perc)) and not no_filter:
        for remove_datakey in datakey:
            #No need to include fill value if filter mask is used
            #bindict[remove_datakey] = -999.5
            bindict['Filter_Mask'] = 0
    # before writing out give total weight
    bindict['Weight_Frac'] = (np.sum(bindict['weight']) / max_weight)
    bindict['weight'] = np.sum(bindict['weight'])
    bindict['Max_weight'] = max_weight
    plotdict = {}
    # if type(datakey)!=list:
    #     plot_keys = [
    #         "Longitude",
    #         "Latitude",
    #         "polyLon",
    #         "polyLat",
    #         "N_Bin_Points",
    #         "N_Oversampled_Points",
    #         "N",
    #         "N_Dates",
    #         "Weight_Frac",
    #         "weight",
    #         "Max_weight",
    #         "sdev",
    #         "median",
    #         "Filter_Mask",
    #         datakey
    #     ]
    # else:
    plot_keys = [
        "Longitude",
        "Latitude",
        "polyLon",
        "polyLat",
        "N_Bin_Points",
        "N_Oversampled_Points",
        "N",
        "N_Dates",
        "Weight_Frac",
        "weight",
        "Max_weight",
        "Filter_Mask"]
    for tmp_key in datakey:
        plot_keys += [
        tmp_key+"_sdev",
        tmp_key+"_median",
        tmp_key]
    for plotkey in plot_keys:
        plotdict[plotkey] = np.array([bindict[plotkey]])
    return plotdict, bindict


def write_dict(outdict, filename, header_order=None):
    '''
    Takes a binned average dictionary and writes it out to a file
    For Validation efforts
    '''
    import pandas as pd
    import os.path

    #for outkey in outdict.keys():
        # print(outdict[outkey])

    plotdirdict = "{}{}".format(currentdir, "plots/")
    # Make directory if it doesn't exist
    if not os.path.exists(plotdirdict):
        os.makedirs(plotdirdict)
    # Path+filename
    completefilename = os.path.join(plotdirdict, filename)

    df = pd.DataFrame(outdict)
    print('ehm1',header_order)
    if header_order:
        df[header_order].to_csv(
            "{}{}".format(
                plotdirdict,
                "{}.csv".format(filename)
            ),
            index=False,
            header=True
        )
        txt = df[header_order].to_string(
            index=False,
            header=True
        )
    else:
        df.to_csv(
            "{}{}".format(
                plotdirdict,
                "{}.csv".format(filename)
            ),
            index=False,
            header=True
        )

        txt = df.to_string(
            index=False,
            header=True
        )
    with open("{}.txt".format(completefilename), 'w') as f:
        f.write(txt)


def queue_average_dates(cmd_list, nthreads, dates):
    '''
    Generate average plot as a jobfile and submit to the queue
    '''

    # Remove extra args from sys.argv[:]
    base_cmd = [
        arg for arg in cmd_list if arg not in [
            '-j', '--jobsub', "--silent", "--nosave", "--show"
        ]
    ]

    if nthreads > 44:
        msg = "Too many requested threads. Capping at 44"
        raise RuntimeWarning(msg)
        nthreads = 44

    # Average plots are really variable in binsize and oversample size, etc
    # therefore I am setting it as the max time until a strong way of
    # Determing the duration is determined.
    # est_time = 21600
    # est_time = 3600
    est_time = 600 * nthreads
    if est_time > 10800:
        est_time = 10800
    #if est_time > 21600:
    #    est_time = 21600
    mem_limit = 2000 * nthreads
    if mem_limit > 200000:
        mem_limit = 200000
    # if mem_limit > 100000:
    #     mem_limit = 100000

    # Join remaining arguments
    cmd_list = [" ".join(base_cmd)]

    # Write Job file
    jobfile = prt.write_jobfile(currentdir, cmd_list, "plot_avg_queue_%s" % ('_'.join(dates)), est_time, ncpus=nthreads,
                                mem=mem_limit)

    # Submit job to queue
    prt.submit_job(jobfile)


def queue_average(cmd_list, nthreads):
    '''
    Generate average plot as a jobfile and submit to the queue
    '''

    # Remove extra args from sys.argv[:]
    base_cmd = [
        arg for arg in cmd_list if arg not in [
            '-j', '--jobsub', "--silent", "--nosave", "--show"
        ]
    ]

    if nthreads > 44:
        msg = "Too many requested threads. Capping at 44"
        raise RuntimeWarning(msg)
        nthreads = 44

    # Average plots are really variable in binsize and oversample size, etc
    # therefore I am setting it as the max time until a strong way of
    # Determing the duration is determined.
    # est_time = 21600
    # est_time = 3600
    est_time = 1200 * nthreads
    if est_time > 10800:
        est_time = 10800
    # if est_time > 21600:
    #     est_time = 21600
    mem_limit = 2000 * nthreads
    if mem_limit > 200000:
        mem_limit = 200000
    # if mem_limit > 100000:
    #     mem_limit = 100000
    # Join remaining arguments
    cmd_list = [" ".join(base_cmd)]

    # Write Job file
    jobfile = prt.write_jobfile(currentdir, cmd_list, "plot_avg_queue_beta2", est_time, ncpus=nthreads, mem=mem_limit)

    # Submit job to queue
    prt.submit_job(jobfile)


# ----------------------------------------------------------------------------
# FUNCTIONS FOR GLOBAL PIXEL PLOTS
# ----------------------------------------------------------------------------


def plot_globe_for_one_day(
        directories,
        lat_range,
        lon_range,
        plots,
        show_flag=False,
        save_flag=True,
        kml_flag=False,
        quality_flags=[4],
        Day_Night_Flag="All",
        nthreads=1,
        average_flag=False
):
    '''
    Accepts a directory of cdf files and plots it. in our case, this is 1 day.
    '''
    import retv_plot_tools as rpt
    from multiprocessing import Pool

    for datakey in plots:

        prt.printMessage(
            "Importing {} directories in {} threads...".format(
                len(directories),
                nthreads
            ),
            "--"
        )

        if nthreads == 1:
            comb_plt_dat = []
            comb_append = comb_plt_dat.append
            [
                comb_append(
                    rpt.unpack_plotdata(directory)
                ) for directory in directories
            ]
        else:
            io_pool = Pool(processes=nthreads)

            comb_plt_dat = io_pool.map(rpt.unpack_plotdata, directories)
            print("done pool, closing")
            io_pool.close()
            print("closed pool, joining")
            io_pool.join()
            print("pool joined")

        plotting_data = rpt.retvplotdat(comb_plt_dat)
        del comb_plt_dat

        # Filter data based off flags
        prt.printMessage("Filtering Data...", "--")
        plotting_data.filterData(
            lat_range=lat_range,
            lon_range=lon_range,
            quality_flags=quality_flags,
            Day_Night_Flag=Day_Night_Flag
        )
        if any([plotting_data[key].size == 0 for key in plotting_data.keys()]):
            print("No data returned from filter, skipping directory.")
            continue
        plotting_data.polyGen()

        # Generate Title
        date = plotting_data['Date'][0]
        #title = "{}  |  {}\nDay/Night: {}".format(
        #    get_datatitle(datakey),
        #    date,
        #    Day_Night_Flag
        #)
        title = make_title_format(
            get_datatitle(datakey),
            plotting_data['Date'],
            Day_Night_Flag)

        tmplonr, tmplatr = ct.latlon_fromdir(directory)
        # set footprint locations
        plotting_data['polyLat'] = [rpt.compute_ellipse_footprint(lat, lon, range, azi, zenith, cris_fovDia)[0] for
                                    lat, lon, range, azi, zenith in
                                    zip(plotting_data['Latitude'], plotting_data['Longitude'],
                                        plotting_data['SatelliteRange'],
                                        plotting_data['SatelliteAzimuthAngle'], plotting_data['SatelliteZenithAngle'])]
        plotting_data['polyLon'] = [rpt.compute_ellipse_footprint(lat, lon, range, azi, zenith, cris_fovDia)[1] for
                                    lat, lon, range, azi, zenith in
                                    zip(plotting_data['Latitude'], plotting_data['Longitude'],
                                        plotting_data['SatelliteRange'],
                                        plotting_data['SatelliteAzimuthAngle'], plotting_data['SatelliteZenithAngle'])]

        # Generate title
        # title = make_title_format(
        #     get_datatitle(datakey),
        #     plotting_data['Data'],
        #     Day_Night_Flag)
        #title = "{}  |  {} to {}\nDay/Night: {}".format(
        #    get_datatitle(datakey),
        #    min(plotting_data['Date']),
        #    max(plotting_data['Date']),
        #    Day_Night_Flag
        #)
        # Unpack new list of data and plot.
        unpack_and_plot(
            plotting_data,
            currentdir,
            datakey,
            Day_Night_Flag=Day_Night_Flag,
            lat_range=lat_range,
            lon_range=lon_range,
            save_flag=save_flag,
            show_flag=show_flag,
            date="{}_{}".format(min(plotting_data['Date']), max(plotting_data['Date'])),
            title=title,
            kml_flag=kml_flag,
            average_flag=average_flag
        )

        print("")
        prt.printMessage("Plot complete", "-=")

    return


def queue_daily_world(
        cmd_list,
        dates,
        lat_range,
        lon_range,
        name="plt_sglr_sgld_queue",
        ncpus=40
):
    '''
    Generate custom region plot commands into individual days to divide among
    CPUs and submit to queue
    '''
    import datetime
    import math

    start = datetime.datetime.strptime(dates[0], "%Y%m%d").date()
    end = datetime.datetime.strptime(dates[1], "%Y%m%d").date() if len(dates) == 2 else (
            start + datetime.timedelta(days=0))
    start_str = start.strftime("%Y%m%d")
    end_str = end.strftime("%Y%m%d")

    base_cmd = [
        arg for arg in cmd_list if arg not in [
            '-j', '--jobsub', "-g", start_str, end_str, "--silent", "--nosave", "--show"
        ]
    ]

    num_regions = math.ceil(max(lat_range) - min(lat_range)) * math.ceil(max(lon_range) - min(lon_range))
    time_per_day_per_region = 150
    time_per_day = num_regions * time_per_day_per_region
    max_time = 1800

    prt.write_and_submit_daily(
        base_cmd,
        dates,
        lat_range,
        lon_range,
        time_per_day,
        max_time=max_time,
        name=name,
        ncpus=ncpus
    )


# ----------------------------------------------------------------------------
# FUNCTIONS FOR SIGNLE REGION DAILY PLOTS
# ----------------------------------------------------------------------------


def plot_sgl_rgn_for_one_day(
        directory,
        lat_range,
        lon_range,
        plots,
        show_flag=False,
        save_flag=True,
        kml_flag=False,
        quality_flags=[4],
        Day_Night_Flag="All",
        average_flag=False
):
    '''
    Accepts a directory of cdf files and plots it. in our case, this is 1 day.
    '''
    import retv_plot_tools as rpt

    for datakey in plots:

        prt.printMessage("Beginning Directory:\n{}".format(directory), "=-")

        # Append date to list for future labels/titles
        date = "{:%Y%m%d}".format(ct.date_fromdir(directory))

        # get plotting data, filter
        plotting_data = rpt.unpack_plotdata(directory)  # , use_pickle=use_pickle)
        plotting_data = rpt.retvplotdat(plotting_data)

        prt.printMessage("Filtering Data...", "--")
        plotting_data.filterData(
            lat_range=lat_range,
            lon_range=lon_range,
            quality_flags=quality_flags,
            Day_Night_Flag=Day_Night_Flag
        )
        plotting_data.sortAscending(key='tot_col')

        if any([plotting_data[key].size == 0 for key in plotting_data.keys()]):
            print("No data returned from filter, skipping directory.")
            continue

        # Generate Title
        date = plotting_data['Date'][0]
        #title = "{}  |  {}\nDay/Night: {}".format(
        #    get_datatitle(datakey),
        #    date,
        #    Day_Night_Flag
        #)
        title = make_title_format(
            get_datatitle(datakey),
            plotting_data['Date'],
            Day_Night_Flag)

        tmplonr, tmplatr = ct.latlon_fromdir(directory)
        # set footprint locations
        plotting_data['polyLat'] = [rpt.compute_ellipse_footprint(lat, lon, range, azi, zenith, cris_fovDia)[0] for
                                    lat, lon, range, azi, zenith in
                                    zip(plotting_data['Latitude'], plotting_data['Longitude'],
                                        plotting_data['SatelliteRange'],
                                        plotting_data['SatelliteAzimuthAngle'], plotting_data['SatelliteZenithAngle'])]
        plotting_data['polyLon'] = [rpt.compute_ellipse_footprint(lat, lon, range, azi, zenith, cris_fovDia)[1] for
                                    lat, lon, range, azi, zenith in
                                    zip(plotting_data['Latitude'], plotting_data['Longitude'],
                                        plotting_data['SatelliteRange'],
                                        plotting_data['SatelliteAzimuthAngle'], plotting_data['SatelliteZenithAngle'])]
        # print(tmp_foot[1])
        # print(np.shape(tmp_foot))
        # raise
        # Unpack new list of data and plot.
        unpack_and_plot(
            plotting_data,
            prt.formatDir(
                os.path.dirname(directory[:-1])
            ),
            datakey,
            Day_Night_Flag=Day_Night_Flag,
            lat_range=tmplatr,
            lon_range=tmplonr,
            save_flag=save_flag,
            show_flag=show_flag,
            date=date,
            title=title,
            kml_flag=kml_flag,
            average_flag=average_flag
        )

        print("")
        prt.printMessage("Plot complete", "-=")

    return


def make_title_format(datakey, dates, Day_Night_Flag):
    import calendar
    """
    Helper function that looks at the difference between the starting and
    ending dates and outputs a title in format that is easiest to read.

    Example:
        Date given: 2018/01/01 2018/01/31
        Format returned: Surface NH3 | Month of Janurary 2018

        Date given: 2018/01/01 2018/01/01
        Format returned: Surface NH3 | 1 Janurary 2018

        Date given: 2018/01/01 2018/12/31
        Format returned: Surface NH3 | Year of 2018

    Args:
        dates(list): list of datetime objects
        datakey(str): type of data obtained

    Returns:
        str: to be used for the title
    """
    print("Creating plot title")
    start_date = min(dates)
    end_date = max(dates)
    if type(start_date) is unicode:
        start_date = start_date.encode('ascii','ignore')
        end_date = end_date.encode('ascii','ignore')
    if type(start_date) is str:
        start_date = datetime.date(int(start_date[:4]),
                                   int(start_date[4:6]),
                                   int(start_date[6:8]))
        end_date = datetime.date(int(end_date[:4]),
                                 int(end_date[4:6]),
                                 int(end_date[6:8]))
    print(start_date)
    last_day_month = calendar.monthrange(end_date.year,
                                         end_date.month)[1]
    if start_date == end_date:
        title = "{data} | {day} {month} {year}".format(
            data=datakey,
            day=end_date.strftime("%d"),
            month=end_date.strftime("%b"),
            year=end_date.strftime("%Y"))

    elif (start_date.month == 1 and start_date.day == 1 and
          end_date.month == 12 and end_date.day == 31):

        if start_date.year == end_date.year:
            title = "{data} | Year of {year}".format(
                data=datakey,
                year=start_date.year)
        else:
            title = "{data} | {s_year} to {e_year}".format(
                data=datakey,
                s_year=start_date.year,
                e_year=end_date.year)

    elif (start_date.day == 1 and end_date.day == last_day_month):

        if start_date.month == end_date.month:
            title = "{data} | Month of {month}, {year}".format(
                data=datakey,
                month=end_date.strftime("%B"),
                year=end_date.strftime("%Y"))

        else:
            title = "{data} | {year} {s_mon} to {e_mon}".format(
                data=datakey,
                year=start_date.year,
                s_mon=start_date.strftime("%b"),
                e_mon=end_date.strftime("%b"))

    else:
        title = "{data} | {s_d}/{s_m}/{s_y} to {e_d}/{e_m}/{e_y}".format(
            data=datakey,
            s_d=start_date.day,
            s_m=start_date.month,
            s_y=start_date.month,
            e_d=end_date.day,
            e_m=end_date.month,
            e_y=end_date.month)

    return (title)#+ "\n Day/Night: {}".format(Day_Night_Flag))


def queue_daily_sgl(
        cmd_list,
        dates,
        lat_range,
        lon_range,
        name="plt_sglr_sgld_queue",
        ncpus=40
):
    '''
    Divide single region/day commands into individual days to divide among CPUs
    and submits to queue
    '''
    import datetime

    start = datetime.datetime.strptime(dates[0], "%Y%m%d").date()
    end = datetime.datetime.strptime(dates[1], "%Y%m%d").date() if len(dates) == 2 else (
            start + datetime.timedelta(days=0))
    start_str = start.strftime("%Y%m%d")
    end_str = end.strftime("%Y%m%d")

    base_cmd = [
        arg for arg in cmd_list if arg not in [
            '-j', '--jobsub', "-g", start_str, end_str, "--silent", "--nosave", "--show"
        ]
    ]

    time_per_day_per_region = 120
    max_time = 1800

    prt.write_and_submit_daily(
        base_cmd,
        dates,
        lat_range,
        lon_range,
        time_per_day_per_region,
        max_time=max_time,
        name=name,
        ncpus=ncpus
    )


# ----------------------------------------------------------------------------
# MAIN FUNCTIONS
# ----------------------------------------------------------------------------
def run_main(cmd_list,
             dates,
             directories,
             lat_range,
             lon_range,
             plots,
             use_combined_files=False,
             save_flag=True,
             show_flag=False,
             kml_flag=False,
             binsize=binsize,
             oversample_binsize=oversample_binsize,
             dict_flag=True,
             wgt_thresh_perc=0.05,
             day_thresh_perc=0.05,
             no_filter=False,
             quality_flags=[4],
             Day_Night_Flag="All",
             world_flag=world_flag,
             average_flag=average_flag,
             nthreads=1,
             plot_by_parts_flag=False,
             to_nc=False,
             from_nc=False,
             remove_outliers=False,
             parts_out_flag=False
             ):
    # maybe in future we need more plots flags?
    print(plots)
    if len(plots)>1:
        plots_fig = "all"
    else:
        plots_fig = plots[0]
    filetail = "{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
        plots_fig,
        binsize,
        oversample_binsize,
        wgt_thresh_perc,
        day_thresh_perc)
    if plot_by_parts_flag:
        launch_parts_jobs(cmd_list, dates, currentdir, lon_range, lat_range, filetail=filetail,
                          Day_Night_Flag=Day_Night_Flag)
        exit()
    prt.printMessage("Plotting Retrieval", "==")
    prt.printMessage("Inputted Arguments", "=-")
    print(cmd_list)
    prt.printMessage("Directories:", " -")
    for directory in directories:
        print(directory)
    prt.printMessage("Latitude Range:", " -")
    print("{} to {}".format(lat_range[0], lat_range[1]))
    prt.printMessage("Longitude Range:", " -")
    print("{} to {}".format(lon_range[0], lon_range[1]))
    # dates

    # If the user has flagged a global plots, averaged over range
    if average_flag:
        plot_average_of_range(dates,
                              directories,
                              lat_range,
                              lon_range,
                              plots,
                              use_combined_files=use_combined_flag,
                              save_flag=save_flag,
                              show_flag=show_flag,
                              kml_flag=kml_flag,
                              binsize=binsize,
                              oversample_binsize=oversample_binsize,
                              dict_flag=dict_flag,
                              wgt_thresh_perc=wgt_thresh_perc,
                              day_thresh_perc=day_thresh_perc,
                              no_filter=no_filter,
                              quality_flags=quality_flags,
                              Day_Night_Flag=Day_Night_Flag,
                              nthreads=nthreads,
                              to_nc=to_nc,
                              from_nc=from_nc,
                              remove_outliers=remove_outliers
                              )
    # Else if the user has asked for single days global plots, no averaging/combining
    elif world_flag:
        plot_globe_for_one_day(
            directories,
            lat_range,
            lon_range,
            plots,
            save_flag=save_flag,
            show_flag=show_flag,
            kml_flag=kml_flag,
            quality_flags=quality_flags,
            Day_Night_Flag=Day_Night_Flag,
            nthreads=nthreads,
            average_flag=average_flag
        )
    # Regular, single region plot
    elif not world_flag and not average_flag:
        for directory in directories:
            plot_sgl_rgn_for_one_day(
                directory,
                lat_range,
                lon_range,
                plots,
                save_flag=save_flag,
                show_flag=show_flag,
                kml_flag=kml_flag,
                quality_flags=quality_flags,
                Day_Night_Flag=Day_Night_Flag,
                average_flag=average_flag
            )
    # check if helper files for multi-part calculations are needed
    if parts_out_flag:
        date = "{}_{}".format(
            min(dates), max(dates)
        )
        filetail = "{}_bin_{}_os_{}_wgt_{}_dthres_{}".format(
            plots[0] if len(plots)==1 else 'all',
            binsize,
            oversample_binsize,
            wgt_thresh_perc,
            day_thresh_perc)
        savename = "{}_{}_{}_{}_{}_DNF_{}_{}_helper.txt".format(
            date,
            ot.regnum2strnum(round(west if not lon_range[0] else lon_range[0], 1)),
            ot.regnum2strnum(round(east if not lon_range[0] else lon_range[1], 1)),
            ot.regnum2strnum(round(south if not lat_range[0] else lat_range[0], 1)),
            ot.regnum2strnum(round(north if not lat_range[0] else lat_range[1], 1)),
            Day_Night_Flag.lower(),
            filetail)
        aa = open(currentdir + "plots/nc_dir/"+savename,'w')
        aa.write('Filler')
        aa.close()
    return


if __name__ == "__main__":

    # Set up arguments
    parser = argparse.ArgumentParser(
        description='Tool to plot a finished retrieval')
    parser.add_argument(
        '-d',
        '--dir',
        type=str,
        nargs='*',
        required=False,
        default=[],
        help='Retrieval directory'
    )
    parser.add_argument(
        '-g',
        type=str,
        nargs='*',
        required=False,
        default=[],
        metavar="YYYYMMDD",
        help='Use inputted date range and run over the global folders.'
    )
    parser.add_argument(
        '--gdir',
        type=str,
        nargs='*',
        required=False,
        default=[],
        help="Global top directory for '-g' searching"
    )
    parser.add_argument(
        '--lat_range',
        type=float,
        nargs=2,
        required=False,
        default=[],
        help='Range of latitudes you wish to observe (for --globaldir only)'
    )
    parser.add_argument(
        '--lon_range',
        type=float,
        nargs=2,
        required=False,
        default=[],
        help='Range of longitudes you wish to observe (for --globaldir only)'
    )
    parser.add_argument(
        '-p',
        '--plot_flag',
        type=str,
        nargs="*",
        required=False,
        default="xretv",
        help='What plots youd like to see. (xretv or tot_col)'
    )
    parser.add_argument(
        '-q',
        '--quality_flags',
        type=str,
        nargs="*",
        required=False,
        default=[4],
        help='What quality flags you would like to include'
    )
    parser.add_argument(
        '-n',
        '--daynight',
        type=int,
        required=False,
        default=2,
        help='0 for night plot, 1 for day plot, 2 for daylight, 3 for no daylight, 4 for all data'
    )
    parser.add_argument(
        '-w',
        '--world',
        action="store_true",
        required=False,
        default=False,
        help='Whether or not you would like the multiple regions plotted on one graph'
    )
    parser.add_argument(
        '-a',
        '--average',
        action="store_true",
        required=False,
        default=False,
        help='Whether or not you would like to perform averaging across these days'
    )
    parser.add_argument(
        '-b',
        '--bin',
        type=float,
        required=False,
        default=binsize,
        help='Main bin size for averaging'
    )
    parser.add_argument(
        '-o',
        '--oversample_binsize',
        type=float,
        required=False,
        default=oversample_binsize,
        help='Oversample bin size for averaging.(How many degrees past edge of bin)'
    )
    parser.add_argument(
        '--wgt',
        type=float,
        nargs=1,
        required=False,
        default=[.05],
        help="Percentage weight threshold (between 0, 1) DEFAULT: 0.05"
    )
    parser.add_argument(
        '--daythres',
        type=float,
        nargs=1,
        required=False,
        default=[.05],
        help="Percentage total days threshold (between 0, 1) DEFAULT: 0.05"
    )
    # parser.add_argument(
    #     '--new_pickle',
    #     action="store_true",
    #     required=False,
    #     default=False,
    #     help='Flag to control whether to generate a new pickle file'
    # )
    parser.add_argument(
        '--use_combined',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control whether to use combined files instead of pickles for reading'
    )
    parser.add_argument(
        '--nosave',
        action="store_false",
        required=False,
        default=True,
        help='Flag to control whether to not save plots'
    )
    parser.add_argument(
        '--kml',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control whether to not make a KML file'
    )
    parser.add_argument(
        '--show',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control whether to show plots'
    )
    parser.add_argument(
        '--dict',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control whether to generate gridded file'
    )
    parser.add_argument(
        '-j',
        '--jobsub',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control whether to use the queue'
    )
    parser.add_argument(
        '-c',
        '--cpus',
        type=int,
        required=False,
        default=44,
        help='Flag to control whether how many CPU it will use on the queue. does not aply to "-a" plots'
    )
    parser.add_argument(
        '-t',
        '--threads',
        type=int,
        required=False,
        default=1,
        help='How many threads/processes to split the binning and averaging into'
    )
    parser.add_argument(
        '--silent',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control whether any output is shown'
    )
    parser.add_argument(
        '--parts',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control a split into monthly parts, and should usually be used for multi-year plots,\n'+
             'this flag is only needed for continent sized regions\n'+
             'and should be used when making multi-year plots.\n'
    )
    parser.add_argument(
        '--parts_out_flag',
        action="store_true",
        required=False,
        default=False
    )
    parser.add_argument(
        '--to_nc',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control if level3 nc files are to be created\n'+
             '(essential to switch on when using --parts)'
    )
    parser.add_argument(
        '--from_nc',
        action="store_true",
        required=False,
        default=False,
        help='Flag to control if plots are read from the level 3 nc files\n'+
             'Only used when rerunning a failed plot'
    )
    parser.add_argument(
        '--no_filter',
        action="store_true",
        required=False,
        default=False,
        help='Flag to override plotting filters\n'+
             'Usually only used by the overhead script that controls the --parts run\n'+
            'When running with --parts you dont have to switch this on'
    )
    parser.add_argument(
        '--outlier',
        action="store_true",
        required=False,
        default=False,
        help='Flag to switch on the removal of outliers'
    )
    # Unpack Args
    args = parser.parse_args()
    dates = args.g[:]
    gdir = args.gdir[0] if args.gdir else ct.get_globaldir()
    plots = args.plot_flag[:] if (
            "list" in str(type(args.plot_flag[:])).lower()
    ) else [args.plot_flag[:]]
    if plots ==['all']:
        plots = [u'xa', u'dt_nh3', u'tsfc', u'xretv', u'tsfc_ap', u'tot_col_meas_error', u'tot_col_total_error', u'tstd_ap', u'tot_col', u'tsfc_for', u'DOF', u'nedt', u'snr']
    print('plots',plots)
    lat_range = [
        min(args.lat_range[:]), max(args.lat_range[:])
    ] if args.lat_range[:] else [None, None]

    lon_range = [
        min(args.lon_range[:]), max(args.lon_range[:])
    ] if args.lon_range[:] else [None, None]

    weight_threshold_percent = args.wgt[0]
    day_threshold_percent = args.daythres[0] if not args.no_filter else 0.00

    # plot flags
    world_flag = args.world
    average_flag = args.average
    binsize = args.bin
    oversample_binsize = args.oversample_binsize
    # use_pickle = not args.new_pickle
    use_combined_flag = args.use_combined
    # print('usepickle', use_pickle)
    # raise ValueError
    quality_flags = args.quality_flags[:]

    if args.daynight == 0:
        Day_Night_Flag = "Night"
    elif args.daynight == 1:
        Day_Night_Flag = "Day"
    elif args.daynight == 2:
        Day_Night_Flag = "DayLight"
    elif args.daynight == 3:
        Day_Night_Flag = "No_DayLight"
    else:
        Day_Night_Flag = "All"
    print('Day_Night_Flag', Day_Night_Flag)
    # Some plots take too long to run interactively, must be run on queue
    queue_flag = args.jobsub

    # Flags for showing, saving
    show_flag = args.show
    save_flag = args.nosave
    kml_flag = args.kml
    dict_flag = args.dict
    plot_by_parts_flag = args.parts
    parts_out_flag = args.parts_out_flag
    from_nc = args.from_nc
    to_nc = args.to_nc
    no_filter = args.no_filter
    remove_outliers = args.outlier
    # If the silent option is flagged, redirect any output to /dev/null
    if args.silent:
        f = open(os.devnull, 'w')
        sys.stdout = f

    # Generate dirs based of user input, and set up whether or not this plot will be a ranged average
    if dates and not args.dir:
        # check for length of series, decide for use of daily/monthly/yearly files
        start = datetime.datetime.strptime(dates[0], "%Y%m%d").date()
        end = datetime.datetime.strptime(dates[1], "%Y%m%d").date() if len(dates) == 2 else (
                start + datetime.timedelta(days=0))

        ndays = (end - start).days
        if ndays >= 28 and use_combined_flag:  # and ndays < 365:
            plot_date_type = 'monthly'
            # Find all directories for these dates in global dir
            directories = [prt.formatDir(dire) + "cdf/" for dire in
                           ct.generate_monthly_dirs(dates, lat_range, lon_range, gdir=gdir)]
        # elif ndays >= 365:
        #     plot_date_type = 'yearly'
        #     # Find all directories for these dates in global dir
        #     directories = [prt.formatDir(dire) + "cdf/" for dire in
        #                    ct.generate_yearly_dirs(dates, lat_range, lon_range, gdir=gdir)]
        else:
            plot_date_type = 'daily'
            # Find all directories for these dates in global dir
            directories = [prt.formatDir(dire) + "cdf/" for dire in
                           ct.generate_dirs(dates, lat_range, lon_range, gdir=gdir)]

    elif dates and args.dir:
        # check for length of series, decide for use of daily/monthly/yearly files
        start = datetime.datetime.strptime(dates[0], "%Y%m%d").date()
        end = datetime.datetime.strptime(dates[1], "%Y%m%d").date() if len(dates) == 2 else (
                start + datetime.timedelta(days=0))

        ndays = (end - start).days
        if ndays >= 28 and use_combined_flag:  # and ndays < 365:
            plot_date_type = 'monthly'
            # Find all directories for these dates in global dir
            directories = [prt.formatDir(dire) + "cdf/" for dire in
                           ct.generate_monthly_dirs_choice_dir(dates, lat_range, lon_range, dire=args.dir, gdir=gdir)]
        # elif ndays >= 365:
        #     plot_date_type = 'yearly'
        #     # Find all directories for these dates in global dir
        #     directories = [prt.formatDir(dire) + "cdf/" for dire in
        #                    ct.generate_yearly_dirs(dates, lat_range, lon_range, gdir=gdir)]
        else:
            plot_date_type = 'daily'
            # Find all directories for these dates in global dir
            directories = [prt.formatDir(dire) + "cdf/" for dire in
                           ct.generate_dirs_choice_dir(dates, lat_range, lon_range, dire=args.dir, gdir=gdir)]


    elif ((not dates) and (args.dir)):
        # If there are no inputted dates, but there are inputted directories
        directories = []
        for directory in list(args.dir[:]):
            directories.extend([dire + "cdf/" for dire in ct.find_dirs(directory)])

    else:
        # If neither, take current direectory and find all directories beneath
        prt.printMessage("Finding directories recursively (this may take a while)", "--")
        directories = [dire + "cdf/" for dire in ct.find_dirs(currentdir)]
    # quick killer
    # raise ValueError
    if queue_flag and not world_flag and not average_flag:
        ncpus = args.cpus
        queue_daily_sgl(sys.argv[:], dates, lat_range, lon_range, name="plt_sglr_sgld_queue", ncpus=ncpus)
    elif queue_flag and world_flag and not average_flag:
        ncpus = args.cpus
        queue_daily_world(sys.argv[:], dates, lat_range, lon_range, name="plt_sglr_sgld_queue", ncpus=ncpus)
    elif queue_flag and average_flag:
        queue_average_dates(sys.argv[:], args.threads, dates)
    else:
        run_main(sys.argv[:],
                 dates,
                 directories,
                 lat_range,
                 lon_range,
                 plots,
                 use_combined_files=use_combined_flag,
                 save_flag=save_flag,
                 show_flag=show_flag,
                 kml_flag=kml_flag,
                 binsize=binsize,
                 oversample_binsize=oversample_binsize,
                 dict_flag=dict_flag,
                 wgt_thresh_perc=weight_threshold_percent,
                 day_thresh_perc=day_threshold_percent,
                 no_filter=no_filter,
                 quality_flags=quality_flags,
                 Day_Night_Flag=Day_Night_Flag,
                 world_flag=world_flag,
                 average_flag=average_flag,
                 nthreads=args.threads,
                 plot_by_parts_flag=plot_by_parts_flag,
                 to_nc=to_nc,
                 from_nc=from_nc,
                 remove_outliers=remove_outliers,
                 parts_out_flag=parts_out_flag
                 )