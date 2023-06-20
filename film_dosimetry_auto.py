#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# file_name.py
# Python 3.8.3
"""
Author: Joseph Bateman
Created: Fri Nov  5 16:01:25 2021
Modified: Fri Nov  5 16:01:25 2021

Description 
-------------

"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from PIL import Image
from astropy.modeling import models, fitting
from scipy.signal import savgol_filter

# class holding fitting and dose measurement methods
class film_analyser:
    # initialise by defining calibration constants
    # and retrieving dose in seperate channels from given file
    def __init__(self, filename):
        # seperate film into rgb channels and return data values
        def get_rgb_d_from_file(filename):
            # open file and convert to Image object in RGB format for analysis
            input_image = Image.open(filename).convert("RGB")
            # convert to numpy array for ease
            input_imarray = np.array(input_image)
            # X_d reweighted because reasons
            X_d = np.array(input_imarray / 65535)
            # transpose carried out to retrieve channels from matrix
            X_d_transpose = np.transpose(X_d)
            # values in r,g,b colour channels retrieved
            rX_d, gX_d, bX_d = X_d_transpose
            return rX_d, gX_d, bX_d

        # method to convert image matrix to dose with calibration constants
        def get_dose(a, b, c, X_d):
            # conversion formula and matrix manipulation to retrieve dose
            dose = np.fliplr(np.rot90((a - c * X_d) / (X_d - b), 3))
            return dose
        
        #ebt3 calibration params up to 25 Gy

        #a_r = 9.08646370e-03
        #b_r = 4.83772917e-04
        #c_r = 3.35869684e+00
        
        #a_g = 1.25584443e-02
        #b_g = 2.13539882e-04
        #c_g = 4.98839944e+00
        
        #a_b = 1.84179296e-02
        #b_b = 3.86486899e-04
        #c_b = 1.10265417e+01

        #ebt3 full calibration params

       # a_r = 1.08054260e-02
       # b_r = 3.96464552e-04
       # c_r = 4.03747782e+00

       # a_g = 1.53583642e-02
       # b_g = 9.04451292e-05
       # c_g = 6.24521113e+00

       # a_b = 2.25123199e-02
       # b_b = 2.75485087e-04
       # c_b = 1.36316506e+01
       
       #mdv3 cnew alibration params

        a_r = 7.21071436e-02
        b_r = 4.67275120e-04
        c_r = 2.38511751e+01

        a_g = 1.11306451e-01
        b_g = 3.18670259e-04
        c_g = 3.79674212e+01

        a_b = 2.32990669e-01
        b_b = 7.33015188e-04
        c_b = 9.54995622e+01
        
        #mdv3 old calibration params

        #a_r = 7.33992768e-02
        #b_r = 4.73098274e-04
        #c_r = 2.47012235e+01

        #a_g = 1.16746035e-01
        #b_g = 2.36675318e-04
        #c_g = 4.08090283e+01

        #a_b = 2.77850507e-01
        #b_b = 4.96471261e-04
        #c_b = 1.17878125e+02

        
        # retrieve red, green and blue image data from file
        rX_d, gX_d, bX_d = get_rgb_d_from_file(filename)
        # retrieve dose in each colour channel
        self.dose_red = get_dose(a_r, b_r, c_r, rX_d)
        self.dose_green = get_dose(a_g, b_g, c_g, gX_d)
        self.dose_blue = get_dose(a_b, b_b, c_b, bX_d)

    # vile function to automatically crop images
    # janky as shit, whoever's reading this pls make it better

    def remove_borders(
        self, dose_profile, left_xlim=80, right_xlim=350, up_ylim=90, down_ylim=350
    ):
        # collapse 2d matrix into 1d matrix of median values
        def collapse_along_axis(matrix):
            # initialise list for holding result
            collapsed_matrix = []
            # loop through matrix and add median of each row to 1D array
            # median better than mean here as it smooths out extrema
            for i in matrix:
                collapsed_matrix.append(np.median(i))
            # return 1D array of medians
            return collapsed_matrix

        # collapse matrix to get median x dose distributions along y axis
        collapsed_dose_along_y = collapse_along_axis(dose_profile)
        # collapse matrix to get median y dose distributions along x axis
        collapsed_dose_along_x = collapse_along_axis(np.transpose(dose_profile))
        # smooth data even further so that only significant spikes remain
        cdy_hat = savgol_filter(collapsed_dose_along_y, 41, 1)
        cdx_hat = savgol_filter(collapsed_dose_along_x, 41, 1)

        # uncomment lines for relevant graphs showing the smoothed fit

        # y_pixels = np.arange(0, len(dose_profile), 1)
        # x_pixels = np.arange(0, len(dose_profile_T), 1)

        # plt.figure()
        # plt.plot(y_pixels, cdy_hat)
        # plt.plot(y_pixels, collapsed_dose_along_y)
        # plt.yscale("log")
        # plt.figure()

        # plt.plot(x_pixels, cdx_hat)
        # plt.plot(x_pixels, collapsed_dose_along_x)

        # define maximum limits for cropping
        # decrease l variables and increase r variables to reduce max cropping
        lx_sp = left_xlim
        rx_sp = right_xlim
        ly_sp = up_ylim
        ry_sp = down_ylim

        # defines pixel limits in x and y
        # represent limits of what respective variables above can be set tp
        x_ltp = 0
        x_rtp = len(collapsed_dose_along_x)
        y_ltp = 0
        y_rtp = len(collapsed_dose_along_y)

        # search for increase in dose moving outwards from centre
        # and label detected edge coordinates
        # starting at max cropping limit

        # -ve x edge

        for i in range(0, lx_sp):
            if cdx_hat[lx_sp - i] < cdx_hat[lx_sp - i - 1]:
                x_ltp = lx_sp - i
                break
        # +ve x edge
        for i in range(rx_sp, len(cdx_hat)-1):
            if cdx_hat[i] < cdx_hat[i + 1]:
                x_rtp = i
                break
        # -ve y edge
        for i in range(0, ly_sp):
            if cdy_hat[ly_sp - i] < cdy_hat[ly_sp - i - 1]:
                y_ltp = ly_sp - i
                break
        # +ve y edge
        for i in range(ry_sp, len(cdy_hat)-1):
            if cdy_hat[i] < cdy_hat[i + 1]:
                y_rtp = i
                break
        # truncate 2d dose profile in y at the detected edges
        dose_profile = dose_profile[y_ltp : y_rtp + 1]
        # truncate 2d dose profile in x at the detected edges
        # save to new variable because numpy doesn't like in-place modification
        new_dose_profile = []
        for i in range(len(dose_profile)):
            new_dose_profile.append(dose_profile[i][x_ltp : x_rtp + 1])
        # convert new list variable to array
        new_dose_profile = np.array(new_dose_profile)
        # return cropped profile
        return new_dose_profile
        # print(len(dose_profile))
        # plt.plot(x_pixels, cdx_hat)
        # plt.plot(x_pixels, collapsed_dose_along_x)
        # plt.yscale("log")
        # plt.figure(figsize=(8, 2.5))
        # plt.subplot(1, 3, 1)
        # plt.imshow(
        #    new_dose_profile, cmap="jet", interpolation="nearest", vmin=0, vmax=10
        # )
        # plt.title("Data")
        # plt.colorbar()

    # fit 2d gaussian to get beam size and return fit
    def fit_2D_gaussian_to_dose(self, dose_profile, filename):
        # get x and y coordinates from dose data
        y, x = np.mgrid[: np.shape(dose_profile)[0], : np.shape(dose_profile)[1]]
        # define Astropy Gaussian model
        g = models.Gaussian2D(x_mean=175, y_mean=195)
        # intialise least squares fitter
        fitter = fitting.LevMarLSQFitter()
        # carry out fit to data
        p = fitter(g, x, y, dose_profile)
        # plot dose, fit and residuals
        plt.figure(figsize=(8, 2.5))

        plt.subplot(1, 3, 1)
        plt.imshow(dose_profile, cmap="jet", interpolation="nearest", vmin=0, vmax=50)
        plt.title(filename + " Data")
        plt.colorbar()
        plt.subplot(1, 3, 2)
        plt.imshow(p(x, y), cmap="jet", interpolation="nearest", vmin=0, vmax=50)
        plt.title("Gaussian Fit")
        plt.colorbar()
        plt.subplot(1, 3, 3)
        plt.imshow(
            dose_profile - p(x, y), cmap="jet", interpolation="nearest", vmin=-2, vmax=2
        )
        plt.title("Residual")
        plt.colorbar()
        # return fit
        return p

    # construct table of fitting parameters
    def make_fit_params_table(self, red_fit, green_fit, blue_fit):
        # create lists of data and labels for ease of looping through
        fit_list = [red_fit, green_fit, blue_fit]
        colours = ["Red", "Green", "Blue"]
        # initialise lists for holding fit data
        amplitudes = []
        x_means = []
        y_means = []
        x_stddevs = []
        y_stddevs = []
        # fill with data from astropy fitting object
        for i in fit_list:
            amplitudes.append(i.amplitude)
            x_means.append(i.x_mean)
            y_means.append(i.y_mean)
            # I presume there's a reason for the weighting factor here?
            x_stddevs.append(i.x_stddev * 25.5 / 300)
            y_stddevs.append(i.y_stddev * 25.4 / 300)
        # merge data into pandas DataFrame object
        fitting_data = pd.DataFrame(
            {
                "Colour": colours,
                "Amplitude": amplitudes,
                "x_mean": x_means,
                "y_mean": y_means,
                "x_stddev": x_stddevs,
                "y_stddev": y_stddevs,
            }
        )
        # return DataFrame object holding fitting data
        return fitting_data


# path must end with a forward slash
# function to loop through and analyse every film image in a directory
def analyse_all_films(path_to_directory_holding_images):
    # call film analyser object to retrieve dose and get gaussian fitting data
    def get_fitting_data_table_and_dose_green(filename):
        # initialise film analyser object
        film = film_analyser(filename)
        # get doses in colour channels
        dose_red, dose_green, dose_blue = film.dose_red, film.dose_green, film.dose_blue
        # crop images
        # add extra arguments to remove_borders() to define cropping limit
        dose_red = film.remove_borders(dose_red)
        dose_blue = film.remove_borders(dose_blue)
        dose_green = film.remove_borders(dose_green)
        # fit and retrieve fitting parameters for each colour
        for i in range(len(filename)):
            if filename[len(filename) - i - 1] == "/":
                filename = filename[len(filename) - i : len(filename) - 5]
                break

        fitted_red = film.fit_2D_gaussian_to_dose(dose_red, filename)
        fitted_green = film.fit_2D_gaussian_to_dose(dose_green, filename)
        fitted_blue = film.fit_2D_gaussian_to_dose(dose_blue, filename)

        # create table holding all fitting data for film
        fitting_data_table = film.make_fit_params_table(
            fitted_red, fitted_green, fitted_blue
        )
        # return table and green dose for further processing
        return fitting_data_table, dose_green

    # return max dose and save to file
    def get_max_dose(fitting_data_table, dose_green):

        # x mean position green channel
        xmean = fitting_data_table.iloc[1, 2] * 25.4 / 300

        # y mean position green channel
        ymean = fitting_data_table.iloc[1, 3] * 25.4 / 300

        dim = 10
        # find max dose
        max_dose = np.round(
            np.mean(
                dose_green[
                    int(fitting_data_table.iloc[1, 3] - dim / 2) : int(
                        fitting_data_table.iloc[1, 3] + dim / 2
                    ),
                    int(fitting_data_table.iloc[1, 2] - dim / 2) : int(
                        fitting_data_table.iloc[1, 2] + dim / 2
                    ),
                ]
            ),
            2,
        )
        # print max dose
        print(max_dose, fitting_data_table.iloc[1, 4], fitting_data_table.iloc[1, 5])
        # save results to .csv file
        with open("/Users/josephbateman/CLEAR/Film/CHUV_plasmids_25_10_22/plasmids1.csv", "a", newline="") as csv_file:
            print(
                max_dose,
                fitting_data_table.iloc[1, 4],
                fitting_data_table.iloc[1, 5],
                filename,
                file=csv_file,
                sep=",",
            )


    # get directory in bytes format
    directory = os.fsencode(path_to_directory_holding_images)
    # begin loop for analysing every file in folder
    for file in os.listdir(directory):

        # decode from bytes into string
        file_str = file.decode("utf-8")
        if file_str[0] != ".":
            filename = file_str
            # add path to directory to filename to point correctly
            file_str = path_to_directory_holding_images + file_str
            # generate and retrieve fitting table and green_dose
            fitting_data_table, dose_green = get_fitting_data_table_and_dose_green(
                file_str
            )
            print("File: " + filename)
            # print fitting tables
            print(fitting_data_table)
            # retrieve and print max_dose
            get_max_dose(fitting_data_table, dose_green)


# main function
def main():
    # carry out full analysis of films in chosen directory
    analyse_all_films("/Users/josephbateman/CLEAR/Film/CHUV_plasmids_25_10_22/run1/")


main()
