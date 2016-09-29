#!/usr/bin/env python
# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 fenc=utf-8
# Author: Yang Tong
# Created: 2016-07-24 10:15 SGT

#Application: This script applies cubic spline interpolation to fit lines.
#             The given x-y data should be put into a file in which there are just
#             two columns: the first one is x and the second y.
#             >>>Multiple lines can be fitted simultaneously if they have the same x grids.
#                And the x-y data for each line should be separated from each other by a empty line
#How to use: >>>python cubic_spline_interpolation.py
#History:
#       created at 2016-07-24 10:15 SGT
#       modified at 2016-9-29 15:08 SGT

from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import re

def read_x_y(filename):
    f = open(filename, 'r')

    x_points, y_points = [], []
    for line in f:
        x_y = re.findall("[-+0-9\.eE]+", line)
        if x_y:
            x, y = float(x_y[0]), float(x_y[1])
            if y >= 10000:
                continue
            if x not in x_points:
                x_points.append(x)
            ind = x_points.index(x)
	    if len(y_points) == ind:
	        y_points.append([])
	    if len(y_points) < ind:
	        raise Exception, "some mistake"
	    y_points[ind].append(y)
    f.close()

    assert len(x_points) == len(y_points), "The number of x data does not match that of y data. Please check!"
    length = len(y_points[0])
    for item in y_points:
        assert len(item) == length, "The provided data is not complete. Please Check!"
    y_ = []
    for ind in range(len(y_points[0])):
        y_.append([item[ind] for item in y_points])
    return (x_points, y_)

def multiple_lines_cubic_spline_interpolation(x, y_s, start, end, accuracy):
    y_ = []
    for y in y_s:
        x_new, y_new = cubic_spline_interpolation(x, y, start, end, accuracy)
        y_.append(y_new)
    return x_new, y_

def cubic_spline_interpolation(x, y, start, end, accuracy):
    tck = interpolate.splrep(x, y, s=0)

    x_new = np.arange(start, end, accuracy*0.1)
    y_new = interpolate.splev(x_new, tck, der=0)
    return x_new, y_new

def plot_and_save_fig(known_x_y, new_x_y, operation):
    '''Plot the polynomial or save the figure:
       code\toperation
       0\tquit
       1\tonly plot
       2\tonly save
       3\tplot & save'''

    if operation == '0':
        return 0

    if "y" in raw_input("Is the given data also included in the figure? (y or n): ").lower():
        for y in known_x_y[1]:
            plt.plot(known_x_y[0], y, 's') #, label="Given Data", markersize=6, linewidth=2)
    for y in new_x_y[1]:
        plt.plot(new_x_y[0], y, "r-") #, label="Cubic Spline", linewidth=2)
    plt.xlim([new_x_y[0][0]*1.1-0.1*new_x_y[0][-1], 1.10*new_x_y[0][-1]-0.1*new_x_y[0][0]])
    plt.legend(loc='best')
    plt.xlabel(raw_input("Enter the x title: "))
    plt.ylabel(raw_input("Enter the y title: "))
    plt.title(raw_input("Enter the title: "))
    if operation in ('2', '3'):
        plt.savefig("Cubic_Spline.png")
    if operation in ('1', '3'):
        plt.show()

def save_interpolated_data(x, y_s, filename):
    f = open(filename, "w")


    for y in y_s:
        for i, j in zip(x, y):
            f.write("%f\t%f\n" % (i, j))
        f.write("\n")
    f.close()


def cal_significant_figure(accuracy, figs=0):
    while accuracy < 1:
        figs += 1
        accuracy *= 10
    return figs

def input_data():
    print("\n>>>Cubic spline interpolation is employed to give an interpolating polynomial<<<\n")
    filename = raw_input("Please enter the filename where each line includes a known x-y pair: ")
    accuracy = float(raw_input("Please enter the accuracy (e.g.0.1, 0.01, 0.001, ...): accuracy="))
    sig_figs = cal_significant_figure(accuracy)
    return filename, accuracy, sig_figs



if __name__ == "__main__":

    filename, accuracy, sig_figs = input_data()
    x_points, y_points = read_x_y(filename)
    x, y = multiple_lines_cubic_spline_interpolation(x_points, y_points, x_points[0], x_points[-1], accuracy)
    if "y" in raw_input("\n>>>Print the min points for each interpolated line up to the accuary give above? (y or n): ").lower():
        output_format = ">>>For curve %d, E_min=%f at %." + str(sig_figs) + "f<<<"
        for ind, y_ in enumerate(y):
            print(output_format % (ind, y_.min(), x[list(y_).index(y_.min())]))

    if "y" in raw_input("\n>>>Save the cubic spline interpolated data? (y or n): ").lower():
        filename = raw_input("Enter the filename to store the interpolated data:")
        save_interpolated_data(x, y, filename)

    if "y" in raw_input("\n>>>Plot or save the figure? (y or n): ").lower():
        print(plot_and_save_fig.__doc__)
        code = raw_input("Enter the code: ").strip(" ")[0]
        plot_and_save_fig([x_points, y_points], [x, y], code)



