"""
Tools for calculating and plotting percentiles.
Often used for plotting resolutions.

Tom Stuttard
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

from builtins import range
from builtins import zip
from future import standard_library

standard_library.install_aliases()
from utils.plotting.ext_modules import *

from utils.maths.hist import get_intermediate_points


def y_percentiles(
    # Data
    x,  # Array of x data points (1D, N elements)
    y,  # Array of y data points (1D, N elements)
    # Steer calc -> y axis
    percentiles,  # List of percentile values, e.g. [50.,90.]
    # Steer calc -> x axis
    # Manual steps case
    xmin=None,  # First x value in the scan
    xmax=None,  # Last x value in the scan
    xwidth=None,  # x width of band to sample at each scan point (optionally can be a function f(x) returning a width as a function of x (useful for power law distributed data)
    nx=None,  # Number of x scan points
    # Binned case
    xbins=None,  # Alternatively just provide a binning
    # Plot steering
    ax=None,  # Provide a matplotlib `ax` if want to plot
    color=None,
    percentile_linestyles=None,
    percentile_labels=None,
    plot_points=False,  # Optionally can plot the percentile values at eahc scan point (can be used to verify the smoothing procedure)
    # Tuning
    min_num_events=1,  # Min num events in a scan step required to calculate the percentile
    use_mean_x_pos=False,  # If True, the x value of the scan point is the mean it the data points in that band (recommended), otherwise use the center of the window
    # Smoothing
    smooth=None,  # Smoothing to apply: None, "interpolate", gaussian_filter"
    smooth_kw=None,  # kwargs to pass to smoothing function
):
    """
    Compute percentile bands for the y dimension of a 2D (x,y) set of data points.
    Plots them too.

    Method is to scan a window along the x axis, and at each point in the scan to 
    compute the y percentile within that band.

    Smoothing procedures are implemented.
    """

    #
    # Check inputs
    #

    # TODO need to put checks here

    #
    # Determine scan
    #

    x_scan_points = None
    x_intervals = []

    # Choose between the two different x axis step definition cases
    binned_mode = xbins is not None

    # Check args
    if binned_mode:
        assert (xmin is None) and (xmax is None) and (xwidth is None) and (nx is None)
    else:
        assert use_mean_x_pos is False
        assert (
            (xmin is not None)
            and (xmax is not None)
            and (xwidth is not None)
            and (nx is not None)
        )

    # Now get the scan points/intervals
    if binned_mode:

        # Intervals are the bin edges
        for i in range(xbins.size - 1):
            x_intervals.append((xbins[i], xbins[i + 1]))

        # Scan points are not required here

    else:

        # Get scan points
        x_scan_points = np.linspace(xmin, xmax, num=nx, endpoint=True)

        # Get intervals
        for x_scan_point in x_scan_points:
            this_xwidth = xwidth(x_scan_point) if callable(xwidth) else xwidth
            x_intervals.append(
                (x_scan_point - (this_xwidth / 2.0), x_scan_point + (this_xwidth / 2.0))
            )

    # Checks
    if not binned_mode:
        assert len(x_intervals) == len(x_scan_points)

    #
    # Calculate percentiles
    #

    # Containers to fill
    y_percentile_values = [[] for p in percentiles]
    x_scan_mean_x_pos = []

    # Loop over intervals
    for i_x, xinterval in enumerate(x_intervals):

        # Get data points in this interval
        interval_mask = (x >= xinterval[0]) & (x < xinterval[1])

        # Check if found enough data points to compute a percentile
        if interval_mask.sum() >= min_num_events:

            # Calc y value for each percentile in this x scan point
            for percentile, y_percentile_curve in zip(percentiles, y_percentile_values):
                y_percentile_curve.append(np.percentile(y[interval_mask], percentile))

            # Store the mean x value (in the data) for this scan point
            # TODO median?
            x_scan_mean_x_pos.append(np.mean(x[interval_mask]))

        # Add dummy values if could not compute percentiles
        else:
            for y_percentile_curve in y_percentile_values:
                y_percentile_curve.append(np.NaN)
            x_scan_mean_x_pos.append(np.NaN)

    # Numpy-ify
    for i, y_percentile_curve in enumerate(y_percentile_values):
        y_percentile_values[i] = np.array(y_percentile_curve)
    x_scan_mean_x_pos = np.array(x_scan_mean_x_pos)

    # Masks
    y_percentile_masks = [
        np.isfinite(y_percentile_curve) for y_percentile_curve in y_percentile_values
    ]

    #
    # Smoothing
    #

    if smooth is not None:

        # Take a copy of the unsmoothed version
        y_percentile_values_unsmoothed = y_percentile_values

        if smooth_kw is None:
            smooth_kw = {}

        # Spline smoothing
        if smooth == "interpolate":
            from scipy.interpolate import UnivariateSpline

            y_percentile_values = []
            for y_percentile_curve, percentile_mask in zip(
                y_percentile_values_unsmoothed, percentile_masks
            ):
                spline = UnivariateSpline(
                    curve_x[percentile_mask],
                    y_percentile_curve[percentile_mask],
                    **smooth_kw
                )
                y_percentile_values.append(spline(curve_x))

        # Filter smoothing
        elif smooth == "gaussian_filter":
            from scipy.ndimage import gaussian_filter1d

            y_percentile_values = []
            for y_percentile_curve, percentile_mask in zip(
                y_percentile_values_unsmoothed, percentile_masks
            ):
                y_percentile_values.append(
                    gaussian_filter1d(y_percentile_curve, **smooth_kw)
                )

        else:
            raise Exception("Unrecognised 'smooth' option '%s'" % smooth)

    #
    # Format outputs
    #

    if binned_mode:

        # If binned mode, use the bin edges as the curve points
        # Align the percentiles wiht the bin edges too
        output_x_values = xbins
        output_y_percentile_values = [np.append(p, p[-1]) for p in y_percentile_values]
        output_y_percentile_masks = [np.append(m, m[-1]) for m in y_percentile_masks]

    else:

        # If non-binned mode, just use the values
        # Provide the mean positions if requested
        output_x_values = x_scan_mean_x_pos if use_mean_x_pos else x_scan_points
        output_y_percentile_values = y_percentile_values
        output_y_percentile_masks = y_percentile_masks

    #
    # Plot
    #

    # Check user provided an x to plot on
    if ax is not None:

        # Default color
        if color is None:
            color = "red"

        # # Default linestyle
        # if percentile_linestyles is None :
        #     percentile_linestyles = ["-"] * len(percentiles)

        # Plot curves
        plot_percentiles(
            ax=ax,
            x=output_x_values,
            y_percentiles=output_y_percentile_values,
            y_percentile_masks=output_y_percentile_masks,
            binned=binned_mode,
            percentiles=percentiles,
            percentile_linestyles=percentile_linestyles,
            percentile_labels=percentile_labels,
            color=color,
        )

        # # Plot as curve
        # for percentile,y_percentile_curve,percentile_mask,percentile_linestyle in zip(percentiles,y_percentile_values,percentile_masks,percentile_linestyles) :
        #     if xbins is None :
        #         ax.plot( curve_x[percentile_mask], y_percentile_curve[percentile_mask], color=color, linestyle=percentile_linestyle, label="%0.3g%%"%percentile, zorder=10 )
        #     else :
        #         xplot = xbins[ np.append(percentile_mask, percentile_mask[-1] ) ]
        #         yplot = np.append(y_percentile_curve[percentile_mask], y_percentile_curve[percentile_mask][-1])
        #         ax.step( xplot, yplot, where="post", color=color, linestyle=percentile_linestyle, label="%0.3g%%"%percentile, zorder=10 )

    #
    # Done
    #

    return output_x_values, output_y_percentile_values, output_y_percentile_masks


def plot_percentiles(
    ax,
    x,
    y_percentiles,
    y_percentile_masks,
    percentiles,
    percentile_linestyles=None,
    percentile_labels=None,
    binned=True,
    color=None,
):

    # Default color
    if color is None:
        color = "red"

    # Defaults
    if percentile_linestyles is None:
        percentile_linestyles = ["-"] * len(percentiles)
    if percentile_labels is None:
        percentile_labels = [ "%0.2g%%"%p for p in percentiles ]

    # Plot as curve
    for p, yp, ypm, pls, plb in zip(
        percentiles, y_percentiles, y_percentile_masks, percentile_linestyles, percentile_labels
    ):
        if binned:
            ax.step(
                x[ypm],
                yp[ypm],
                where="post",
                color=color,
                linestyle=pls,
                label=plb,
                zorder=10,
            )
        else:
            ax.plot(
                x[ypm],
                yp[ypm],
                color=color,
                linestyle=pls,
                label=plb,
                zorder=10,
            )
