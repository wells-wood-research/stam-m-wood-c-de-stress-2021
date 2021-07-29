# 0. Loading packages--------------------------------------------------------------
import numpy as np
from matplotlib.patches import Ellipse
from scipy.stats import multivariate_normal


def get_cov_ellipse(cov, centre, nstd, **kwargs):
    """
    Return a matplotlib Ellipse patch representing the covariance matrix
    cov centred at centre and scaled by the factor nstd.

    """

    # Find and sort eigenvalues and eigenvectors into descending order
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # The anti-clockwise angle to rotate our ellipse by
    vx, vy = eigvecs[:, 0][0], eigvecs[:, 0][1]
    theta = np.arctan2(vy, vx)

    # Width and height of ellipse to draw
    width, height = 2 * nstd * np.sqrt(eigvals)
    return Ellipse(
        xy=centre, width=width, height=height, angle=np.degrees(theta), **kwargs
    )


def plot_cov_ellipses(data, nstd, **kwargs):
    """
    Plot matplotlib Ellipses patch representing the covariance matrix
    cov centred at centre and scaled by the factor nstd.

    """

    # Calculating the mean
    mean = np.mean(data[["pca_dim0", "pca_dim1"]], axis=0).tolist()

    # Calculating the covariance
    cov = np.cov(data[["pca_dim0", "pca_dim1"]], rowvar=0)

    # Fitting the gaussian
    # normal_dist = multivariate_normal(mean=mean, cov=cov)

    # Plotting ellipse
    plot = get_cov_ellipse(
        cov=cov,
        centre=(mean[0], mean[1]),
        nstd=nstd,
        fc=None,
        linestyle="--",
        linewidth=1,
        color="black",
        fill=False,
        alpha=0.7,
    )

    return plot
