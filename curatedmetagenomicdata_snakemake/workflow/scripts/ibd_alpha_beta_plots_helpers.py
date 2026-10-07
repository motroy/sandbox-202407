"""Shared plotting helpers for the IBD figures."""
import numpy as np
from matplotlib.patches import Ellipse


def draw_ellipse(df, group_name, ax, color):
    """Draw a 1-SD covariance ellipse (first two columns) for one condition."""
    g = df[df["condition"] == group_name].iloc[:, :2]
    mean, cov = g.mean(axis=0).to_numpy(), g.cov().to_numpy()
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * np.sqrt(vals)
    ax.add_patch(Ellipse(xy=mean, width=width, height=height, angle=theta, color=color, alpha=0.2))
