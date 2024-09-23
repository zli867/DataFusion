import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import colors

def plotPolygons(polygon_list, ax, color, crs=None, linewidth=1.0):
    for current_polygon in polygon_list:
        if current_polygon.geom_type == "MultiPolygon":
            for geom in current_polygon.geoms:
                xs, ys = geom.exterior.xy
                if crs is not None:
                    xs, ys = crs(xs, ys)
                ax.plot(xs, ys, color, linewidth=linewidth)

        else:
            xs, ys = current_polygon.exterior.xy
            if crs is not None:
                xs, ys = crs(xs, ys)
            ax.plot(xs, ys, color, linewidth=linewidth)


def get_conf_intercept(alpha, lr, X, y):
    """
    Returns (1-alpha) 2-sided confidence intervals
    for sklearn.LinearRegression coefficients
    as a pandas DataFrame
    """
    coefs = np.r_[[lr.intercept_], lr.coef_]
    X_aux = np.zeros((X.shape[0], X.shape[1] + 1))
    X_aux[:, 1:] = X
    X_aux[:, 0] = 1
    dof = -np.diff(X_aux.shape)[0]
    mse = np.sum((y - lr.predict(X)) ** 2) / dof
    var_params = np.diag(np.linalg.inv(X_aux.T.dot(X_aux)))
    t_val = stats.t.isf(alpha / 2, dof)
    gap = t_val * np.sqrt(mse * var_params)
    return {
        "coeffs": coefs,
        "lower": coefs - gap,
        "upper": coefs + gap
    }


def discrete_cmap(invertals, base_color_scheme="Spectral_r", low_value_color='withe', high_value_color='purple'):
    cmap = plt.get_cmap(base_color_scheme, invertals)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = colors.to_rgba("white")
    cmap = colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap.set_over(color=high_value_color, alpha=1.0)
    # cmap.set_under(color=low_value_color, alpha=1.0)
    return cmap


def plotComparisonIntercept(X, Y, bin_size, ax, fig):
    # scatter
    #histogram definition
    N = len(X)
    bins = [bin_size, bin_size] # number of bins
    # histogram the data
    hh, locx, locy = np.histogram2d(X, Y, bins=bins)
    # Sort the points by density, so that the densest points are plotted last
    Z = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(X, Y)])
    idx = Z.argsort()
    # normalize z2
    x2, y2, z2 = X[idx], Y[idx], Z[idx] / N
    z2_min = np.min(z2[z2 > 0])
    cax = ax.scatter(x2, y2, c=z2, cmap='jet', marker='.', s=2, norm=colors.LogNorm(vmin=z2_min, vmax=z2.max()))
    divider = make_axes_locatable(ax)
    d_cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(cax, cax=d_cax, extend='both')
    cb.ax.set_title('Density', fontsize=10)
    # reg
    X = X.reshape(-1, 1)
    lin_model = LinearRegression().fit(X, Y)
    line_x = np.linspace(0, np.max(X)).reshape(-1, 1)
    line_y = lin_model.predict(line_x)
    # Confident intervals
    cf = get_conf_intercept(0.05, lin_model, X, Y)
    r2_score = lin_model.score(X, Y)
    reg_txt = 'y={:.2f}x+{:.2f}, $R^2$ = {:.2f}'.format(cf["coeffs"][1], cf["coeffs"][0], r2_score)
    # ax.text(0.8, 0.2, reg_txt, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.plot(line_x, line_y, 'r')
    return reg_txt


def discrete_pos_neg_cmap(invertals, base_color_scheme="seismic", neutral_color='white', low_value_color='black', high_value_color='purple'):
    invertals = (invertals // 2) * 2 + 1
    cmap = plt.get_cmap(base_color_scheme, invertals)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    cmaplist[invertals // 2] = colors.to_rgba(neutral_color)
    cmap.set_over(color=high_value_color, alpha=1.0)
    cmap.set_under(color=low_value_color, alpha=1.0)
    return cmap
