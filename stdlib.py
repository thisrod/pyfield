# Standard imports for an environment like Matlab

# This doesn't import any interpolation routines: Python keeps reinventing that wheel, going from square to pentagonal to triangular.  See scipy.ndimage.interpolation and scipy.interpolate

from numpy import \
	abs, amax, amin, all, allclose, arange, argmax, around, array, asarray, ascontiguousarray, \
	ceil, concatenate, diag, diagflat, diagonal, dot, empty, eye, \
	floor, indices, inf, isnan, linspace, load, log, log10, logspace, \
	logical_not, logical_or, \
	mean, meshgrid, nan, ndarray, newaxis, nonzero, \
	ones, pi, prod, ptp, \
	sign, sqrt, squeeze, tensordot, zeros, \
	exp, sin, cos, tan
from numpy.random import rand, randn, seed
from numpy.fft import fft, ifft, fftn, ifftn, fftfreq, fftshift
from numpy.linalg import eig, norm, qr

# work around Pyenv breaking the OSX backend for Matplotlib
import matplotlib; matplotlib.use("Agg")
from matplotlib.pyplot import colorbar, contour, contourf, figure, imshow, legend, loglog, plot, savefig, semilogy, subplot, text, title, xlabel, ylabel, xlim, ylim, xticks, yticks

def sciform(xs):
	"return a collection like xs, with numeric elements replaced by strings.  Also an exponent term"
	
	# bug: sciform((5, 9.7)) prints "5", "10"
	# better way: shift so smallest number is o(1)
	#	round
	#	shift so largest number is o(1)
	#	print so that smallest number has one s.f.
	
	M = int(floor(max(log10(abs(x)) for x in xs)))
	m = int(floor(min(log10(abs(x)) for x in xs)))
	return [("%." + str(M-m) + "f") % (x/10**M) for x in xs], \
		"" if M == 0 else "*1e" + str(M)

seed(666)		# repeatability