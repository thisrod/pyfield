# methods with no docstring are described in README.md

from stdlib import *
from numpy import savez
from os.path import isfile
from scipy.ndimage.interpolation import map_coordinates

_albls = ['_s', '_x', '_y', '_z']

class Grid:

	#
	# display
	#
		
	def __str__(self):
		return "<grid " + " in R^" + str(self.dim()) + \
			" shape " + str(self.shape) + \
			" over " + str([tuple(b) for b in self.bounds()]) + ">"
			
	def __repr__(self):
		return str(self)
		
	#
	# construction
	#
	
	# the method __init__ is responsible for casting instance variables to the right type.
		
	@classmethod
	def from_axes(cls, *axes):
		axes = [array(q).flatten() for q in axes]
		N = array([q.size for q in axes])
		assert (N>=2).all()
		origin = [q[0] for q in axes]
		return cls(
			shape=N,
			h=[ptp(q)/(n-1) for n, q in zip(N, axes)],
			U=eye(len(N)),
			p=origin,
			o=origin)
			
	@classmethod
	def delta(cls, x, h):
		return cls(shape=(1,), h=[h], U=[1], p=[x], o=[x])
	
	@classmethod
	def default(cls):
		ffile = load('fields.npz')
		return cls.from_axes(*[ffile[q] for q in _albls])

	def __init__(self, shape, o, p, h, U):
		# instance variables explained in geometry.tex
		self.U = array(U, dtype=float)
		assert allclose(dot(self.U.T, self.U), eye(self.rank()))
		self.shape = tuple(int(n) for n in shape)
		assert len(self.shape) == self.rank()
		assert all(n>0 for n in self.shape)
		self.h = array(h, dtype=float).reshape((self.rank(),))
		assert (self.h>0).all()
		self.p = array(p, dtype=float).reshape((self.rank(),))
		self.o = array(o, dtype=float).reshape((self.dim(),))
	
	#
	# basic properties
	#
	
	def rank(self):
		"return the rank of arrays of samples"
		return self.U.shape[1]
		
	def dim(self):
		"return the dimension of the common coordinate space"
		return self.U.shape[0]
		
	def axes(self):
		return [p+h*arange(n) for p, h, n in zip(self.p, self.h, self.shape)]
		
	def bounds(self):
		n = array(self.shape)
		return array([p-0.5*h, p+(n+0.5)*h])

	def __len__(self):
		return prod(self.shape)


	#
	# utility
	# 
	
	def _vectors(self, A):
		"return h, p,and  o, reshaped to broadcast over A"
		tail = (1,)*(A.ndim-1)
		return self.h.reshape((self.rank(),)+tail), \
			self.p.reshape((self.rank(),)+tail), \
			self.o.reshape((self.dim(),)+tail)
			
	def _clone(self, **args):
		"return a Grid like self, but with some variables changed"
		d = dict([(x, getattr(self, x)) for x in ['shape', 'h', 'o', 'p', 'U']])
		d.update(args)
		return Grid(**d)
	
	def _trunc(self, axs):
		"return a subgrid of self, projected onto the axes numbered in axs"
		return self._clone(shape=[self.shape[i] for i in axs],
			h = self.h[axs],
			p = self.p[axs],
			U = self.U[:,axs])


	#
	# projection and cartesian products
	#
			
	def __mul__(self, other):
		shape = self.shape + other.shape
		h = concatenate((self.h, other.h))
		p = concatenate((self.p, other.p))
		if self.dim() == other.dim() and \
			allclose(dot(self.U.T, other.U), zeros((self.rank(), other.rank()))):
			return Grid(shape=shape, h=h, p=p,
				o = dot(self.U, dot(self.U.T, self.o)) + \
					dot(other.U, dot(other.U.T, other.o)),
				U = concatenate((self.U, other.U), axis=1))
		else:
			U = zeros((self.dim()+other.dim(), self.rank()+other.rank()))
			U[:self.dim(),:self.rank()] = self.U
			U[self.dim():,self.rank():] = other.U
			return Grid(shape=shape, h=h, p=p,
				o=concatenate((self.o, other.o)),
				U=U)
				
	def __getitem__(self, ix):
	
		def bound(i, n, m):
			if i is None:
				i = m
			if i<0:
				i += n
			assert 0 <= i and i <= n
			return i
			
		if type(ix) is not tuple:
			ix = (ix,)
		assert len(ix) == self.rank()	# numpy missing indices n.y.i.
		
		# is there a simple test whether ix[x] has an integer type?
		daxs = [i for i in range(self.rank()) if
			ix[i] is not None and type(ix[i]) is not slice]
		pt = array([ix[i] for i in daxs])
		D = self._trunc(daxs)
		D = D._clone(shape=(1,)*len(daxs), p=D.w(i=pt), o=D.W(i=pt))
		
		naxs = [i for i in range(self.rank()) if type(ix[i]) is slice]
		assert all([ix[i].step is None for i in naxs])
		low = array([bound(ix[i].start, self.shape[i], 0) for i in naxs])
		high = array([bound(ix[i].stop, self.shape[i], self.shape[i]) for i in naxs])
		N = self._trunc(naxs)
		N = N._clone(shape=high - low, p=N.w(i=low), o=N.W(i=low))
		
		return D*N
				
			
	#
	# indices and coordinates
	#
		
	def w(self, i=None, W=None):
		assert i is None or W is None
		if i is not None:
			i = array(i)
			h, p, o = self._vectors(i)
			return p + h*i
		if W is not None:
			W = array(W)
			h, p, o = self._vectors(W)
			return p + tensordot(self.U.T, W-o, 1)		
		else:
			return Field(self.w(i=self.i()), self)
		
	def W(self, i=None, w=None):
		assert i is None or w is None
		if i is not None:
			i = array(i)
			h, p, o = self._vectors(i)
			return o + tensordot(self.U, h*i, 1)
		if w is not None:
			w = array(w)
			h, p, o = self._vectors(w)
			return o + tensordot(self.U, w-p, 1)
		else:
			return Field(self.W(i=self.i()), self)
		
	def i(self, w=None, W=None):
		assert w is None or W is None
		if w is not None:
			w = array(w)
			h, p, o = self._vectors(w)
			return (w-p)/h
		if W is not None:
			W = array(W)
			h, p, o = self._vectors(W)
			return tensordot(self.U.T, W-o, 1)/h
		else:
			return Field(indices(self.shape), self)
		
	def ww(self, w=None):
		"inertia tensors about the grid origin for a unit mass at each point of r, or the grid by default.  unlear if this is useful for rank greater than three."
		if w is None:
			w = self.w()
		else:
			w = ascontiguousarray(w)
		P = -array(w)[newaxis,::]*array(w)[:,newaxis,::]
		# hack strides to address the diagonal plane of P
		Pd = ndarray(buffer=P, dtype=P.dtype, \
			shape=w.shape, strides=((3+1)*P.strides[1],)+ P.strides[2:])
		Pd += (w**2).sum(0)
		if type(w) is Field:
			return Field(P, w.abscissae)
		else:
			return P
		
	
	#
	# transformation
	#
	
	def shifted(self, new_origin):
		return Grid(p=self.p-new_origin, o=self.o, h=self.h, shape=self.shape, U=self.U)
		
	def translated(self, dr):
		assert False	# n.y.i.
		
	def rotated(self, V, centre=None):
		if centre is None:
			centre = zeros(self.rank())
		assert V.shape == 2*(self.dim(),)
		Rc = self.W(w=centre)
		return self._clone(o=Rc+dot(V, self.o-Rc), U=dot(V,self.U))
	
	# loading and saving to file
		
	def be_default(self):
		assert False	# bit rot
		ftab = dict(load('fields.npz'))
		ftab.update(dict(zip(_albls, self.axes)))
		savez('fields.npz', **ftab)
		
	#
	# comparision
	#
		
	def __eq__(self, other):
		assert False	# testing floats for equality is a sin
		
	def __neq__(self, other):
		return not self.__eq__(other)
		
	def close_points(self, other):
		"test if the points of self and other have similar common coordinates"
		return self.shape == other.shape and \
			all([allclose(getattr(self, x), getattr(other, x))
				for x in ['h', 'o', 'U']])
	
	def close_grid(self, other):
		"test if the points of self and other have similar grid coordinates"
		return self.shape == other.shape and \
			all([allclose(getattr(self, x), getattr(other, x))
				for x in ['h', 'p']])
		
	def close(self, other):
		return self.close_points(other) and self.close_grid(other)
		
	#
	# Fourierology
	#
	
	def reciprocal(self):
		"""wavenumbers at which dfft is sampled
		
the reciprocal grid size is always odd.  has same orientation as self
"""
		ns = array([n+(n+1)%2 for n in self.shape])
		h=2*pi/(self.shape*self.h)
		return Grid(shape=ns, h=h, U=self.U, p=-h*((ns-1)//2), o=zeros(self.dim()))
		
	def kspace(self):
		"""quick and dirty wavenumbers for q&d spectral methods"""
		
		return array(meshgrid(
			*[2*pi*fftfreq(q.size, ptp(q)/(q.size-1)) for q in self.axes()],
			indexing='ij'))
		
	#
	# integration
	#
		
	def S(self, ordinates):
		# integrate: should allow axes to be picked
		# the method ndarray.sum returns an ndarray, not a Field,
		# so this doesn't trigger shape checks.
		return prod(self.h)*ordinates.sum(tuple(range(-self.rank(),0)))
		
	#
	# misc
	#
		
	def blank(self):
		return Field(empty(self.shape), self)
		
	def cover(self, other):
		"return a grid like me, with extra axes prepended to cover the bounds of other"
		# for the axes of other: if it has a significant orthogonal component, prepend it with step geometric_mean(other.h) and origin derived from projecting the axes of other onto it.
		assert False
		
		
	
# For now, load and store read and write the whole file every time a
# saved field changes. The other way would be a singleton object that
# cached the file as a dict, and automatically saved the fields on exit; doing that quickly is beyond my Python
# ability.

def load_field(lbl):
	ffile = load('fields.npz')
	return Field(ffile[lbl], label=lbl)
	

class SampledField(ndarray):
	
	# FIXME binary ufuncs should check the grids are the same
	
	def __new__(cls, ordinates, abscissae=Grid.default(), label=None):
		obj = asarray(ordinates).view(cls)
		obj.abscissae = abscissae
		assert obj.shape[-obj.abscissae.rank():] == obj.abscissae.shape
		if label is not None:
			obj.label = label
		return obj
				
	def __array_finalize__(self, obj):
		if obj is None: return
		self.abscissae = getattr(obj, 'abscissae', None)
		self.coords = getattr(obj, 'coords', None)
		if self.abscissae:
			assert self.shape[-self.abscissae.rank():] == self.abscissae.shape
			
	def __str__(self):
		return '<not yet implemented>'
		
	def __repr__(self):
		return '<not yet implemented>'

	def __getitem__(self, ixs):
		if any([type(i) is slice for i in ixs]):
			return Field(self.view(ndarray).__getitem__(ixs),
				self.abscissae.__getitem__(ixs))
		else:
			return(self.view(ndarray).__getitem__(ixs))
	
	#
	# grid methods
	#
	
	def r(self):
		return self.abscissae.w()
	def R(self):
		return self.abscissae.W()
	def i(self):
		return self.abscissae.i()
	def rr(self):
		return self.abscissae.ww()
	
	def S(self):
		"should allow dimensions to be specified"
		return self.abscissae.S(self)
	
	#
	# loading and saving
	#
	
	def save(self, label=None):
		if label:
			self._label = label
		assert self.abscissae == Grid.default()
		ftab = dict(load('fields.npz'))
		ftab[self._label] = self.ordinates
		savez('fields.npz', **ftab)
		
	#
	# Fourier transforms
	#
	
	def fft(self):
		assert False	# n.y.i.
		ft = fftn(self.view(ndarray))
		ft = fftshift(ft)
		if len(ft) % 2 == 0:
			ft = concatenate((ft[:1]/2, ft[1:], ft[:1]/2))
		# document why this satisfies Rayleigh's theorem
		ft /= sqrt(2*pi)
		F = Field(ft, self.abscissae.reciprocal())
		F.coords = self.abscissae
		return F
	
	def ifft(self):
		return Field(ifftn(self.view(ndarray)), self.abscissae.reciprocal())
		
	#
	# differentiation
	#
	
	def expdsq(self, l):
		"""return exp(l*del^2) self.
		
minimum useful spectral method.
"""

		k = self.abscissae.kspace()
		dvec = exp(l*(k**2).sum(axis=0))
		return Field(ifftn(dvec*fftn(self.view(ndarray))), self.abscissae)


	#
	# interpolation
	#
	
	def sampled(self, abscissae):
		
		# extend my delta axes to two equal samples, at ±h/2
		eself = array(self)
		d = zeros(self.ndim)
		for i in range(self.ndim):
			if eself.shape[i] == 1:
				eself = concatenate((eself, eself), axis=i)
			d[i] = -self.abscissae.h[i]/2
		egrid = self.abscissae._clone(shape=eself.shape,
			p=self.abscissae.p - d,
			o=self.abscissae.o - d)
			
		# interpolate me on an extension of abscissae
		S = abscissae.cover(self)
		f = map_coordinates(eself, egrid.i(W=S.W()), cval=nan)
			
		# if any value along an epsilon axis was not extrapolated, 
		# fill in the nans with zeros.
		eaxs = tuple(range(S.rank()-abscissae.rank()))
		inbnds = logical_not(isnan(f).all(axis=eaxs))
		f[isnan(f) & inbnds] = 0
		
		# integrate over epsilon axes
		return Field(f.sum(axis=eaxs), abscissae)
				
	def setsamples(self, fld):
		# what should happen if some of my epsilon axes
		# are normal or delta axes in fld?
		f = fld.sampled(self.abscissae)
		# copy non-nans back (but what if the nan was interpolated?)
		ix = logical_not(isnan(f))
		self[ix] = f[ix]


	#
	# misc
	#

	def support(self, cut=0):
		"return a subgrid of abscissae, on which ordinates exceed cutoff"
		A = array(self)
		A[isnan(A)] = -inf
		bools = array(nonzero(A>cut))
		return self.abscissae[tuple(slice(m,n) for m, n in zip(bools.min(axis=1), bools.max(axis=1) + 1))]
