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
		return "<grid " + " in R^" + str(self.dim) + \
			" shape " + str(self.shape) + \
			" over " + str([tuple(b) for b in self.bounds]) + ">"
			
	def __repr__(self):
		return str(self)
		
	#
	# construction
	#
	
	# the method __init__ is responsible for casting instance variables to the right type.

	def __new__(cls, shape, o, p, h, U):
		if shape == ():
			return object.__new__(NullGrid)
		else:
			return object.__new__(cls)
		
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
	def delta(cls, x, h=0):
		return cls(shape=(1,), h=[h], U=[[1]], p=[x], o=[x])
	
	@classmethod
	def default(cls):
		ffile = load('fields.npz')
		return cls.from_axes(*[ffile[q] for q in _albls])

	def __init__(self, shape, o, p, h, U):
		# instance variables explained in geometry.tex
		self.U = array(U, dtype=float)
		r = self.rank
		assert allclose(dot(self.U.T, self.U), eye(r))
		N = array(shape, dtype=int)
		self.shape = tuple(N)
		assert N.size == r
		assert (N > 0).all()
		self.h = array(h, dtype=float).reshape((r,))
		assert (self.h >= 0).all()
		assert logical_or(N == 1, self.h > 0).all()
		self.p = array(p, dtype=float).reshape((r,))
		self.o = array(o, dtype=float).reshape((self.dim,))
	
	#
	# basic properties
	#
	
	def __getattr__(self, name):
		if name in Grid._getters:
			return Grid._getters[name](self)
		else:
			assert False
	
	_getters = {
		'rank': lambda self: self.U.shape[1],
		'dim': lambda self: self.U.shape[0],
		'axes': lambda self: [Grid(shape=self.shape[i:i+1], o=self.o,
			h=self.h[i:i+1], p=self.p[i:i+1], U=self.U[:,i:i+1])
			for i in range(self.rank)],
		'bounds': lambda self: array([self.p - 0.5*self.h, 
			self.p + (array(self.shape)-0.5)*self.h]) }

	def __len__(self):
		return prod(self.shape)


	#
	# utility
	# 
	
	def _vectors(self, A):
		"return h, p,and  o, reshaped to broadcast over A"
		tail = (1,)*(A.ndim-1)
		return self.h.reshape((self.rank,)+tail), \
			self.p.reshape((self.rank,)+tail), \
			self.o.reshape((self.dim,)+tail)
			
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
	# along and through
	#
	
	def along(self, other):
		return self._clone(o=other.o, U=other.U[:,:self.rank])
		
	def through(self, x):
		n = self.dim-self.rank
		U, R = qr(concatenate((self.U,eye(self.dim)), axis=1))
		return self * \
			Grid.from_axes(*([[0]]*n)).rotated(U[:, self.rank:]).translated(x)


	#
	# projection and cartesian products
	#
			
	def __mul__(self, other):
		shape = self.shape + other.shape
		h = concatenate((self.h, other.h))
		p = concatenate((self.p, other.p))
		if self.dim == other.dim and \
			allclose(dot(self.U.T, other.U), zeros((self.rank, other.rank))):
			return Grid(shape=shape, h=h, p=p,
				o = dot(self.U, dot(self.U.T, self.o)) + \
					dot(other.U, dot(other.U.T, other.o)),
				U = concatenate((self.U, other.U), axis=1))
		else:
			U = zeros((self.dim+other.dim, self.rank+other.rank))
			U[:self.dim,:self.rank] = self.U
			U[self.dim:,self.rank:] = other.U
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
		assert len(ix) == self.rank	# numpy missing indices n.y.i.
		
		# build the result one axis at a time
		section = Grid(shape=(), o=self.o, p=[], h=[], U=empty((self.dim,0)))
		for q, i in zip(self.axes, ix):
			if type(i) is slice:
				assert i.step is None
				low = bound(i.start, len(q), 0)
				high = bound(i.stop, len(q), len(q))
				section = section*q._clone(
					shape=(high - low,), 
					o=q.o+low*q.U[:,0]*q.h[0],
					p=q.p+low*q.h)
					
		return section
				
			
	#
	# indices and coordinates
	#
	
	# special case for null grids
		
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
			return SampledField(self.w(i=self.i()), self)
		
	def W(self, i=None, w=None):
		assert i is None or w is None
		assert self.rank == self.dim	# points are columns otherwise
		if i is not None:
			i = array(i)
			h, p, o = self._vectors(i)
			return o + tensordot(self.U, h*i, 1)
		if w is not None:
			w = array(w)
			h, p, o = self._vectors(w)
			return o + tensordot(self.U, w-p, 1)
		else:
			return SampledField(self.W(i=self.i()), self)
		
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
			return SampledField(indices(self.shape), self)
		
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
		if type(w) is SampledField:
			return SampledField(P, w.abscissae)
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
			centre = zeros(self.rank)
		assert V.shape == 2*(self.dim,)
		Rc = self.W(w=centre)
		return self._clone(o=Rc+dot(V, self.o-Rc), U=dot(V,self.U))
	
	# loading and saving to file
		
	def be_default(self):
		assert False	# bit rot
		ftab = dict(load('fields.npz'))
		ftab.update(dict(zip(_albls, self.axes)))
		savez('fields.npz', **ftab)
		
	#
	# comparison
	#
	
	def spans(self, other):
		assert self.dim == other.dim
		if other.rank == 1:	# meshgrid barfs
			box = other.bounds
		else:
			box = meshgrid(*other.bounds.T, indexing='ij')
		vertices = other.W(array(box))
		P = dot(self.U, self.U.T)
		return allclose(tensordot(P, vertices, axes=1), vertices)
		
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
		return Grid(shape=ns, h=h, U=self.U, p=-h*((ns-1)//2), o=zeros(self.dim))
		
	def kspace(self):
		"""quick and dirty wavenumbers for q&d spectral methods"""
		
		return array(meshgrid(
			*[2*pi*fftfreq(len(q), q.h[0]) for q in self.axes],
			indexing='ij'))
		
	#
	# integration
	#
		
	def S(self, ordinates):
		# integrate: should allow axes to be picked
		# the method ndarray.sum returns an ndarray, not a Field,
		# so this doesn't trigger shape checks.
		return prod(self.h)*ordinates.sum(tuple(range(-self.rank,0)))
		
	#
	# misc
	#
		
	def blank(self):
		return SampledField(empty(self.shape), self)
		
	def cover(self, other):
		"return a grid like me, with extra axes prepended to cover the bounds of other"
		# for the axes of other: if it has a significant orthogonal component, prepend it with h as the maximum of other.h, size the sum of other.size, and move my grid origin along the new axis until my centre aligns with other's grid origin.
		# this is bodgy
		h = other.h.max()
		n = sum(other.size)
		U = self.U
		o = self.o
		xaxs = 0
		for i in range(other.rank):
			prp = dot(U, dot(U.T, other.U[:,i]))
			if not allclose(prp, zeros(self.dim)):
				prp = prp/norm(prp)
				U = concatenate((prp.reshape((-1,1)), U), axis=1)
				o += (dot(prp, other.o-o)-h*n)*prp
				xaxs += 1
		return self._clone(shape=(n,)*xaxs+self.shape,
			h=concatenate((array([h]*xaxs), self.h)),
			p=concatenate((zeros(xaxs), self.p)),
			o=o,
			U=U)
			

class NullGrid(Grid):
	"""this overrides methods where numpy fails to treat shape 0 arrays the way we want.
	
a field over this is constant everywhere.  this is consistent with Numpy, where zero-dimensional arrays have a single value.
"""
	
	def __init__(self, shape, o, p, h, U):
		# this is called when Grid.__new__ returns a NullGrid.
		self.shape = tuple(shape)
		self.h = array(h)
		self.p = array(p)
		self.o = array(o)
		self.U = array(U)
		assert self.h.size == 0
		assert self.p.size == 0
		assert self.U.shape[1] == 0
		assert self.U.shape[0] == self.o.size
		
	def __str__(self):
		return '<null Grid>'
		
	def W(self):
		assert False	# points aren't defined
		
		
	
# For now, load and store read and write the whole file every time a
# saved field changes. The other way would be a singleton object that
# cached the file as a dict, and automatically saved the fields on exit; doing that quickly is beyond my Python
# ability.

def load_field(lbl):
	ffile = load('fields.npz')
	return SampledField(ffile[lbl], label=lbl)
	

class SampledField(ndarray):
	
	# FIXME broadcasting
	
	def __new__(cls, ordinates, abscissae, label=None):
		obj = asarray(ordinates).view(cls)
		obj.abscissae = abscissae
		assert type(abscissae) is NullGrid or \
			obj.shape[-obj.abscissae.rank:] == obj.abscissae.shape
		if label is not None:
			obj.label = label
		return obj
				
	def __array_finalize__(self, obj):
		if obj is None: return
		self.abscissae = getattr(obj, 'abscissae', None)
		self.coords = getattr(obj, 'coords', None)
		if self.abscissae:
			assert type(self.abscissae) is NullGrid or \
			self.shape[-self.abscissae.rank:] == self.abscissae.shape
			
	def __str__(self):
		return '<not yet implemented>'
		
	def __repr__(self):
		return '<not yet implemented>'

	def __getitem__(self, ixs):
		if any([type(i) is slice for i in ixs]):
			return SampledField(self.view(ndarray).__getitem__(ixs),
				self.abscissae.__getitem__(ixs))
		else:
			return(self.view(ndarray).__getitem__(ixs))
	
	#
	# abscissae
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
		F = SampledField(ft, self.abscissae.reciprocal())
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
		return SampledField(ifftn(dvec*fftn(self.view(ndarray))), self.abscissae)


	#
	# interpolation
	#
	
	def sampled(self, abscissae):

		# this only implements the cases we need. the full
		# version, that integrates over oblique columns in a
		# rectangular grid, is hard. integrating the
		# band-limited interpolant might be feasible.
		
		if self.abscissae.rank == self.abscissae.dim and \
			self.abscissae.dim == abscissae.dim and \
			abscissae.rank == abscissae.dim:
			return self._sampled_normal(abscissae)
		else:
			return self._sampled_special(abscissae)
			
	def _sampled_normal(self, abscissae):
	
		return SampledField(
			map_coordinates(self, self.abscissae.i(W=abscissae.W()),
				cval=nan),
			abscissae)
			
	def _sampled_special(self, abscissae):
		
		# this only handles the normal -> 1e case
		S =  self.abscissae
		assert S.rank == S.dim
		assert S.dim == abscissae.dim
		assert abscissae.dim - abscissae.rank == 1
		
		# sample with my axis in place of the missing one
		
		eax = argmin([norm(x) for x in dot(self.U.T, abscissae.U)])
		eabsc = Grid(shape=S.shape[eax:eax+1],
			U=S.U[:,eax:eax+1],
			h=S.h[eax:eax+1],
			p=S.p[eax:eax+1],
			o=S.o) * abscissae
		f = self.sampled(eabsc)
		
		# integrate over the missing axis
		
		return SampledField(S.h[eax]*array(f).sum(axis=0), abscissae)
				
	def setsamples(self, fld):
	
		# this only handles assigning delta slices on the first axis
		S, R = self.abscissae, fld.abscissae
		assert S.rank == 4
		assert S.dim == 4
		assert allclose(S.U, R.U)
		assert fld.shape[0] == 1
		assert fld.shape.count(1) == 1
		u = S.U[:,0]
		assert allclose(dot(S.U[:,1:].T, S.o), dot(S.U[:,1:].T, R.o))
		assert allclose(S.h[1:], R.h[1:])
		
		# find index range we're looking for
		I = 1 + int((dot(u, R.o-S.o)+R.h[0]/2)/S.h[0])
		i = 1 + int((dot(u, R.o-S.o)-R.h[0]/2)/S.h[0])
		if I == i+1:
			self[i:I,:,:,:] = fld

	
	#
	# plotting
	#
	
	def section_positive(self):
		self._section('gray')
	
	def section_negative(self):
		self._section('gray_r')
		
	def _section(self, cm):
		x, y = [q for q in self.abscissae.axes if len(q)>1]
		window =  (x*y).bounds.T.flatten()
		imshow(squeeze(array(self)).T, interpolation='nearest',
			origin='lower',
			extent=tuple(window), cmap=get_cmap(cm))
		

	#
	# misc
	#

	def support(self, cut=0):
		"return a subgrid of abscissae, on which ordinates exceed cutoff"
		A = array(self)
		A[isnan(A)] = -inf
		bools = array(nonzero(A>cut))
		return self.abscissae[tuple(slice(m,n) for m, n in zip(bools.min(axis=1), bools.max(axis=1) + 1))]
