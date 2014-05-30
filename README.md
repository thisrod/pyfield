Field library
===

by Rodney E. S. Polkinghorne

A field is a scalar, vector or tensor quantity that depends on position and time.  This Python library exports a type `Field`, an `ndarray` containing samples of an array-valued quantity, that also records the points at which the samples were taken.  These points form a rectangular grid in R<sup>n</sup>.  Such an array can integrate, differentiate and Fourier transform itself, generate samples of white noise at its points, and so on.

The library aims to remove the accidental complexity from computing with fields, in order to spare scientific programers from bookeeping, to prevent large classes of bugs from occuring at all, and to allow the remaining code to address physical problems, and the remaining bugs to be removed by physical nous.

The current version of the library assumes that the quantity being sampled is an array, whose components do not vary as the coordinates change underneath them.  Future versions might understand how the components of vectors and tensors should transform in different coordinate systems.  

The library code is a draft.  Many user actions that ought to generate a meaningful exception instead cause an assertion to fail, and sometimes that assertion has a comment "n.y.i.", meaning not yet implemented.  However, this document should describe exactly what the user is allowed to do, and what the library should do in response, and the code should either do that or fail deliberately.  Please report exceptions to that.


Grids
---

The mechanics of recording coordinates are done by a class `Grid`, representing a rectangular grid of points in R<sup>n</sup>, along with a system of coordinates.  The grid need not be aligned with the usual axes or start at the origin: we can represent a skew plane in space, and things like that.

At `Grid` is usually constructed from its axes, as

	S = Grid.from_axes(0.1*arange(3), pi+0.5*arange(4))

The `*` operator on `Grid` is a cartesian product, so the same grid can be constructed with

	x, y = Grid(0.1*arange(10)), Grid(pi+0.5*arange(20))
	S = x*y
	
However, the points of `S` now lie in the plane, while `x` and `y` comprise points on the line, with no record of which corresponds to the first axis of `S` and which to the second.  We'd be better off constructing `S` the first way, and projecting out the axes with

	x, y = S.axes()
	
Now `x` and `y` lie in the plane, and `x*y`, `y*x` and `S` are the same.

FIXME describe the precise rules for * later on

FIXME You can't do S.W() on a low-rank grid, instead do do S.through(point).W()

An array of coodinates for all grid points in R<sup>n</sup> is given by

	S.W()

and coordinates for the point with a known index, say `(1, 2)`, by

	S.W(i=(1, 2))

A `Grid` constructed like this is parallel to the axes of R<sup>n</sup>, but more general grids can be derived by rotation, specified by an orthogonal matrix.

	U = array([[cos(1), sin(1)], [-sin(1), cos(1)]])
	T = S.rotated(U)

As well as rotating grids, we can translate them.

	S.translated((0, -pi))

is a grid similar to `S`, but starting from the origin of the plane.

TODO Implement the full set of MetaFont transforms.

TODO Grids in the Argand plane, complex h, shape is a Gaussian integer

As well as the common coordinates returned by `W()`, each grid has a system of grid coordinates, returned by `w()`.  These are preserved by translations and rotations.

	S.w()
	S.rotated(U).w()
	S.translated((0, -pi)).w()

However, the method `shifted` translates the origin of the grid coordinates

	S.shifted((1,1)).w()

The methods `W` and `w` can take an array of indices or of coordinates

	allclose(S.w(), S.w(i=indices(S.shape)))
	allclose(S.w(), S.w(W=S.W()))

The method `i` calculates indices from coordinates

	S.i(w=S.w())
	S.i(W=S.W())

These methods return a `Field` if called with no arguments or with a `Field`, and an `ndarray` if passed an `ndarray`.

Grids and Fields have bounds, that extend half a grid step past the extreme points at each edge.

	one.bounds()

This has shape 2*one.rank().


Fields
---

A `Field` is constructed from a `Grid` and a set of data sampled on it, as follows

	one = Field(ones(S.shape), S)

This is an `ndarray`, so we can do things like

	allclose(one + exp(1j*pi*one), 0*one)

In fact, the array `S.W()`, which can also be expressed as `one.R()`, is a field on `S`, so

	r0 = (one.R()*one).S()/one.S()

is a complicated way to find the centre of mass of a rectangle.  The method `S` computes the integral of a `Field`.

When a field has more than one component, we need a convention about the order that the axes go in when the samples are represented as an array.  The convention is

other components*axis components*sample points.

So, for example, samples of the gradient of the wavefunction of a spin 1/2 particle on a 4*5*6 grid would be stored in an array of shape 2*3*4*5*6.  The 2 is the spin components, and the 3 is the x, y and z components of the gradient vector.

Fields can be interpolated on other grids

	a = one.sampled(T)

Some points in `a` and `b` have been extrapolated; these samples have the value `nan`.  The corresponding assignment operation is `setsamples`

	a[:,:] = 2
	one.setsamples(a)
	
The values of `one` outside the bounds of `T` are nearly preserved.  The actual algorithm for setsamples is as follows.  It samples the lvalue on the grid of the rvalue, subtracts the result from the rvalue samples, interpolates that on the lvalue, and adds that to the lvalue samples.  In the case that the grids coincide, this reduces to assignment.

Does setting samples on disjoint parts of a field commute?

We can get an array of results by adding slices at close times to a field with a larger step.  Does assigning the slices to the field do the same thing?


Degenerate grids
---

So far, we have seen grids that cover a rectangle in the plane, a cube in space, or generally an n-prism in n-space.  The library is more general than this.  However, there are two special cases we have not considered.  These are a grid whose dimension is greater than its rank, and a grid that has only one point along some of its axes.  Or both: these are not mutally exclusive.

A grid with only one point along an axis is treated fairly simply.  A field over such a grid is sampled normally along the axes with multiple points.  On the degenerate axis, is varies as sinc(2&pi;x/h), with the first zeros at &pm;h.  Even with one point, the grid still has a spacing, for this purpose.  If this field is sampled on a nondegenerate grid, the sinc function will actually be sampled.  This means that taking slices with spacing h, on grids of width h, then adding them back together, will reconstruct the bandwidth limited interpolant of the field that was sliced.

Such grids can be constructed by subscription, e.g., S[2:3, :].  A subscription S[2, :] would create a low-rank grid: subscripting a field in this way constructs a constant field whose values are a slice instead of an integral.

No degenerate grid can be constructed using `Grid.from_axes`.  There is no way to infer a reasonable step from an axis with a single point, and this is required even in degenerate grids.  However, `Grid.delta(x, h)` constructs a 1D delta grid at coordinate x with step h.

A grid with rank less than its dimension is constructed by supplying `None` as a subscript.  We already saw this when `x` and `y` were constructed from `S`.  A field over this is treated as a constant over the missing axes in the common space.  For example, it might be a function of space in a dynamical simulation.  Another example is a field whose value is one of the grid coordinates, which is constant over the other coordinates.

There are some rules about such grids and sampling.  These ensure that delta grids and low-rank grids are round trip compatible: if we sample from one to the other and back again, we get the same field.  The rules when a delta field or a constant field is sampled are simple: they're treated as a delta function or a constant.  It is an error to sample a delta function on a grid that is not degenerate over the delta axes.  When a delta field is sampled on a delta grid, the interpolated value is the field value if the grid lies within the bounds of the field-delta grids still have a step along the degenerate axis, and the bounds are computed as usual, a box of width h about the single point.  Otherwise, the interpolated value is zero.  When an ordinary field is sample on a delta grid, a section through the field is interpolated.

Sampling on a low-rank grid integrates over missing axes.

Sampling from a low-rank grid assumes constant over missing axes.

Subscripting a field causes it to be interpolated on the grid derived from subscripting the abscissae with the same way.  For example, integration over the first axis is done by sampling on a low-rank grid

	f[None, :]


Sampling
---

The method `sampled` has an inverse, `setsamples`.

	timeslice = Field(q, Grid.from_axes([t])*R)
	results.setsamples(timeslice)

These bear a similar relation to `__getitems__` and `__setitems__`.  


Plotting
---

Fields of rank 1 and 2 have some methods to assist plotting.  These methods label the axes with the grid coordinates.

A rank 1 field is plotted against its axis as follows

	exp(x).plot()
	
If `x` has multiple components, these are plotted as separate lines.

A rank 2 field can be plotted in greyscale with

	(x*y).positive()
	(x*y).negative(interpolation='nearest')
	
and variations on this.  Positive and negative have the same meaning as in photography.


Spectral methods
---

Fourier transforms are complicated, and I haven't worked out all the details yet.

The Fourier transform of a field is returned by the method fft

	delta = ones.fft()

This is a SpectralField, wheras ones is a SampledField; these are both subtypes of ndarray.  The inverse is of course

	allclose(ones, delta.ifft())

The wavenumbers for which the elements of delta are coefficients are found by

	delta.k()

So we can take a gradient as

	zero = (1j*delta.k()*ones.fft()).ifft()

which is equivalent to

	zero = ones.D()

A SpectralField records the bounds of the field from which it was tranformed.  These, along with the Grid of wavenumbers, allow ifft to reconstruct the original grid.

How should even-sized grids be handled?  Should the reciprocal grid always be odd, with the last term possibly split between +f and -f, so that real fields are interpolated with real functions?


Efficiency goals
---

The intent is that `Grid` operations should be done once, when problems are set up, whereas `Field` operations might be called in inner loops.  Therefore, only `Field` is optimised for speed: `Grid` is optimised for maintainability and numerical correctness.  The exception is trivial grid operations, where the grids involved represent the same points.  These should be shortcut, so that users can reason in terms of equality and ignore identity.


History and motivation
---

This library was inspired by XMDS by Greg Collecutt.  With a decade of hindsight, I think it is clear that part of this system is very helpful to users, but other parts aren't flexible enough for general use, and are too complicated for general users to extend.  I've tried to copy the former and omit the latter.  Once this library is imported into Python, I hope that the rest will be fairly straightforward to implement.

The names of the methods `w` and `W` are intended to suggest a coordinate or a frequency, but be ambigous between the two.


Things to change
---

It would be nice to use the fielab idea, where subscripts in parentheses refer to coordinates and those in brackets refer to indices.  In Python, however, we can't use slice syntax in function calls, so that wouldn't be practical.  Users should refer to coordinates far more often than they refer to indices, so subscripts should be grid coordinates.  There might be trouble ensuring that a grid actually contains points, but `[w-0.5*h, w+0.5*h]` will do that.

The syntax for delta grids and low-rank grids can be as follows: `R[:, None]` is a rank 1 grid, but `R[:, pi]` is a delta grid.

The bounds of a field extend by half a grid step past each point, forming a box of size shape*h.  This is the only meaning of the grid step of a delta field.  Extrapolation is allowed inside the bounds.  A delta field is extrapolated as a constant, so that it can be resampled on a delta grid that is numerically close to its sampling plane.

Idea: Field is an abstract superclass, with subclasses for different representations, such as a sampled field and its fourier coefficients.  All of these record the bounds of the field, which, along with the grid of wavenumbers, allow the sampling grid to be reconstructed.


Limitations of Python
---

At some point, I might extend this library to a language.  Among other things, that would allow the following limitations of Python to be overcome:

Only brackets can take slice notation, so we have somewhat ugly methods sampled() and setsamples().  The language would use parentheses for these, where the index is a grid or the assigned value a field, and would use brackets for indices.  This is consistent with function calls: exp(x) can be read either as "the values of the function exp on the grid x" or "the field exp interpolated on the grid x", but these are the same thing.  It isn't clear that indirect interpolation, f(x) where f and x are both fields, is generally useful.
