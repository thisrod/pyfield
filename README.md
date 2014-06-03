Field library
===

by Rodney E. S. Polkinghorne

A field is a scalar, vector or tensor quantity that depends on position and time.  This Python library exports a `SampledField` type, being an `ndarray` containing samples of an array-valued quantity, that also records the points at which the samples were taken.  These points form a rectangular grid in R<sup>n</sup>.  Such an array can integrate, differentiate and Fourier transform itself, generate samples of white noise at its points, and so on.  Other representations, such as `SpectralField`, represent the same data in different ways.  Each of these can be regarded as an expansion over a different basis.

The library aims to remove the accidental complexity from computing with fields, in order to spare scientific programers from bookeeping, to prevent large classes of bugs from occuring at all, and to allow the remaining code to address physical problems, and the remaining bugs to be removed by physical nous.

The current version of the library assumes that the quantity being sampled is an array, whose components do not vary as the coordinates change underneath them.  Future versions might understand how the components of vectors and tensors should transform in different coordinate systems.  

The library code is a draft.  Many user actions that ought to generate a meaningful exception instead cause an assertion to fail, and sometimes that assertion has a comment "n.y.i.", meaning not yet implemented.  However, this document should describe exactly what the user is allowed to do, and what the library should do in response, and the code should either do that or fail deliberately.  Please report exceptions to that.

Concept: the limited bandwidth interpolant of f_i=sin(kn), -N<n<N, f_i=0 otherwise.  Do these span the same space as the sincs?  Yes, because the interpolants can be written as linear combinations of the sincs.

Spectral derivatives might go nuts at the edges of a grid, due to Gibbs effect.  But this isn't an issue if the field goes smoothly to zero at the edges: it can be extended to a periodic field that is everywhere differentiable.

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

	S[1,2].W()

A `Grid` constructed like this is parallel to the axes of R<sup>n</sup>, but more general grids can be derived by rotation, specified by an orthogonal matrix.

	U = array([[cos(1), sin(1)], [-sin(1), cos(1)]])
	T = S.rotated(U)

As well as rotating grids, we can translate them.

	S.translated((0, -pi))

is a grid similar to `S`, but starting from the origin of the plane.

TODO Implement the full set of MetaFont transforms.

TODO Grids in the Argand plane, complex h, shape is a Gaussian integer

As well as the common coordinates returned by `W()`, each grid has a system of grid coordinates.  These can be obtained by casting the grid to an ndarray.

	array(S)
	array(S.rotated(U))
	array(S.translated((0, -pi)))

Grid coordinates are preserved by translations and rotations.  The method `shifted` translates the origin of the grid coordinates while leaving the common coordinates fixed

	S.shifted((1,1)).w()

The method `i` calculates indices from coordinates

	S.i()
	S.i(w=array(S))
	S.i(W=S.W())

These methods return a `Field` if called with no arguments or with a `Field`, and an `ndarray` if passed an `ndarray`.

A grid has notional bounds, that extend half a grid step past the extreme points at each edge.

	S.bounds()

This has shape 2*one.rank().  This is most useful in telling the inverse fourier transform where to put the reconstructed field.


Fields
---

A `SampledField` is constructed from a `Grid` and a set of data sampled on it, as follows

	one = SampledField(ones(S.shape), S)

This is an `ndarray`, so we can do things like

	allclose(one + exp(1j*pi*one), 0*one)

In fact, the array `S.W()`, which can also be expressed as `one.R()`, is a field on `S`, so

	r0 = (one.R()*one)[None]/one[None]

is a complicated way to find the centre of mass of a rectangle, in common coordinates.  Indexing a field with `None` computes its integral.

When a field has more than one component, we need a convention about the order that the axes go in when the samples are represented as an array.  The convention is

other components*axis components*sample points.

So, for example, samples of the gradient of the wavefunction of a spin 1/2 particle on a 4*5*6 grid would be stored in an array of shape 2*3*4*5*6.  The 2 is the spin components, and the 3 is the x, y and z components of the gradient vector.


Sampling
---

Fields can be interpolated on other grids

	a = one.sampled(T)

In general, with special relevance to sampling, a field is treated as the limited bandwidth interpolant of the samples.  This is equivalent to adding up sinc functions at every grid point.  Extrapolated values take the values of the sincs outside the grid, which asymptote to zero with at least the reciprocal of the distance from the grid.

An important special case is a grid with shape 1 along some dimension.  Such a grid has a step, and can be constructed with `Grid.delta(x, h)`, or by indexing with `R[0:1,:,:]`.  On the degenerate axis, it varies as sinc(2&pi;x/h), with the first zeros at &pm;h.  This is consistent with the general rule that fields are limited bandwidth interpolants of their samples.  It means that taking slices with spacing h, on grids of width h, then adding them back together, will reconstruct a bandwidth limited interpolant of the original field.

TODO implement `along` to make deriving an oblique delta grid less bodgy.

A grid with one point and a step of zero is treated as a sinc in the limit of zero width; this is the default of `Grid.delta(x)` and `Grid.from_axes([x])`.  This will give zero when sampled on any other grid; however, its values can still be plotted or extracted by indexing.

The corresponding assignment operation is `setsamples`

	a[:,:] = 2
	one.setsamples(a)
	
The algorithm to do this is as follows.  The target of the assignment, the field `one` is sampled on `T`, the grid of the assigned value `a`.  This gives a set of samples on `T`, which are subtracted from the samples of `a`.  This difference is interpolated back on `R`, the grid of `one`, and the result added to the samples of `one`.

In the case that the grids coincide, this reduces to assignment.  There are some other desirable properties we should investigate.  Does setting samples on parts of a field commute, when the bounds of those parts don't overlap?  Also, an array of results can be constructed by adding slices at integration timesteps to an accumulator field with a larger step.  Does assigning the slices to the accumulator do the same thing, regardless of its initial value?


Degenerate grids
---

Grids are still more general: the dimension may exceed the rank.  Such low-rank grids can be constructed by subscription, e.g., S[2, :].  Fields over such a grid are constant over the missing axes.  Subscripting a field in this way constructs a field whose values are a slice of the subscripted field at the subscripts.

A grid with rank less than its dimension can also be constructed by supplying `None` as a subscript.  We already saw this when `x` and `y` were constructed from `S`.  For example, it might be a function of space in a dynamical simulation.  Another example is a field whose value is one of the grid coordinates, which is constant over the other coordinates.

Sampling on a low-rank grid integrates over missing axes.

Subscripting a field with `None` causes it to be sampled on the grid derived from subscripting the abscissae with the same way.  For example, integration over the first axis is done by

	f[None, :]


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

A coloured phase plot of a complex field, as in ?, is made by

	exp(z).argand()


Spectral methods
---

Fourier transforms are complicated, and I haven't worked out all the details yet.

The Fourier transform of a field is returned by the method fft

	delta = ones.fft()

This is a SpectralField, wheras ones is a SampledField; these are both subtypes of ndarray.  The inverse is of course

	allclose(ones, delta.ifft())

An array of the wavenumbers for which the elements of delta are coefficients are found by

	delta.k()

So we can take a gradient as

	zero = (1j*delta.k()*ones.fft()).ifft()

which is equivalent to

	zero = ones.D()

A SpectralField records the orientation, bounds, and grid origin of the field from which it was tranformed.  These, along with the Grid of wavenumbers, allow ifft to reconstruct the original grid.

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
