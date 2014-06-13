Field library
===

by Rodney E. S. Polkinghorne

A field is a scalar, vector or tensor quantity that depends on position and time.  This Python library exports a `SampledField` type, being an `ndarray` containing samples of an array-valued quantity, that also records the points at which the samples were taken.  These points form a rectangular grid in R<sup>n</sup>.  Such an array can integrate, differentiate and Fourier transform itself, generate samples of white noise at its points, and so on.  Other representations, such as `SpectralField`, represent the same data in different ways.

The classes can be regarded as expansions of the same field over different bases.  A `SampledField` expands over sinc functions, a `SpectralField` expands over complex exponentials, and so on.  These bases are indexed by a `Grid`.

The library aims to remove the accidental complexity from computing with fields, in order to spare scientific programers from bookeeping, to prevent large classes of bugs from occuring at all, and to allow the remaining code to address physical problems, and the remaining bugs to be removed by physical nous.

The current version of the library assumes that the quantity being sampled is an array, whose components do not vary as the coordinates change underneath them.  Future versions might understand how the components of vectors and tensors should transform in different coordinate systems.  

The library code is a draft.  Many user actions that ought to generate a meaningful exception instead cause an assertion to fail, and sometimes that assertion has a comment "n.y.i.", meaning not yet implemented.  However, this document should describe exactly what the user is allowed to do, and what the library should do in response, and the code should either do that or fail deliberately.  Please report exceptions to that.

Concept: the limited bandwidth interpolant of fi=sin(kn), [-N,N], fi=0 otherwise.  Do these span the same space as the sincs?  Yes, because the interpolants can be written as linear combinations of the sincs.


Grids and fields
---

The simplest representation of a field is an array of the values it takes at points on a rectangular grid.  The library implements this with the class `SamplingGrid`.  This represents a grid of points in R<sup>n</sup>, and a system of coordinates aligned with the grid axes.  These need not be aligned with the usual axes of R<sup>n</sup>, and the origin need not be a point of the grid.

The superclass constructor `Grid` returns a `SamplingGrid`.  The usual way to construct one is from arrays of the points on its axes, as

	S = Grid(0.1*arange(3), pi+0.5*arange(4))

The `*` operator acts as a cartesian product on `Grid`, so the same grid can be constructed with

	x, y = Grid(0.1*arange(10)), Grid(pi+0.5*arange(20))
	S = x*y
	
However, the points of `S` now lie in the plane, while `x` and `y` comprise points on the line.  This prevents the computer from associating `x` with the first axis of `S` and `y` with the second.  For that reason, it is better to construct the grid the first way, then project out `x` and `y` with

	x, y = S.axes
	
Now `x` and `y` lie in the plane, so `x*y` and `y*x` both reconstruct `S`.  

A `Field` is an array of samples on a grid.  It can be constructed as

	g = sin(Field(x))
	
Here, `Field(x)` constructs a field whose value is the grid coordinate x.  The coordinate grid `x` is stored as `g.abscissae`.  It would be nice to say, `sin(x)`, but implementing that in Python is too complicated.  A `Field` is an `ndarray`, and can be constructed from a `Grid` and an array with the same shape

	g = Field(x, zeros(x.shape))

It would also be nice to say, `Field(x) + Field(y)`, but the terms are one-dimensional arrays with different shapes, and Numpy is jealous of its broadcasting rules.  Instead, we have to resample on the grid `S`, that spans `x` and `y`

	g = Field(x)[S] + Field(y)[S]
	
Subscripting a `Field` with a `Grid` causes the field to be interpolated on the points of the grid.

The construction

	q = Field(S)

returns a `2*S.shape` array.  The plane `q[0,:,:]` is the same as `Field(x)[S]`: it has the shape of the grid `S`, and contains the x coordinates of its points.  Likewise, `q[0,:,:]` has the `y` coodinates.

The grids seen so far are the simplest kind possible.  They all start at the origin, and are aligned with the axes of the common coordinate space R<sup>n</sup>.  More general grids can be constructed by geometric transformations, that maintain a rectangular grid, and leave each point at the same grid coordinates, but move them together in the common space.

	U = array([[cos(1), sin(1)], [-sin(1), cos(1)]])
	T = S.rotated(U).translated([5, pi])

The only non obvious point about translations is that the argument is specified in common coordinates, not in grid coordinates.  Rotations take a unitary matrix, that transforms a unit vector along an axis of the old grid to one along the corresponding axis of the rotated grid.


Grid and common coordinates
---

At this point, having constructed the grid `T`, we can construct such things as

	Field(T)[S]

and these require us to think more closely about the coordinates of each point in the grid `T` and their relation to coordinates in the common space.  The grid coordinates of the point `T[1,0]` are `[0.1, pi]`, the distances along the grid axes from the grid origin.  Its common coordinates are some rotation and translation of this vector.  The expression `Field(T)[S]` will interpolate coordinates in the grid `T` for each point of the grid `S`; the grid `S` will not actually have points at these coordinates, except by coincidence.

It is possible to shift the origin of the grid coordinates of `T` away from the point `T[0,0]` with index zero.  This is done with 

	T.shifted([-1, 2.5])

This has its origin such that the grid coordinates of the index origin, `T.shifted(a)[0,0]`, are the vector `a`.


Interpolation
---

It is obvious what coordinates should be returned by `Field(T)[S]`, but to do this more generally, for example `sin(Field(x))[x.translated(0.5)]`, we need some methods of interpolation and extrapolation: the sines of the coordinates, sampled at the points of `x`, must be extended to some notional function.  Different types of grid do this differently, but the ones considered so far use band limited interpolation.  A sinc function is placed at each grid point, with widths along each axis such that it has zeros at the other points.  These sinc functions are combined with coefficients equal to the sample values.  This gives a sensible result when the sampled function is smooth on the scale of the grid step, and decays almost to zero at the edges.  Otherwise, a band-limited interpolant will oscillate wildly where the sampled function changes sharply.

This form of interpolation is used consistently, whenever a `Field` needs to be treated as a continuous function.  For example, if we say

	f = Field(x)[S] + Field(y)[S]
	g = f[x]

the library interprets `g` as the function of `x` obtained by integrating `f` along the `y` axis.  It obtains this by integrating the sinc interpolant, and sampling the result at the points of `x`.


Projection grids
---

This brings up a point that we've been glossing over so far.  The grids `S`, `x` and `y` all lie in the plane, as indicated by their `dimension` attribute

	assert S.dimension == 2
	assert x.dimension == 2

However, the points of `S` are arranged along 2 axes, while those of `x` and `y` are collinear.  This can be reported as follows

	assert S.rank == 2
	assert x.rank == 1

We will refer to grids such as `x` as projection grids, because a field subscripted with one is integrated over the missing axes.  Conversely, a `Field` such as `sin(Field(x))` is interpreted as constant over the `y` axis: a special case of this rule is `Field(x)[S]`.  This allows outer products to be taken easily, as in

	sin(Field(x))[S]*exp(-Field(y)**2)[S]

A "point" of a projection grid is better thought of as a line, plane, or general slice along the missing axes.  Field samples are constant over this plane, and sampling integrals are taken over it.  One consequence of this is that the grid can be translated parallel to these slices with no effect, and the component of the grid origin that is orthogonal to them is the only meaningful part.

A geometric method applies particularly to projection grids.  

	y.along(S)

returns a grid with points at the same grid coordinates as those of `y`, but with the origin of those coordinates the same as that of `S`, and oriented along the first axis of `S`.  In general, the axes of the result are oriented along the first `y.rank` axes of `S`.

If two projection grids have the same dimension, and their axes are all numerically orthogonal, their cartesian product lies in the same common space as them.  So, for example, the product of `S.axes` is always `S`.  Otherwise, the product lies in the cartesian product of the common spaces, and the axes of each grid keep their orientation in the subspace corresponding to the original common space.  The axes of the product are ordered in the obvious way, from left to right through the factors.


Section grids
---

It is possible to have a grid with only one point along a certain axis, for example `S[0:1,:]`.  This is distinct from a projection grid: wheras `y` is conceptually a set of lines parallel to the `x` axis, this is a particular set of points on these lines.  This grid has a notional step along its first axis, the same as the step of `S`.  This is nonzero, so samples on the grid are extrapolated to a sinc function along the axis; this follows the general rule, except that there is only one sinc function.

Suppose `R` is a grid covering the same space as `S`, but sampled less densely, and `f` is a field over `S`. Then f is the sum of the fields `f[n:n+1,:]`.  Also, resampling distributes over addition: if the fields `f[n:n+1,:][R]` are added up, the sum is equal to `f[R]`.  This provides a convenient way to interpolate the solution to a differential equation on a grid coarser than the timestep.

If the step is zero, the sinc function has zero width.  We will call this a section grid.  A continuous function extrapolated from it is zero almost everywhere, and the integral of such a function is zero.  The only case in which resampling this field gives a nonzero result is when the sampling grid is another section grid, that is spanned by the original one.

The result of `through` is a section grid.   It can also be constructed with `Grid.delta(h)`, which can be aligned to an existing grid with `along`.

In general, `through` returns a grid that has been extended to full rank, i.e. rank equal to dimension, for which the specified point lies in the span.  The span of a grid is the set of points that can be reached from its origin by travelling in the directions of its axes and, in the case of a projection grid, along its slices.  For anything but a section grid, this is the entire space.  The inclusion of one span in another can be tested with the predicate `R.spans(S)`.


Plotting
---

Fields of rank 1 and 2 have some methods to assist plotting.  These methods label the axes with the grid coordinates.  These methods can also be used on section grids, where all but 1 or 2 axes have zero step.

A rank 1 field is plotted against its axis as follows

	exp(x).plot()
	
If `x` has multiple components, these are plotted as separate lines.

A rank 2 field can be plotted in greyscale with

	(x*y).plot_positive()
	(x*y).plot_negative(interpolation='nearest')
	
and variations on this.  Positive and negative have the same meaning as in photography.

A coloured phase plot of a complex field, as in whatsisname, is made by

	exp(z).plot_argand()


Fourier transforms
---

A `Field` started off as an array of samples.  Then we defined a standard way to interpolate these, so that it defined a continuous function.  This function is a linear combination of a set of sinc functions defined by the grid, with coefficients equal to the samples in the field.  We mentioned that there are other kinds of grids, that construct continuous functions in different ways.  All of these ways are linear combinations with the `Field` array as the coefficients: what changes is the basis functions defined by the `Grid`.

Besides the `SamplingGrid` that we've been using until now, the next most common case is the `ReciprocalGrid`.  The basis functions in this case are complex exponentials.  More precisely, the basis functions of `S.reciprocal()` are the extrapolated functions that would result from sampling complex exponenials on the points of `S`, with wavenumbers such that they were periodic across it.  Therefore, a field over `S` represents the same function as one over `S.reciprocal()` if the elements of one are the discrete Fourier transform of the samples of the other.  If `f` is a field over `x`, then `f[x.reciprocal()]` is its Fourier transform.

The coordinates of a reciprocal grid are the wavenumbers of the complex exponentials.  So we take a spectral derivative with

	k = x.reciprocal()
	ftick = (-1j*Field(k)*f[k])[x]

This can be done with the method

	ftick = f.D()
	
In general, this returns a vector field that is the gradient of `f`.  If `f` is already a vector field, `f.D()` returns a field of its Jacobian matrix.  Every type of grid has a way to compute the derivative of its interpolant, though this might not always be as efficient as taking a Fourier transform.

A `SamplingGrid` has a method `bounds()`, that returns the coordinates half a step past the end of each axis, as a 2*rank array.  The corresponding `ReciprocalGrid` remembers these; along with the wavenumbers, they should uniquely specify a box in R<sup>n</sup> and a grid in it.

There are some complications to this, such as how you handle the extreme wavenumber of a grid with an even number of points.  These will be specified later.


TODO Grids in the Argand plane, complex h, shape is a Gaussian integer

TODO specify `spans`

TODO document bounds, and how they are preserved by reciprocal grids.

TODO gaussian noise

TODO extension to Hilbert space, grids with variational parameters, derivative wrt those parameters.


Subscripting
---

This is getting complicated.  There is a different way of factoring the classes, which would simplify some things.  A `Grid` represents a basis for a space of functions of R<sup>n</sup>, and a `Field` is simply an array of coefficients of those functions.  Resampling would mean expanding over the new basis.  (Basis isn't the right word, because it isn't complete.)  So, for example, a `SamplingGrid` represents a basis of sinc functions, and a `ReciprocalGrid` represents the bounded complex exponentials that span the same space.  In this picture, it makes total sense that a Fourier grid remembers its bounds!  This pushes even more of the machinery into the grid classes.

This also allows some things to be simpler and more general: for instance, a Fourier transform is a resampling on a `ReciprocalGrid`, and we could construct a cartesian product of a `ReciprocalGrid` on some axes and `SamplingGrid` on the others, perhaps to represent wavenumbers as a function of time.  The model even extends to wavelet and Gabor frames, coherent states and Fock states!  But what about finite difference formulae, and other things that aren't easily expressed in terms of bases?

In this model, there are two ways to subscript a field.  We might want to extract certain coefficients, which we could do with an index, or an array of indices.  The semantics of this are the same as for ndarrays.  We might also want to expand over a different basis, which can be done by subscripting with the new grid.

Subscripting a field with `None` causes it to be sampled on the grid derived from subscripting the abscissae with the same way.  For example, integration over the first axis is done by

	f[None, :]

Grids can be subscripted in a variety of ways.  Section means a grid with one point and step 0, while projection means a low-rank grid.

<table>
<tr><th>subscript type</th><th>effect</th></tr>
<tr><td>integer</td><td>section by index</td></tr>
<tr><td>integer slice</td><td>subgrid by index</td></tr>
<tr><td>singleton "</td><td>sinc subgrid</td></tr>
<tr><td>float</td><td>section by grid coordinates</td></tr>
<tr><td>float slice</td><td>subgrid by grid coordinates, same step, symmetric in interval defined by slice</td></tr>
<tr><td>ends nearly equal</td><td>sinc subgrid at coordinate</td></tr>
<tr><td>Grid</td><td>God knows</td></tr>
<tr><td>1D array</td><td>grid subscripted with array elements</td></tr>
<tr><td>ndarray</td><td>array of grids, with lists of subscripts along last axis of array</td></tr>
</table>

It is an error to subscript with a mixture of ints and floats.  Note that fully subscripting a grid gives a zero-dimensional section.  Subscripting a field with this gives a field with a value at one point, which, under the semantics of ndarras, is nearly equivalent to the value.

An important special case is a grid with shape 1 along some dimension.  Such a grid has a step, and can be constructed with `Grid.delta(x, h)`, or by indexing with `R[0:1,:,:]`.  On the degenerate axis, it varies as sinc(2&pi;x/h), with the first zeros at &pm;h.  This is consistent with the general rule that fields are limited bandwidth interpolants of their samples.  It means that taking slices with spacing h, on grids of width h, then adding them back together, will reconstruct a bandwidth limited interpolant of the original field.

A _section grid_ with one point and step zero is treated as a sinc of zero width.  This is the default returned by `Grid.delta(x)` and `Grid.from_axes([x])`.  Sampling a section on a grid whose points are close to its span interpolates in the span.  Otherwise, sampling from a section gives zero: this is always the case when a section is sampled on a grid with a higher dimensional span.

The corresponding assignment operation is `setsamples`

	a[:,:] = 2
	one.setsamples(a)
	
The algorithm to do this is as follows.  The target of the assignment, the field `one` is sampled on `T`, the grid of the assigned value `a`.  This gives a set of samples on `T`, which are subtracted from the samples of `a`.  This difference is interpolated back on `R`, the grid of `one`, and the result added to the samples of `one`.

In the case that the grids coincide, this reduces to assignment.  There are some other desirable properties we should investigate.  Does setting samples on parts of a field commute, when the bounds of those parts don't overlap?  Also, an array of results can be constructed by adding slices at integration timesteps to an accumulator field with a larger step.  Does assigning the slices to the accumulator do the same thing, regardless of its initial value?


History and motivation
---

This library was inspired by XMDS by Greg Collecutt.  With a decade of hindsight, I think it is clear that part of this system is very helpful to users, but other parts aren't flexible enough for general use, and are too complicated for general users to extend.  I've tried to copy the former and omit the latter.  Once this library is imported into Python, I hope that the rest will be fairly straightforward to implement.

The names of the methods `w` and `W` are intended to suggest a coordinate or a frequency, but be ambigous between the two.


Limitations imposed by Numpy
---

This is a set of notes on things that are currently too complicated to be worth implementing, that would become practical if Numpy were extended in certain ways.

There is no easy way to override the way `ufunc` broadcasts over arrays.  Since `Grid` knows the geometry, in principle it could figure out which axes are aligned and which are orthogonal.  In practice, it isn't too much work to resample on a grid that spans everything.  Instead of `a + b`, we can say `a[R] + b[R]`.  This statisfies Python's taste for explicitness.

It would be nice for `sin(x)` to be a `Field` when `x` is a `Grid`.  The `__array__` special method comes tantalisingly close, but there's no way to preserve the `Field` subclass.  This just adds some line noise, with sin(Field(x)).