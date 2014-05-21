Field library
===

by Rodney E. S. Polkinghorne

A field is a scalar, vector or tensor that depends on position and time.  The `Field` exported by this Python library is an `ndarray` containing samples of a field, that also records the points in R<sup>n</sup> at which the samples were taken.  Such an array can integrate, differentiate and Fourier transform itself, generate samples of white noise at its points, and so on.  The library aims to remove the accidental complexity from computing with fields, in order to spare scientific programers from bookeeping, to prevent large classes of bugs from occuring at all, and to allow the remaining code to address physical problems, and the remaining bugs to be removed by physical nous.

The mechanics of this are done by a class `Grid`.  This is a abstract representation of a rectangular grid of points in R<sup>n</sup>, along with a system of grid coordinates.  The points must form a rectangular grid, but this does not have to be aligned with the usual axes or start at the origin.  We can represent a skew plane in space, and things like that.

At `Grid` is constructed from its axes, as

	R = Grid(0.1*arange(3), pi+0.5*arange(4))

An array of coodinates for all grid points in R<sup>n</sup> is given by

	R.R()

and coordinates for the point with a known index, say `(1, 2)`, by

	R.R(i=(1, 2))

The `*` operator on `Grid` is a cartesian product, so the same grid can be constructed with

	x, y = Grid(0.1*arange(10)), Grid(pi+0.5*arange(20))
	R = x*y

A `Grid` constructed like this is parallel to the axes of R<sup>n</sup>, but more general grids can be derived by rotation, specified by an orthogonal matrix.

	U = array([[cos(1), sin(1)], [-sin(1), cos(1)]])
	R.rotated(U)

A `Field` is constructed from a `Grid` and a set of data sampled on it, as follows

	one = Field(ones(R.shape), R)

This is an `ndarray`, so we can do things like

	allclose(one + exp(1j*pi*one), 0*one)

In fact, the array `R.R()`, which can also be expressed as `one.R()`, is a field on `R`, so

	r0 = (one.R()*one).S()/one.S()

is a complicated way to find the centre of mass of a rectangle.  The method `S` computes the integral of a `Field`.


As well as rotating grids, we can translate them.

	S = R.translated((0, -pi))

is a grid similar to `R`, but starting from the origin of the plane.  Fields can be interpolated on other grids

	a = one.sampled(S)
	b = one.sampled(R.rotated(U))

Some points in `a` and `b` have been extrapolated; these samples have the value `nan`.

As well as the common coordinates returned by `R()`, each grid has a system of grid coordinates, returned by `r()`.  These are preserved by translations and rotations.

	R.r()
	R.rotated(U).r()
	R.translated((0, -pi)).r()

However, the method `shifted` translates the origin of the grid coordinates

	R.shifted((1,1)).r()

The methods `R` and `r` can take an array of indices or of coordinates

	allclose(R.r(), R.r(i=indices(R.shape)))
	allclose(R.r(), R.r(R=R.R()))

The method `i` calculates indices from coordinates

	R.i(r=R.r())
	R.i(R=R.R())

These methods return a `Field` if called with no arguments or with a `Field`, and an `ndarray` if passed an `ndarray`.


Spectral methods
---

Fourier transforms are complicated, and I haven't worked out all the details yet.

0. Are the original and the transform one `Field`, or two?

2. How should even-sized grids be handled?  Should the reciprocal grid always be odd, with the last term possibly split between +f and -f, so that real fields are interpolated with real functions?

1. How does the FT of a field remember the original grid, so that the inverse transform can translate back to the right place?  If even grids have odd shaped transforms, how does the inverse transform remember the shape of the original field?


Low rank grids
---

We said at the start that a `Grid` need not be aligned with the usual axes, or start at the origin.  There is a further generalisation: the dimension of a `Grid` may exceed its rank.  For example, a rank 2 grid might have points lying in R<sup>3</sup>, on the plane x+y+z=1.  This is the most general case.

Sampling on a one-point grid interpolates a section.

Sampling from a one-point grid multiplies by a delta function at the point.

Sampling on a low-rank grid integrates over missing axes.

Sampling from a low-rank grid assumes constant over missing axes.

The above rules make sampling between one-point and low-rank grids an identity.


Efficiency goals
---

The intent is that `Grid` operations should be done once, when problems are set up, whereas `Field` operations might be called in inner loops.  Therefore, only `Field` is optimised for speed: `Grid` is optimised for maintainability and numerical correctness.  The exception is trivial grid operations, where the grids involved represent the same points.  These should be shortcut, so that users can reason in terms of equality and ignore identity.


History and motivation
---

This library was inspired by XMDS by Greg Collecutt.  With a decade of hindsight, I think it is clear that part of this system is very helpful to users, but other parts aren't flexible enough for general use, and are too complicated for general users to extend.  I've tried to copy the former and omit the latter.  Once this library is imported into Python, I hope that the rest will be fairly straightforward to implement.

Many details were inspired by the books of Nick L. Trefethen, the clarity of whose thoughts and programs I greatly admire.

The names of the methods `w` and `W` are intended to suggest a coordinate or a frequency, but be ambigous between the two.


Things to change
---

It would be nice to use the fielab idea, where subscripts in parentheses refer to coordinates and those in brackets refer to indices.  In Python, however, we can't use slice syntax in function calls, so that wouldn't be practical.  Users should refer to coordinates far more often than they refer to indices, so subscripts should be grid coordinates.  There might be trouble ensuring that a grid actually contains points, but `[w-0.5*h, w+0.5*h]` will do that.

The syntax for delta grids and low-rank grids can be as follows: `R[:, None]` is a rank 1 grid, but `R[:, pi]` is a delta grid.