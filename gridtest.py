from field import *

x, y, z, s = [Grid.from_axes(q) for q in
	[arange(2), 0.3*arange(3), pi+arange(5), [0, 17]]]
N = Grid(shape=(), h=[], p=[], o=[1,2,3], U=empty((3,0)))
R = x*y*z
T = s*R
f = SampledField(R.blank(), R);  f[:,:,:] = 1

# check assignment worked
assert f[0,0,0] == 1

# check index and coordinate functions

e = array([[1,1,1],[1,3,2]]).T
assert all(type(x) is Field for x in [f.i(), f.r(), f.R(), f.rr()])
assert all(type(x) is ndarray for x in [R.i(w=e), R.w(i=e), R.W(i=3), R.ww(w=e)])
assert allclose(R.w(), R.w(i=indices(R.shape)))
assert allclose(R.w(W=e), e)

assert allclose(e, R.i(W=R.W(i=e)))
assert allclose(e, R.W(i=R.i(W=e)))
# fill the round-trip tests out

# test rotation
U = array([[cos(pi/6), -sin(pi/6), 0], [sin(pi/6), cos(pi/6), 0], [0, 0, 1]])
origin = R.W(w=(0,0,0)).reshape((3,1,1,1))
assert allclose(R.rotated(U).W() - origin, tensordot(U, R.W() - origin, 1))

# test subgrids
Rs = R[0:2, 1:2, 1:3]
assert allclose(Rs.w(), array(R.w())[:, 0:2, 1:2, 1:3])
assert allclose(Rs.W(), array(R.W())[:, 0:2, 1:2, 1:3])
assert Rs.close(R[0:2, 1:2, 1:-2])

# test low rank
assert R.close(R[:,:,0]*R[0,0,:])
assert R.close(R[:,0,0]*R[0,:,0]*R[0,0,:])
# assert y.close(R[pi,:,exp(1)])

# test subfields
assert (array(f[:,0,0]) == array(f)[:,0,0]).all()

# FFT tests
assert allclose(2*pi*abs(fftfreq(x.shape[0], x.h[0])).max(),
	x.reciprocal().w().max())
assert allclose(2*pi*abs(fftfreq(y.shape[0], y.h[0])).max(),
	y.reciprocal().w().max())
assert allclose(2*pi*abs(fftfreq(z.shape[0], z.h[0])).max(),
	z.reciprocal().w().max())
assert allclose(R.rotated(U).reciprocal().w(), R.reciprocal().rotated(U).w())
assert R.reciprocal().close(x.reciprocal()*y.reciprocal()*z.reciprocal())

g = SampledField(ones(x.shape), x)
G = g.fft()

# Rayleigh's theorem
assert allclose((g**2).S(), (G**2).S())

# support
assert f.support().close(f.abscissae)

# possible bugs
# interpolation on low-rank grids that are close but unequal