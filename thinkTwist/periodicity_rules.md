# Periodicity Rules For Twisted Quotient Lattice

Boundary conditions:
- (i + Nx, j) ~ (i, j)
- (i, j + Ny) ~ (i + Ns, j)

Direction vectors:
- x = (1, 0)
- y = (0, 1)
- y_minus_x = (-1, 1)

Closed-form periods:
- p_x = Nx
- p_y = Nx*Ny / gcd(Nx, Ns)
- p_y_minus_x = Nx*Ny / gcd(Nx, Ny + Ns)

## Rule For Triangle Side Lengths From Periodicities

Let the local triangle side lengths be:
- lx = |x|
- ly = |y|
- lr = |y-x|

If you want target half-orbit physical scales (Dx, Dy, Dr), the matching rule is:
- (p_x/2)*lx = Dx
- (p_y/2)*ly = Dy
- (p_y_minus_x/2)*lr = Dr

So:
- lx = 2*Dx / p_x
- ly = 2*Dy / p_y
- lr = 2*Dr / p_y_minus_x

For square-continuum targets Dx:Dy:Dr = 1:1:sqrt(2), lengths are fixed up to one overall scale s:
- lx = s
- ly = (p_x/p_y)*s
- lr = sqrt(2)*(p_x/p_y_minus_x)*s

Feasibility condition (a real triangle must exist):
- |lx - ly| <= lr <= lx + ly

Equivalent scale-free condition for the square-target case:
- |1 - p_x/p_y| <= sqrt(2)*p_x/p_y_minus_x <= 1 + p_x/p_y

If this fails, no Euclidean triangle can exactly realize those target directional scales.

Notes:
- Ns is effectively modulo Nx.
- Each formula returns an integer period in [1, Nx*Ny].
- A numerical verifier script is provided in thinkTwist/verify_direction_periodicity_rules.py.
