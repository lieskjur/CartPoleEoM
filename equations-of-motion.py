#
#                             //
#                   \        //
#                   /\      //
#                  /  \    //
#                 /    \  //
#               r_p     \//
#               /       >< m_p
#              / |-th->//
#             /  |    //
#           \/   |   //
#            \   |  //
#             \  | //
#   |   |------\-|//------|
#   |   |       \//       |
#   g   |  m_c  ( )       |
#   |   |                 |
#   V   |=================|
# _________OOO___|___OOO_________
#                |
# //|-----s----->|

from sympy import symbols, Function, Matrix, sin, cos, diff, latex, simplify
from lagrangian2eom import massMatrix, velocityTerms, potentialTerms

t = symbols("t")
s = Function("s")(t)
theta = Function("theta")(t)

g, mass_c, l, mass_p, moi_p = symbols("g m_c l m_p I_p")

# generalized coordinates
q = Matrix([s, theta])

# Kinetic energy
dxp = diff(s, t) + l * diff(theta, t) * cos(theta)
dyp = -l * diff(theta, t) * sin(theta)

kin = 0.5 * (mass_c * diff(s, t) ** 2 + mass_p * (dxp**2 + dyp**2))

# Potential energy
pot = g * mass_p * l * cos(theta)

# Quantities
M = massMatrix(kin, q, t)
c = velocityTerms(kin, q, t)
tau_p = potentialTerms(pot, q)

# Printout
print("T &=", latex(simplify(kin)), "\\\\")
print("V &=", latex(simplify(pot)), "\\\\")
print("M &=", latex(simplify(M)), "\\\\")
print("c &=", latex(simplify(c)), "\\\\")
print("\\tau_p &=", latex(simplify(tau_p)))
