# Nonzero-Neumann wall-data Poisson MMS

This verifies prescribed outward-normal derivative data for the scalar
electric-potential reference solve. It is not a MUG result.

Manufactured solution:
- `phi = (x - 1/2)^2 + (z - 1/2)^2 - mean(phi)`.
- `laplacian(phi) = 4`.
- outward `dphi/dn = 1` on all four rectangle walls.
- mean(phi) is fixed to zero as the Neumann gauge.

Result:
- n = 65.
- max error = 3.220e-15.
- rms error = 1.533e-15.
- solved mean(phi) = -1.245e-16.

Red-team interpretation:
- This checks scalar wall-gradient plumbing only.
- It does not yet set the wall gradient from an actual `u x B` field.
- It does not include Robin wall conductance or material interfaces.
