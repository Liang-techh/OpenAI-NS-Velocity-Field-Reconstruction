# Equation (5.6) regular Omega source provenance

Source: OpenAI, *Finite Time Blowup for Navier–Stokes*, Section 5, Eqs. (5.2), (5.6), and the paragraph immediately following (5.6). The repository remains pinned to official Lean commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd` for formal cross-checks.

The paper defines, at coefficient order `k`,

`Omega_k = T_{0,k} V_k + sum_{i+j=k} [ V_i (d_X V_j - V_j/(2X)) + U_i Z_{0,j} V_j ] - 2 X d_XX V_k - Z^[2]_{0,k-1} V_{k-1}`,

with `T_{a,n}=T_{a+lambda_n}`, `Z_{a,n}=Z_{a+lambda_n}`, `Z^[2]_{a,n}=Z_{a+lambda_n-D} Z_{a+lambda_n}`, `lambda_n=2nh`, and negative coefficient indices interpreted as zero.

Equation (5.2) gives `V_j = X v_j` with smooth `v_j`. Dividing Eq. (5.6) by `X` therefore has a regular axis extension. `background_recurrence.py` implements that quotient directly from the second jet of `v_j`; it does not numerically divide an evaluated `Omega_k` by `X`.

The regularized terms used by the implementation are:

- `T_b(Xv)/X = L^-1 [ (1-b)v + D eta v_eta + X v_X ]`.
- `Z_b(Xv)/X = L^-1 [ 2 eta (b-1)v + d v_eta - 2 eta X v_X ]`.
- `V_i (V_j' - V_j/(2X))/X = v_i (v_j/2 + X v_j')`.
- `-2X V_k''/X = -2(2 v_k' + X v_k'')`.
- The shifted axial term is evaluated by writing `Z_b(Xv)=Xw` and applying the same regular `Z` formula to `w`, analytically differentiating `w`; no finite-difference derivative is used in production.

Tests independently compare these regular expressions with the unreduced Eqs. (4.2)/(5.6) away from the axis. The shifted double-`Z` test deliberately uses finite differences only in the *test oracle* to avoid duplicating the production derivative algebra. Tests also evaluate the regular quotient at `X=0` and compare with its small-positive-`X` limit.

Truth boundary: this is a paper-derived pointwise recurrence row and a genuine constructor for the regular source `Omega_k/X` once the lower-order coefficient jets are known. It does **not** solve Eqs. (5.3)-(5.4), does not yet produce the coefficient profiles from the leading field, and does not upgrade Stage 2 beyond `formal-structure`.
