# Physics Gates

These gates are mandatory before final delivery.

## Global
- [ ] Potential and energy are distinct.
- [ ] SI conversion is explicit.
- [ ] `EC - EV = Eg`.
- [ ] No NaN/Inf in valid outputs.
- [ ] Deterministic common energy reference.

## PN
- [ ] `Vbi = VT ln(NA ND / ni^2)`.
- [ ] `VD = Vbi - VA`.
- [ ] `W`, `xp`, `xn` correct.
- [ ] `NA*xp = ND*xn`.
- [ ] `Emax` identities agree.
- [ ] Potential drop equals `VD`.
- [ ] `VA=0` gives flat EF.
- [ ] `VA!=0` does not give a global equilibrium EF.
- [ ] `EFn_n - EFp_p = VA` in eV.
- [ ] No fabricated depletion-region quasi-Fermi profile.
- [ ] `VD<=0` is rejected.

## BJT
- [ ] `VBE = VB - VE`.
- [ ] `VBC = VB - VC`.
- [ ] Both junction barriers validated.
- [ ] Four operating regions classified.
- [ ] Near-zero bias tolerance exists.
- [ ] Base depletion overlap detected.
- [ ] No fabricated spatial quasi-Fermi transport.

## MOS
- [ ] x=0 is Si/SiO2 interface.
- [ ] x>0 is semiconductor bulk direction.
- [ ] bulk potential is zero.
- [ ] `psi_s > 0` bends bands downward.
- [ ] `VFB = PhiMS - Qox/Cox`.
- [ ] Level-1 VT equation validated.
- [ ] accumulation/depletion/inversion classification validated.
- [ ] Level-2 solver is bounded.
- [ ] solver residual/convergence are exposed.
- [ ] silicon bands are not drawn inside oxide.
