## Most tests deliberately use tiny simulations to keep the suite fast. Test
## the user-facing finite-locus warning explicitly rather than repeating it in
## every fixture.
options(rBahadur.warn_small_m = FALSE)
