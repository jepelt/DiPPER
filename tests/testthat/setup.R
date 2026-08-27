# Is CmdStan available, so that DiPPER models can actually be fitted?
#
# Set DIPPER_NO_CMDSTAN=true to force the CmdStan-free code paths.
has_cmdstan <- function() {
    !identical(Sys.getenv("DIPPER_NO_CMDSTAN"), "true") &&
        instantiate::stan_cmdstan_exists()
}
