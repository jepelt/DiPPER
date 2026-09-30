# Is CmdStan available, so that DiPPER models can actually be fitted?
#
# Set DIPPER_NO_CMDSTAN=true to force the CmdStan-free code paths.
has_cmdstan <- function() {
    !identical(Sys.getenv("DIPPER_NO_CMDSTAN"), "true") &&
        instantiate::stan_cmdstan_exists()
}


# Are the DiPPER Stan models available as compiled executables?
has_dipper_models <- function() {

    if (!has_cmdstan()) {
        return(FALSE)
    }

    mod <- tryCatch(
        instantiate::stan_package_model(name = "dipper_dp_asym",
                                        package = "DiPPER"),
        error = function(e) NULL
    )
    if (is.null(mod)) {
        return(FALSE)
    }

    exe <- tryCatch(mod$exe_file(), error = function(e) character(0))
    length(exe) > 0 && nzchar(exe[1]) && file.exists(exe[1])
}
