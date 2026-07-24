# Validate every successful solve and conservation record in a HICAR adjoint
# projection log.  A malformed or non-finite numeric token is a hard failure;
# awk's implicit string-to-number conversion must not turn it into zero.

function is_number(value) {
    return value ~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([Ee][+-]?[0-9]+)?$/
}

function value_after(label,    i) {
    for (i = 1; i <= NF; ++i) {
        if ($i == label) return $(i + 1)
    }
    return ""
}

/HICAR native FGMRES[+]line:/ {
    ++solver_count
    residual = value_after("relative_residual=")
    if (!is_number(residual) || residual + 0 > 1.0e-5) {
        printf "Invalid solver record at line %d: %s\n", NR, $0 > "/dev/stderr"
        bad = 1
    }
}

/HICAR adjoint conservation: relative_Bq=/ {
    ++conservation_count
    constraint = value_after("relative_Bq=")
    target = value_after("target=")
    if (!is_number(constraint) || !is_number(target) ||
        constraint + 0 > target + 0 || target + 0 > 2.0e-5) {
        printf "Invalid conservation record at line %d: %s\n", NR, $0 > "/dev/stderr"
        bad = 1
    }
}

/wind solve rejected|adjoint projection rejected|native FGMRES[+]line failed/ {
    printf "Rejected solver state at line %d: %s\n", NR, $0 > "/dev/stderr"
    bad = 1
}

END {
    if (solver_count == 0) {
        print "No successful FGMRES record found" > "/dev/stderr"
        bad = 1
    }
    if (conservation_count == 0) {
        print "No adjoint conservation record found" > "/dev/stderr"
        bad = 1
    }
    exit bad
}
