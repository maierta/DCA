#!/usr/bin/awk -f
# Generates the tc/preconfigured input files from the templates input_sp.json.in and input_tp.json.in.
# Run this script from tutorials/tc/.

BEGIN {
    # Temperature series from high to low.
    n_temps = 9
    temps[1]  = 1
    temps[2]  = 0.75
    temps[3]  = 0.5
    temps[4]  = 0.25
    temps[5]  = 0.125
    temps[6]  = 0.1
    temps[7]  = 0.09
    temps[8]  = 0.08
    temps[9]  = 0.07

    # Pre-computed beta = 1/T values matching the existing preconfigured files.
    betas[1]  = 1
    betas[2]  = 1.33333
    betas[3]  = 2
    betas[4]  = 4
    betas[5]  = 8
    betas[6]  = 10
    betas[7]  = 11.1111
    betas[8]  = 12.5
    betas[9]  = 14.2857

    # Temperatures that also need a two-particle (tp) input file.
    has_tp[6] = 1
    has_tp[7] = 1
    has_tp[8] = 1
    has_tp[9] = 1

    # Fixed physical parameters.
    dens = 0.95
    hubbardU = 6
    vec1 = "[2, 0]"
    vec2 = "[0, 2]"

    for (i = 1; i <= n_temps; i++) {
        t = temps[i]

        # Previous temperature for self-energy restart ("zero" for the highest T).
        prev = (i == 1) ? "zero" : temps[i - 1]

        # DCA iterations: 8 for the first temperature, 6 otherwise.
        iters_sp = (i == 1) ? 8 : 6

        dir = "preconfigured/T=" t
        cmd = "mkdir -p " dir
        system(cmd)

        # Generate input_sp.json.
        sp_in = "input_sp.json.in"
        sp_out = dir "/input_sp.json"
        while ((getline line < sp_in) > 0) {
            gsub(/BETA/, betas[i], line)
            gsub(/DENS/, dens, line)
            gsub(/HUBBARDU/, hubbardU, line)
            gsub(/VEC1/, vec1, line)
            gsub(/VEC2/, vec2, line)
            gsub(/CURRENT_TEMP/, t, line)
            gsub(/PREVIOUS_TEMP/, prev, line)
            # Fix the initial self-energy for the highest temperature.
            gsub(/"\.\/T=zero\/dca_sp\.hdf5"/, "\"zero\"", line)
            gsub(/ITERS/, iters_sp, line)
            print line > sp_out
        }
        close(sp_in)
        close(sp_out)

        # Generate input_tp.json for temperatures that need it.
        if (has_tp[i]) {
            tp_in = "input_tp.json.in"
            tp_out = dir "/input_tp.json"
            while ((getline line < tp_in) > 0) {
                gsub(/BETA/, betas[i], line)
                gsub(/DENS/, dens, line)
                gsub(/HUBBARDU/, hubbardU, line)
                gsub(/VEC1/, vec1, line)
                gsub(/VEC2/, vec2, line)
                gsub(/CURRENT_TEMP/, t, line)
                # PREVIOUS_TEMP is not used in the tp template, but substitute
                # with the current temperature for safety.
                gsub(/PREVIOUS_TEMP/, t, line)
                gsub(/ITERS/, 1, line)
                print line > tp_out
            }
            close(tp_in)
            close(tp_out)
        }
    }
}
