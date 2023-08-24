#pragma once
#include <stat_types.hpp>
#include <psrk.hpp>

namespace thermostat {

    template<typename Callable>
    inline auto regular_falsi(Callable zerofunction, const double& Xstart_min, const double& Xstart_max, const double& Delta_Allowed, const double& Xmin_allowed, const double& Xmax_allowed, int& Max_iterations, int& Iterations) {

        auto Xmin = Xstart_min;
        auto Xmax = Xstart_max;
        auto Xmin_alt = 0.0;
        auto Xmax_alt = 0.0;
        auto F_Xmin_alt = 0.0;
        auto F_Xmax_alt = 0.0;
        auto F_Xmin = 0.0;
        auto F_Xmax = 0.0;
        auto Stop_min = 0;
        auto Stop_max = 0;
        auto Int_Errorflag = 0;
        bool root_in = false;
        bool Root_Found = false;
        bool additional_checks = false;
        auto Count_enlarge = 0;
        auto Count_reduce = 0;
        auto result = 0.0;
        auto Xakt = 0.0;
        auto DERIV = 0.0;
        auto F_Xakt = 0.0;
        auto Deviation_X = 0.0;
        auto Xroot = 0.0;
        auto Xtest = 0.0;
        auto  F_Xtest = 0.0;
        auto xmin = 0.0;
        auto xmax = 0.0;
        auto xakt = 0.0;

        if (Xmin < Xmin_allowed) { Int_Errorflag = 1; }
        if (Xmin < Xmin_allowed) { Xmin = Xmin_allowed; }
        if (Xmax > Xmax_allowed) { Int_Errorflag = Int_Errorflag + 1; }
        if (Xmax > Xmax_allowed) { Xmax = Xmax_allowed; }
        if (Int_Errorflag == 1) { Int_Errorflag = 0; }
        if (Int_Errorflag == 2) { Int_Errorflag = 1; }

        if (Int_Errorflag == 0) {
            F_Xmin = zerofunction(Xmin);
        }

        if (Int_Errorflag == 0) {
            F_Xmax = zerofunction(Xmax);
        }

        if ((F_Xmin * F_Xmax) < 0) { root_in = true; }

        // main loop
        while (!(root_in) && Int_Errorflag == 0 && Count_enlarge <= 10) {
            Xmin_alt = Xmin;
            Xmax_alt = Xmax;
            F_Xmin_alt = F_Xmin;
            F_Xmax_alt = F_Xmax;
            if (Stop_min == 0) { Xmin = Xmin - (Xstart_max - Xstart_min) * 0.5; }
            if (Stop_max == 0) { Xmax = Xmax + (Xstart_max - Xstart_min) * 0.5; }
            Count_enlarge = Count_enlarge + 1;
            if ((Xmin - Xmin_allowed) < -1e-11) { Int_Errorflag = 1; }
            if ((Xmax - Xmax_allowed) > 1e-11) { Int_Errorflag = 1; }
            if ((Int_Errorflag == 0) && (Stop_min == 0)) { F_Xmin = zerofunction(Xmin); }
            if ((Int_Errorflag == 0) && (Stop_min == 0)) { F_Xmax = zerofunction(Xmax); }
            if ((abs(F_Xmin_alt) < abs(F_Xmin)) && ((F_Xmin_alt * F_Xmin) > 0)) { Stop_min = 1; }
            if (Stop_min == 1) { Xmin = Xmin_alt; }
            if (Stop_min == 1) { F_Xmin = F_Xmin_alt; }
            if ((abs(F_Xmax_alt) < abs(F_Xmax)) && ((F_Xmax_alt * F_Xmax) > 0)) { Stop_max = 1; }
            if (Stop_max == 1) { Xmax = Xmax_alt; }
            if (Stop_max == 1) { F_Xmax = F_Xmax_alt; }
            if ((F_Xmin * F_Xmax) < 0) { root_in = true; }
        }

        // If no root was found by enlarging the interval, there may be two roots within the initial interval.
        // To check out this posibility, the original interval is subdivided starting in the midle of the original interval.
        // The size of the new interval starts with 1 / 10 of the original interval and is enlarged then.

        while (!(root_in) && (Int_Errorflag == 0) && (Count_reduce < 9)) {
            Count_reduce = Count_reduce + 1;
            Xmin = (Xstart_max + Xstart_min) / 2.0 - (Xstart_max - Xstart_min) / 20.0 * Count_reduce;
            Xmax = (Xstart_max + Xstart_min) / 2.0 + (Xstart_max - Xstart_min) / 20.0 * Count_reduce;
            //Check, whether Xmin and Xmax get out of the allowed range
            if (Xmin < Xmin_allowed) { Int_Errorflag = 1; }
            if (Xmax > Xmax_allowed) { Int_Errorflag = 1; }
            if (Int_Errorflag == 0) { F_Xmin = zerofunction(Xmin); }
            if (Int_Errorflag == 0) { F_Xmax = zerofunction(Xmax); }
            //Check, whether a root is within the new limits
            if ((F_Xmin * F_Xmax) < 0.0) { root_in = true; }
        }

        //If no root could be found in the original interval, by enlarging or by reducing of the interval the
        //routine returns with Errorflag = 2

        if (!(root_in)) { result = -2; }

        if (result == 0) { // no error occured 


            while (!(Root_Found) && (Iterations <= Max_iterations)) {

                Iterations = Iterations + 1; //Count iteration steps
                auto Deriv = (F_Xmax - F_Xmin) / (Xmax - Xmin);  //Calculate the "derivative" for the current interval;

                if (abs(Deriv) > 1E-12) {
                    Deviation_X = F_Xakt / Deriv;//Deviation_X is the distance to the linear approximated root from the akt position on the x - axis
                    //SH 10 / 2018: If derivative is large make additional checks(see line 300) if root found or not
                    if (abs(Deriv) > 1E10) {
                        additional_checks = true;
                    }
                    else {
                        additional_checks = false;
                    }
                }
                else {
                    Deriv = signum(Deriv);
                }

                //Calculate Xakt starting from the X value with the smaller function value
                if (abs(F_Xmin) < abs(F_Xmax)) {
                    // Calculate Xakt starting from the lower limit of the interval
                    if (abs(Deriv) > 0 && Iterations % 5 != 0) {
                        Xakt = Xmin + abs(F_Xmin / Deriv);
                        // Use this algorithm only three out of four times to avoid problems in regions  // with strong curvature
                    }
                    else {
                        Xakt = (Xmin + Xmax) / 2.;
                    }
                }
                else {
                    if (abs(Deriv) > 0 && Iterations % 5 != 0) {            // Calculate Xakt starting from the upper limit of the interval
                        Xakt = Xmax - abs(F_Xmax / Deriv);
                    }                                 // Use this algorithm only three out of four times to avoid problems in regions
                    else {                                                               // with strong curvature
                        Xakt = (Xmin + Xmax) / 2.0;
                    }
                }


                F_Xakt = zerofunction(Xakt);


                //Check whether the root has been found
                //The criterion Delta_allowed is applied to the deviation in X

                if (abs(Deriv) > 1E-12) {
                    Deviation_X = F_Xakt / Deriv;
                    if (abs(Deriv) > 1E10) {
                        additional_checks = true;
                    }
                    else {
                        additional_checks = false;
                    }
                }
                else {
                    Deriv = signum(Deriv);
                }

                if (((abs(Deviation_X) < Delta_Allowed) && (!(additional_checks))) || ((abs(Deviation_X) < Delta_Allowed) && (additional_checks) && ((abs((xakt - xmin) / xakt) < Delta_Allowed) || (abs((xakt - xmax) / xakt) < Delta_Allowed)))) {
                    Root_Found = true;
                    Xroot = Xakt + Deviation_X;  // Xroot usually is a better guess for the root then Xakt
                }
                else  //Root has not yet been found, new interval has to be defined
                {
                    if (abs(Deviation_X) < ((Xmax - Xmin) / 10.0)) {

                        if (((Xakt - Xmin) < (Xmax - Xakt)) && ((Xakt - Xmin) > 0.0)) { Deriv = (F_Xakt - F_Xmin) / (Xakt - Xmin); }
                        if (((Xakt - Xmin) > (Xmax - Xakt)) && ((Xmax - Xakt) > 0.0)) { Deriv = (F_Xmax - F_Xakt) / (Xmax - Xakt); }

                        if (abs(Deriv) > 1E-12) {
                            Deviation_X = F_Xakt / Deriv;
                        }
                        else {
                            Deriv = signum(Deriv);
                        }

                        if ((F_Xmin * F_Xakt) < 0.0) {
                            Xtest = Xakt - 2.0 * abs(Deviation_X);
                            if (Xtest < ((Xmin + Xakt) / 2.0)) { Xtest = (Xmin + Xakt) / 2.0; }
                            F_Xtest = zerofunction(Xtest);
                            if ((F_Xtest * F_Xakt) < 0.0) { Xmin = Xtest; }
                            if ((F_Xtest * F_Xakt) < 0.0) { F_Xmin = F_Xtest; }
                        }
                        else if ((F_Xmax * F_Xakt) < 0.0) {
                            Xtest = Xakt + 2.0 * abs(Deviation_X);
                            if (Xtest > ((Xmax + Xakt) / 2.0)) { Xtest = (Xmax + Xakt) / 2.0; }
                            F_Xtest = zerofunction(Xtest);
                            if ((F_Xtest * F_Xakt) < 0.0) { Xmax = Xtest; }
                            if ((F_Xtest * F_Xakt) < 0.0) { F_Xmax = F_Xtest; }
                        }
                    }

                    if ((F_Xmin * F_Xakt) < 0.0) {
                        Xmax = Xakt;
                        F_Xmax = F_Xakt;
                    }
                    else if ((F_Xmax * F_Xakt) < 0.0) {
                        Xmin = Xakt;
                        F_Xmin = F_Xakt;
                    }
                }

            } // end while
        } //end if

        if ((Count_enlarge + Count_reduce) > 0) { Xroot = -1000; }
        if (!(Root_Found)) { Xroot = -1000; }

        return Xroot;

    }


    template <typename TypeT, typename TypeX, typename TypeK>
    inline auto rache_rich(psrk& fldi, TypeK& K_val, TypeX& x_known, TypeT& vapfrac, int errval) {
        auto ncomp = x_known.size();
        auto sumx_bubble = 0.0;
        auto sumx_dew = 0.0;
        auto sumx_2phase = 0.0;
        std::vector<double> x_vap(ncomp);
        std::vector<double> x_liq(ncomp);
        bool found = false;
        auto frac_type = 0;
        auto Delta_allowed = 1E-8;
        auto frac_min_allowed = -1E-16;
        auto frac_min = -1E-16;
        auto frac_max = 0.5;
        auto frac_max_allowed = 1.0 + 1E-16;
        auto Max_iterations = 50;
        auto frac = 0.0;

        auto rache_diff = [frac_type, ncomp, K_val, x_known](double frac) {
            auto  rac_func = 0.0;
            if (frac_type == 0) {
                for (size_t i = 0; i < ncomp; i++) {
                    rac_func = rac_func + x_known[i] * (K_val[i] - 1.0) / (1.0 - frac + frac * K_val[i]);
                }
            }
            else {
                for (size_t i = 0; i < ncomp; i++)
                {
                    rac_func = rac_func + x_known[i] * (K_val[i] - 1.0) / (frac + (1.0 - frac) * K_val[i]);
                }
            }
            return rac_func;
        };

        for (size_t i = 0; i < ncomp; i++) {
            sumx_bubble = sumx_bubble + x_known[i] * K_val[i];
            sumx_dew = sumx_dew + x_known[i] / K_val[i];
            sumx_2phase = sumx_2phase + x_known[i] * (K_val[i] - 1.0) / (K_val[i] + 1.0);
        }

        if (sumx_bubble < 1.0 + 1E-15) {
            vapfrac = 0;
            for (size_t i = 0; i < ncomp; i++) {
                x_vap[i] = x_known[i] * K_val[i] / sumx_bubble;
            }
            x_liq = x_known;
            found = true;
        }

        if (sumx_dew < 1.0 + 1E-15) {
            vapfrac = 0;
            for (size_t i = 0; i < ncomp; i++) {
                x_liq[i] = x_known[i] / K_val[i] / sumx_dew;
            }
            x_vap = x_known;
            found = true;
        }

        //The mixture is in the 2 phase region with 0.5 < vapfrac < 1.
        //Instead of vapfrac the liquid fraction alpha is used, since vapfrac near unity might cause round off errors(see Michelsen& Mollerup)
        if (!(found)) {
            //The liquid fraction is chosen as independent variable
            if (sumx_2phase > 0.0) {
                frac_type = 1;
            }
            else {
                frac_type = 0;
            }

            auto Iterations = 0;
            frac = 0.0;
            auto res = regular_falsi(rache_diff, frac_min, frac_max, Delta_allowed, frac_min_allowed, frac_max_allowed, Max_iterations, Iterations);
            if (errval == 3) { errval = 0; };
            if (errval != 0) { errval = -2111; }
            if (frac_type == 1) {
                vapfrac = 1.0 - frac;
            }
            else {
                vapfrac = frac;
            }

            sumx_dew = 0.0;
            sumx_bubble = 0.0;
            for (size_t i = 0; i < x_known.size(); i++) {
                x_vap[i] = x_known[i] * K_val[i] / (1.0 - vapfrac + vapfrac * K_val[i]);
                sumx_dew = sumx_dew + x_vap[i];
                x_liq[i] = x_known[i] / (1.0 - vapfrac + vapfrac * K_val[i]);
                sumx_bubble = sumx_bubble + x_liq[i];
            }

            for (size_t i = 0; i < x_vap.size(); i++)
            {
                x_vap[i] = x_vap[i] / sumx_dew;
                x_liq[i] = x_liq[i] / sumx_bubble;
            }


        }
        return std::make_tuple(x_vap, x_liq, vapfrac);
    }
    template <typename TypeT, typename TypeP, typename TypeX>
    inline auto PTX_startvals_pT(psrk& fldi, TypeT& T, TypeP p, TypeX& x_known, TypeT& vapfrac, int errval) {

        std::vector<double>  K_val(x_known.size());
        for (size_t i = 0; i < x_known.size(); i++) {
            K_val[i] = fldi.PC[i] / p * exp(5.3730 * (1.0 + fldi.w[i]) * (1.0 - fldi.TC[i] / T));
        }

        return  rache_rich(fldi, K_val, x_known, vapfrac, errval);

    }
    template <typename TypeT, typename TypeP>
    inline auto PhaseDet(psrk& fldi, TypeT& T, TypeP& p, const std::vector<double>& x) {

        auto errval = 0;
        auto vapfrac = 0.0;
        std::vector<TypeT> x_vap, x_liq2, x_liq;
        tie(x_vap, x_liq, vapfrac) = PTX_startvals_pT(fldi, T, p, x, vapfrac, errval);
        return 0.0;
    }

}