#pragma once

namespace thermostat {

    bool has_any_digits(const std::string& s) {
        return std::any_of(s.begin(), s.end(), ::isdigit);
    }

    template <typename T>
    typename std::enable_if<std::is_unsigned<T>::value, int>::type
        inline constexpr signum(T const x) {
        return T(0) < x;
    }

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, int>::type
        inline constexpr signum(T const x) {
        return (T(0) < x) - (x < T(0));
    }

    enum phase { VAPOR , LIQUID };
    enum prop_t { PVT, VLE, VE, SND };
    enum VLE_PHASE { BOTH, LIQ, VAP, NONE };
    enum CALC_STATUS { SUCCESS, FAILED };

    std::map<prop_t, std::string>  prop_t_str = { {PVT,"PVT"},{SND,"SND"} ,{VE,"VE"} ,{VLE,"VLE"} };

    struct compare_settings {
        std::string model_type;
        std::string model_path;
        std::string srk_path;
        std::string BIPcollection;
        std::string depcollection;
        std::string coolprop_root;
        std::vector<std::string> mix;
        std::vector<std::string> dev_files;
        std::string dev_dir;
        std::string combination;
    };

    struct fluid_info {
        std::vector<double> tc, pc, dc, w;
    };

    struct vle_collection {
        std::vector<double> pvap, pliq, xL_0, xV_0;
        std::vector<Eigen::ArrayXd> rhovecL, rhovecV;
    };

    // Collection of relative and relative absolut deviations
    struct dev_collection {
        double rel_abs;
        double rel;
        dev_collection(double rel_abs, double rel) : rel_abs(rel_abs), rel(rel) {};
        dev_collection() {};
    };

    struct data_point {
        double T, P, D, SND;
        Eigen::ArrayXd x;
        std::string author, mix_name;
        int nr;
        int id;
        double dense_v, dense_l, snd_l, snd_v;
        std::vector<double> dense_v_pure;
        std::vector<double> dense_l_pure;
        double x1, x2, y1, y2;
        double ve_v, ve_l;

        // relative and abolsut deviations
        dev_collection deviation;
        dev_collection deviation_vle_liq_p, deviation_vle_vap_p;
        dev_collection deviation_vle_liq_x, deviation_vle_vap_x;

        VLE_PHASE VLE_P;
        CALC_STATUS STATUS;
        CALC_STATUS STAT_VLE_LIQ_P, STAT_VLE_VAP_P;
        CALC_STATUS STAT_VLE_LIQ_X, STAT_VLE_VAP_X;
        data_point() : deviation(0.0, 0.0), deviation_vle_liq_p(0.0, 0.0), deviation_vle_vap_p(0.0, 0.0), deviation_vle_liq_x(0.0, 0.0), deviation_vle_vap_x(0.0, 0.0) {};
    };
}