#pragma once

#include <teqp/types.hpp>
#include <teqp/derivs.hpp>
#include <teqp/algorithms/VLE.hpp>
#include <teqp/models/multifluid_ancillaries.hpp>
#include <teqp/models/cubics.hpp>
#include <teqp/cpp/teqpcpp.hpp>
#include <teqp/json_builder.hpp>
#include <teqp/cpp/deriv_adapter.hpp>
#include <teqp/ideal_eosterms.hpp>
#include <teqp/models/multifluid.hpp>

#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream> 
#include <string>
#include <algorithm>

#include <thermostat/stat_algorithms.hpp>
#include <thermostat/stat_types.hpp>
#include <thermostat/psrk.hpp>


using namespace thermostat;
using namespace teqp;
using namespace std;
using namespace Mie;

namespace thermostat {

    template<typename T, typename TargetT>
    vector<int> findItems(vector<T> const& v, TargetT& target) {
        vector<int> indices;
        auto it = v.begin();
        while ((it = find_if(it, v.end(), [&](T const& e) { return e == target; }))
            != v.end())
        {
            indices.push_back(distance(v.begin(), it));
            it++;
        }
        return indices;
    }

    inline auto getUniqueValues(const vector<double>& inputVector) {
        vector<double> uniqueVector;
        vector<double> sortedVector = inputVector;
        sort(sortedVector.begin(), sortedVector.end());
        sortedVector.erase(unique(sortedVector.begin(), sortedVector.end()), sortedVector.end());
        uniqueVector = sortedVector;
        return uniqueVector;
    }

    template <typename ModelType>
    inline auto check_phase(const ModelType& model, const double& T, const double& rho, const Eigen::ArrayXd& molefrac, const phase& phase_est) {

        // Check if pressure is positive (phase doesnt matter)
        auto is_correct_phase = false;
        auto p = 0.0; auto dpdd = 0.0; auto d2pdd2 = 0.0;
        auto ar01 = 0.0; auto ar02 = 0.0; auto ar03 = 0.0;
        ar01 = model->get_Ar01(T, rho, molefrac);
        ar02 = model->get_Ar02(T, rho, molefrac);
        ar03 = model->get_Ar03(T, rho, molefrac);
        p = (1.0 + ar01) * T * rho * model->get_R(molefrac);
        dpdd = (1.0 + 2.0 * ar01 + ar02) * model->get_R(molefrac) * T;
        d2pdd2 = (2.0 * ar01 + 4.0 * ar02 + ar03) * model->get_R(molefrac) * T / rho;


        if (p > 0.0 && dpdd > 0.0 && d2pdd2 < 0.0 && phase_est == VAPOR) {
            is_correct_phase = true;
        }

        if (!(is_correct_phase) && p > 0.0 && dpdd > 0.0 && d2pdd2 > 0.0 && phase_est == LIQUID) {
            is_correct_phase = true;
        }

        return is_correct_phase;
    }

    template <typename ModelType>
    tuple<double, double, double, double> get_rho_bounds(const ModelType& model, const double& T, const double& rho_est, const Eigen::ArrayXd& molefrac) {

        auto rhomax = 1.100 * rho_est;
        auto rhomin = 0.900 * rho_est;
        auto rhomax_allowed = 2.0 * rho_est;
        auto rhomin_allowed = 0.70 * rho_est;
        auto ptest = 0.0;
        ptest = (1.0 + model->get_Ar01(T, rho_est, molefrac)) * T * rho_est * model->get_R(molefrac);

        if (ptest < 1E-12) {
            rhomax = (rhomin_allowed + rhomax_allowed) / 2.0;
            rhomin = rhomax * 0.80;
        }
        else if (isnan(ptest)) {
            rhomax = rhomax * 0.70;
            rhomin = rhomax * 0.50;
            rhomin_allowed = rhomin * 0.10;
        }
        return  make_tuple(rhomin, rhomax, rhomin_allowed, rhomax_allowed);
    }

    template <typename ModelType>
    tuple<double, double, double, double> get_rho_bounds_extended(const ModelType& model, const double& T, const double& rho_est, const Eigen::ArrayXd& molefrac) {

        auto rhomax = 2.00 * rho_est;
        auto rhomin = 0.500 * rho_est;
        auto rhomax_allowed = 3.0 * rho_est;
        auto rhomin_allowed = 0.2 * rho_est;
        auto ptest = 0.0;
        ptest = (1.0 + model->get_Ar01(T, rho_est, molefrac)) * T * rho_est * model->get_R(molefrac);

        if (ptest < 1E-12) {
            rhomax = (rhomin_allowed + rhomax_allowed) / 2.0;
            rhomin = rhomax * 0.80;
        }
        else if (isnan(ptest)) {
            rhomax = rhomax * 0.70;
            rhomin = rhomax * 0.50;
            rhomin_allowed = rhomin * 0.10;
        }
        return  make_tuple(rhomin, rhomax, rhomin_allowed, rhomax_allowed);
    }

    template <typename ScalarRho>
    auto get_calc_status(const ScalarRho& rho1, const ScalarRho& rho2) {
        auto status = CALC_STATUS::SUCCESS;
        if (rho1 < 0.0 && rho2 < 0.0) {
            status = CALC_STATUS::FAILED;
        }
        return status;
    }

    template <typename ModelType, typename ScalarT, typename ScalarP>
    auto get_rho_from_t_and_p(const ModelType& model, psrk& psrk_s, const ScalarT& T, const ScalarP& P, const Eigen::ArrayXd& molefrac) {


        auto rho_est = 0.0;
        using tdx = TDXDerivatives<ModelType, ScalarT, decltype(molefrac)>;

        // Assume phase is liquid
        phase  phase = LIQUID;

        // Get estimate from pressure srk model
        auto R = model->get_R(molefrac);
        rho_est = psrk_s.get_density(T, P, molefrac, R, phase);


        // Define deviation pressure function
        auto p_diff = [T, P, molefrac, &model](double rho) {
            auto A01 = 0.0;
            A01 = model->get_Ar01(T, rho, molefrac); // model->get_Ar01(T,rho,molefrac)
            return P - (1.0 + A01) * T * rho * model->get_R(molefrac);
        };


        // Initialize the iteration parameters
        auto iter = 0;
        auto max_iter = 50;

        // Get bounds for regular falsi
        auto b = get_rho_bounds(model, T, rho_est, molefrac);
        auto rho_l = 0.0;
        rho_l = regular_falsi(p_diff, get<0>(b), get<1>(b), 1E-8, get<2>(b), get<3>(b), max_iter, iter);
        auto is_correct = check_phase(model, T, rho_l, molefrac, phase);

        if (!is_correct) {
            rho_l = -1000;
        }

        phase = VAPOR;
        rho_est = psrk_s.get_density(T, P, molefrac, model->get_R(molefrac), phase);
        b = get_rho_bounds(model, T, rho_est, molefrac);
        auto rho_v = 0.0;
        rho_v = regular_falsi(p_diff, get<0>(b), get<1>(b), 1E-8, get<2>(b), get<3>(b), max_iter, iter);
        is_correct = check_phase(model, T, rho_v, molefrac, phase);
        if (!(is_correct)) {
            b = get_rho_bounds_extended(model, T, rho_est, molefrac);
            rho_v = regular_falsi(p_diff, get<0>(b), get<1>(b), 1E-8, get<2>(b), get<3>(b), max_iter, iter);
            is_correct = check_phase(model, T, rho_v, molefrac, phase);
            if (!is_correct) {
                rho_v = -1000;
            }
        }
        return make_tuple(rho_v, rho_l, get_calc_status(rho_v, rho_l));
    }

    auto get_dev(const double& one, const double& two) {
        auto dev = one - two;
        auto per_dev = 100 * (dev) / one;
        if (isnan(abs(dev))) { dev = -1000; };
        if (isnan(abs(per_dev))) { per_dev = -1000; };
        return dev_collection(abs(per_dev), per_dev);
    }

    auto find_iso_idx(vector<double>& vec, double& key) {
        auto itr = find(vec.begin(), vec.end(), key);
        return distance(vec.begin(), itr);
    }

    auto find_closest(vector<double>& vec, double& value) {
        assert(!vec.empty());
        auto it = min_element(vec.begin(), vec.end(), [value](double a, double b) {
            return abs(value - a) < abs(value - b);
            });
        assert(it != vec.end());

        return distance(vec.begin(), it);
    }


    inline auto divide_vec(const vector<double>& v, const double& v1) {
        vector<double> results;
        transform(v.begin(), v.end(), back_inserter(results), [&v1](double v) { return v / v1; });
        return results;
    }


    template <typename ModelT, typename RHOT>
    auto get_x_dev(ModelT& model, double& T, double& p, vector<double>& p_curve, RHOT& rhovecL, RHOT& rhovecV) {


        // Get nearest pressure in vle curve to get start values 
        auto idx = find_closest(p_curve, p);

        auto rhovecL_tmp = rhovecL[idx];
        auto rhovecV_tmp = rhovecV[idx];
        // Calculate composition with pressure specified from data point
        auto result = model->mix_VLE_Tp(T, p, rhovecL_tmp, rhovecV_tmp);

        auto xL_0 = result.rhovecL[0] / result.rhovecL.sum();
        auto xV_0 = result.rhovecV[0] / result.rhovecV.sum();

        return make_tuple(xL_0, xV_0);
    }

    template <typename ModelT, typename FracType, typename RHOT>
    auto get_p_sat(ModelT& model, FracType& x_exp, double& T, vector<double>& x_curve, RHOT& rhovecL, RHOT& rhovecV) {

        auto idx = find_closest(x_curve, x_exp);

        // Get nearest composition on vle curve to get start values 
        auto rhovecL_tmp = rhovecL[idx];
        auto rhovecV_tmp = rhovecV[idx];

        auto x = (Eigen::ArrayXd(2) << x_exp, 1.0 - x_exp).finished(); // Mole fractions in the liquid phase (to be kept constant)
        auto [return_code, rhovecLnew, rhovecVnew] = model->mix_VLE_Tx(T, rhovecL[idx], rhovecV[idx], x, 1e-10, 1e-8, 1e-10, 1e-8, 10);


        auto rhoL = rhovecLnew.sum();
        auto rhoV = rhovecVnew.sum();

        auto derivs = model->get_Ar01(T, rhoL, x);
        auto RT = model->get_R(x) * T;
        double pL = rhoL * RT * (1.0 + derivs);

        CALC_STATUS status;

        // check return code if the calculation converged
        switch (return_code)
        {
        case(VLE_return_code::xtol_satisfied):
            status = CALC_STATUS::SUCCESS;
            break;
        case(VLE_return_code::functol_satisfied):
            status = CALC_STATUS::SUCCESS;
            break;
        default:
            status = CALC_STATUS::FAILED;
            break;
        }

        return make_tuple(status, pL);
    }



    class data_set {

    private:
        vector<data_point> points;
        vector<string>     head;
        prop_t prop;
        vector<double> isotherms;
        vector<vle_collection> vle;
        map<int, string> references;
        vector<int> biblio;
    public:

        auto get_reference(int& id) {
            return references[id];
        }

        auto get_points() {
            return points;
        }

        auto get_prop() {
            return prop;
        }

        void insert_point(data_point& point) {
            points.push_back(point);
        }

        void insert_points(vector<data_point> points_in) {
            points = points_in;
        }

        auto get_biblios() {
            return biblio;
        }

        auto get_temperatures() {
            vector<double> T;
            transform(points.begin(), points.end(), back_inserter(T),
                [](data_point const& p) -> double { return p.T; });
            return T;
        }

        auto get_pressures() {
            vector<double> T;
            transform(points.begin(), points.end(), back_inserter(T),
                [](data_point const& p) -> double { return p.P; });
            return T;
        }

        auto get_x0() {
            vector<double> T;
            transform(points.begin(), points.end(), back_inserter(T),
                [](data_point const& p) -> double { return p.x[0]; });
            return T;
        }

        auto get_x0_vle() {
            vector<double> T;
            transform(points.begin(), points.end(), back_inserter(T),
                [](data_point const& p) -> double { return p.x1; });
            return T;
        }

        auto get_y0_vle() {
            vector<double> T;
            transform(points.begin(), points.end(), back_inserter(T),
                [](data_point const& p) -> double { return p.y1; });
            return T;
        }

        auto get_unique_x0() {
            vector<double> T;
            transform(points.begin(), points.end(), back_inserter(T),
                [](data_point const& p) -> double { return p.x[0]; });
            return getUniqueValues(T);
        }


        auto get_ids() {
            vector<double> T;
            transform(points.begin(), points.end(), back_inserter(T),
                [](data_point const& p) -> double { return p.id; });
            return T;
        }

        auto get_points_by_idx(vector<int>& idx) {
            vector<data_point> Result(idx.size());
            auto p = points;
            transform(idx.begin(), idx.end(), Result.begin(), [p](size_t pos) {return p[pos]; });
            return Result;
        }

        void evaluate_head() {
            int aut_start = 11; // author info stats in line 12
            int counter = 0;
            string ref;
            auto nr = 0;
            auto id = 0;

            for (auto line : head) {
                if (counter > aut_start) {
                    stringstream ss;
                    ss << line;
                    ss >> nr >> id >> ref;
                    ref = line.substr(line.find(ref), line.length());
                    references.insert({ id,ref });
                    if (find(biblio.begin(), biblio.end(), id) == biblio.end()) {
                        biblio.push_back(id);
                    }

                    ref = "";
                }
                counter += 1;
            }
        }

        // get dev file PVT for mixtures and add to data points
        void get_pvt_mix_line(string& line, bool& swaped) {
            stringstream ss;
            data_point p;
            double x1 = 0.0; double x2 = 0.0;
            string name_mix, add_info;
            ss << line;

            if (swaped) { ss >> p.P >> p.D >> p.T >> x2 >> x1; }
            else { ss >> p.P >> p.D >> p.T >> x1 >> x2;; }
            stringstream rest_line;
            if (line.length() > 6 && x1 > 0.0 && x2 > 0.0) {
                add_info = line.substr(51, line.length());
                rest_line << add_info.substr(7, add_info.length());
                p.author = add_info.substr(0, 7);
                rest_line >> p.id >> p.mix_name >> p.nr;
                p.x = (Eigen::ArrayXd(2) << x1, x2).finished();
                if (p.P > 0.0 || (p.D > 0.0 && p.T > 0.0)) {
                    points.push_back(p);
                }
            }

        }

        void get_snd_mix_line(string& line, bool& swaped) {
            stringstream ss;
            data_point p;
            double x1 = 0.0; double x2 = 0.0;
            string name_mix, add_info;
            ss << line;
            if (swaped) { ss >> p.P >> p.T >> p.SND >> x2 >> x1; }
            else { ss >> p.P >> p.T >> p.SND >> x1 >> x2; }

            stringstream rest_line;
            if (line.length() > 6) {
                add_info = line.substr(62, line.length());
                rest_line << add_info.substr(7, add_info.length());
                p.author = add_info.substr(0, 7);
                rest_line >> p.id >> p.mix_name >> p.nr;
                p.x = (Eigen::ArrayXd(2) << x1, x2).finished();
                if (p.P > 0.0 || (p.SND > 0.0 && p.T > 0.0)) {
                    points.push_back(p);
                }
            }
        }

        void get_vle_mix_line(string& line, bool& swaped) {
            stringstream ss;
            data_point p;
            double x1 = 0.0; double x2 = 0.0;
            string name_mix, add_info;
            ss << line;
            auto x1_ = 0.0;
            auto x2_ = 0.0;
            auto y1_ = 0.0;
            auto y2_ = 0.0;

            if (swaped) {
                ss >> p.P >> p.T >> p.x2 >> p.x1 >> p.y2 >> p.y1;
            }
            else {
                ss >> p.P >> p.T >> p.x1 >> p.x2 >> p.y1 >> p.y2;
            };


            if (abs(p.x1) > 1E-12 && abs(p.x2) > 1E-12 && abs(p.y1) > 1E-12 && abs(p.y2) > 1E-12) {
                p.VLE_P = VLE_PHASE::BOTH;
            }
            else if (abs(p.x1) > 1E-12 && abs(p.x2) > 1E-12 && abs(p.y1) < 1E-12 && abs(p.y2) < 1E-12)
            {
                p.VLE_P = VLE_PHASE::LIQ;
            }
            else if (abs(p.x1) < 1E-12 && abs(p.x2) < 1E-12 && abs(p.y1) > 1E-12 && abs(p.y2) > 1E-12)
            {
                p.VLE_P = VLE_PHASE::VAP;
            }
            else
            {
                p.VLE_P = VLE_PHASE::NONE;
            }

            stringstream rest_line;
            // IGNORE PURE FLUID POINTS
            if (line.length() > 6 && !(p.VLE_P == VLE_PHASE::NONE)) {
                add_info = line.substr(51, line.length());
                rest_line << add_info.substr(7, add_info.length());
                p.author = add_info.substr(0, 7);
                rest_line >> p.id >> p.mix_name >> p.nr;
                p.x = (Eigen::ArrayXd(2) << x1, x2).finished();
                if (p.P > 0.0 || (p.SND > 0.0 && p.T > 0.0)) {
                    points.push_back(p);
                }
            }
        }

        // read dev file by path and assign data
        auto read_dev_file(string& path, bool& swaped) {
            ifstream file;
            string line;
            file.open(path);
            bool head_end = false;

            prop = PVT;

            if (path.find("PVT") != string::npos) {
                prop = PVT;
            }
            else if (path.find("VE") != string::npos) {
                prop = VE;
            }
            else if (path.find("SND") != string::npos) {
                prop = SND;
            }
            else if (path.find("VLE") != string::npos) {
                prop = VLE;
            }


            while (getline(file, line)) {
                if (line.find("@PCL5") != string::npos) {
                    head.push_back(line);
                    head_end = true;
                }
                else if (!(head_end)) {
                    head.push_back(line);
                }
                else // numbers
                {
                    if (prop == PVT) {
                        get_pvt_mix_line(line, swaped);
                    }
                    else if (prop == VE) {
                        get_pvt_mix_line(line, swaped);
                    }
                    else if (prop == SND) {
                        get_snd_mix_line(line, swaped);
                    }
                    else if (prop == VLE) {
                        get_vle_mix_line(line, swaped);
                    }
                }
            }
            file.close();

            evaluate_head();

            // get isotherms
            vector<double> temperatures, pressures;
            transform(points.begin(), points.end(), back_inserter(temperatures),
                [](data_point const& p) -> double { return p.T; });

            // get corresponding pressures
            transform(points.begin(), points.end(), back_inserter(pressures),
                [](data_point const& p) -> double { return p.P; });

            //FOR VLE GET ALL ISOTHERMS 
            // FOR THE OTHER PROPERTIES ONLY IF pressure is indicated with -1 or -2
            for (size_t i = 0; i < temperatures.size(); ++i) {
                if (count(isotherms.begin(), isotherms.end(), temperatures.at(i)) == 0) {
                    if (prop == VLE) {
                        isotherms.push_back(temperatures.at(i));
                    }
                    else if (!(prop == VLE) && pressures[i] < 0.0) {
                        isotherms.push_back(temperatures.at(i));
                    }


                }
            }

            return prop_t_str[prop];
        }

        // get the p-x relation for all isotherms in the dev file
        template <typename ModelType, typename PureType, typename AncType>
        void get_binary_vle(ModelType& mixture, PureType& pure, AncType& anc, psrk& helper, string& dev_dir) {

            // also check tc values of pures
            auto [tc, dc] = pure[0]->solve_pure_critical(helper.TC[0], helper.DC[0]);

            string vle_out = dev_dir + "/binary_vle.csv";
            ofstream file_vle;
            file_vle.open(vle_out);
            string line;
            line = "p;t;x;symbolsnr;Name;error";
            file_vle << line << endl;
            for (auto t : isotherms) {
                if (t < tc) {
                    //auto [rho_vap, rho_liq] = get_vle_pure_start(get<0>(pure), t, 0);
                    auto rhoV = anc.rhoV(t), rhoL = anc.rhoL(t);
                    auto vle_pure_1 = pure[0]->pure_VLE_T(t, rhoL, rhoV, 50);
                    auto rhovecV0 = (Eigen::ArrayXd(2) << vle_pure_1(1), 0.0).finished();
                    auto rhovecL0 = (Eigen::ArrayXd(2) << vle_pure_1(0), 0.0).finished();
                    auto res = mixture->trace_VLE_isotherm_binary(t, rhovecL0, rhovecV0);

                    // dump the json file for each isotherm
  /*                  string  dump_s;
                    dump_s = res.dump();
                    string s = to_string(t);
                    string name = "vle_" + s + ".json";
                    ofstream out(name);
                    out << dump_s;
                    out.close();*/

                    auto vec_packed = [](const auto& res, const auto& key) {
                        vector<double> x;
                        for (auto res_point : res) {
                            auto val = res_point[key];
                            x.push_back(val);
                        }
                        return x;
                    };

                    auto eigen_pack = [](const auto& res, const auto& key) {
                        vector<Eigen::ArrayXd> x;
                        for (auto res_point : res) {
                            vector<double> val = res_point[key];
                            x.push_back(toeig(val));
                        }
                        return x;
                    };

                    vle_collection vle_i;
                    vle_i.pvap = vec_packed(res, "pV / Pa");
                    vle_i.pliq = vec_packed(res, "pL / Pa");
                    vle_i.xL_0 = vec_packed(res, "xL_0 / mole frac.");
                    vle_i.xV_0 = vec_packed(res, "xV_0 / mole frac.");
                    vle_i.rhovecL = eigen_pack(res, "rhoL / mol/m^3");
                    vle_i.rhovecV = eigen_pack(res, "rhoV / mol/m^3");
                    vle.push_back(vle_i);

                    // write the vle to plotter file

                    for (size_t i = 0; i < vle_i.pliq.size(); i++) {
                        line = "";
                        line = format("{0:10.6f};{1:10.4f};{2:10.6f};{3:10d};{4:10s};{5:10d}", vle_i.pliq[i] / 1E6, t, vle_i.xL_0[i], 100, "This work", 0);
                        file_vle << line << endl;
                    }
                    for (size_t i = 0; i < vle_i.pliq.size(); i++) {
                        line = "";
                        line = format("{0:10.6f};{1:10.4f};{2:10.6f};{3:10d};{4:10s};{5:10d}", vle_i.pvap[i] / 1E6, t, vle_i.xV_0[i], 101, "This work", 0);
                        file_vle << line << endl;
                    }

                }
            }
            file_vle.close();
        }

        // Function for PVT mixture stuff
        // function to calculate the density from temperature and pressure
        template <typename ModelType>
        void pvt_fill(const ModelType& model, psrk& helper) {
            for (size_t i = 0; i < points.size(); i++) {
                auto p_exp = points[i].P * 1E6;
                if (points[i].P == -1.0) { // bubble point density
                    auto idx = find_iso_idx(isotherms, points[i].T);
                    auto pL = 0.0;
                    tie(points[i].STAT_VLE_LIQ_P, pL) = get_p_sat(model, points[i].x[0], points[i].T, vle[idx].xL_0, vle[idx].rhovecL, vle[idx].rhovecV);
                    tie(points[i].dense_v, points[i].dense_l, points[i].STATUS) = get_rho_from_t_and_p(model, helper, points[i].T, pL, points[i].x);
                    points[i].P = pL / 1E6;
                }
                else
                {
                    tie(points[i].dense_v, points[i].dense_l, points[i].STATUS) = get_rho_from_t_and_p(model, helper, points[i].T, p_exp, points[i].x);
                }

            }
        }

        template <typename ModelType, typename PureT>
        void ve_fill(const ModelType& model, const PureT& pures, psrk& helper) {
            for (size_t i = 0; i < points.size(); i++) {
                auto p_exp = points[i].P * 1E6;
                tie(points[i].dense_v, points[i].dense_l, points[i].STATUS) = get_rho_from_t_and_p(model, helper, points[i].T, p_exp, points[i].x);
            }
            auto x_pure = (Eigen::ArrayXd(1) << 1.0).finished();
            for (size_t i = 0; i < points.size(); i++) {
                auto p_exp = points[i].P * 1E6;
                auto [rho1_v, rho1_l, status_1] = get_rho_from_t_and_p(pures[0], helper, points[i].T, p_exp, x_pure);
                auto [rho2_v, rho2_l, status_2] = get_rho_from_t_and_p(pures[1], helper, points[i].T, p_exp, x_pure);
                if (status_1 == CALC_STATUS::SUCCESS) {
                    points[i].ve_v = 1.0 / points[i].dense_v - points[i].x(0) / rho1_v - points[i].x(1) / rho2_v;
                }
                else
                {
                    points[i].ve_v = -1000.0;
                }
                if (status_1 == CALC_STATUS::SUCCESS) {
                    points[i].ve_l = 1.0 / points[i].dense_l - points[i].x(0) / rho1_l - points[i].x(1) / rho2_l;
                }
                else
                {
                    points[i].ve_l = -1000.0;
                }
            }
        }

        template <typename ModelType, typename IdealT>
        void snd_fill(const ModelType& model, const IdealT& ideal, psrk& helper) {
            // First calculate density for all state points
            for (size_t i = 0; i < points.size(); i++) {
                auto p_exp = points[i].P * 1E6;
                tie(points[i].dense_v, points[i].dense_l, points[i].STATUS) = get_rho_from_t_and_p(model, helper, points[i].T, p_exp, points[i].x);
            }
            // get speed of sound for each phase
            auto snd = [&model, ideal](const double& T, const double& rho, const  Eigen::ArrayXd& molefrac, const auto& mw) {
                //using tdx = TDXDerivatives<ModelType, double, vector<double>>;
                auto A01 = model->get_Ar01(T, rho, molefrac);
                auto A02 = model->get_Ar02(T, rho, molefrac);
                auto A11 = model->get_Ar11(T, rho, molefrac);
                auto A20 = model->get_Ar20(T, rho, molefrac);

                auto wih = AlphaCallWrapper<AlphaWrapperOption::idealgas, decltype(ideal)>(ideal);
                using tdxi = TDXDerivatives<IdealT, double, decltype(molefrac)>;
                auto A20_id = 0.0;
                A20_id = tdxi::get_Agenxy<2, 0, ADBackends::autodiff>(wih, T, rho, molefrac);
                auto srq = (1.0 + A01 - A11) * (1.0 + A01 - A11);
                return sqrt((1.0 + 2.0 * A01 + A02 - srq / (A20_id + A20)) * model->get_R(molefrac) / mw * T);
            };
            for (size_t i = 0; i < points.size(); i++) {
                auto mw = 0.0;
                mw = helper.get_mw(points[i].x);
                points[i].snd_l = snd(points[i].T, points[i].dense_l, points[i].x, mw);
                points[i].snd_v = snd(points[i].T, points[i].dense_v, points[i].x, mw);
            }
        }

        // Selector for which data type should be calculated
        template <typename ModelType, typename PureT, typename IdealT>
        void perform_dev_calc(const ModelType& model, const PureT& pures, const IdealT& ideal, psrk& helper, bool& swaped, string& file_out) {
            switch (prop)
            {
            case PVT:
                pvt_fill(model, helper);
                calc_dev_pvt();
                write_pvt_dev(file_out, swaped);
                break;
            case VE:
                ve_fill(model, pures, helper);
                calc_dev_ve();
                write_pvt_dev(file_out, swaped);
                break;
            case SND:
                snd_fill(model, ideal, helper);
                calc_dev_snd();
                write_pvt_dev(file_out, swaped);
                break;
            case VLE:
                calc_dev_vle(model, pures, helper);
                break;
            default:
                break;
            }
        }


        inline auto get_vle_dev(CALC_STATUS& STATUS, double& prop1, double& prop2) {
            dev_collection dev(0.0, 0.0);
            switch (STATUS)
            {
            case(CALC_STATUS::SUCCESS):
                dev = get_dev(prop1, prop2);
                break;
            default:
                break;
            }
            return dev;
        }


        template <typename ModelT, typename PureT>
        void calc_dev_vle(const ModelT& model, const PureT& pures, psrk& helper) {

            vector<double> pc(2);
            auto molefrac = (Eigen::ArrayXd(1) << 1.0).finished();

            // get critical pressure of the pure fluids
            for (auto i = 0; i < 2; i++) {
                auto [tc, dc] = pures[i]->solve_pure_critical(helper.TC[i], helper.DC[i]);
                pc[i] = (1.0 + pures[i]->get_Ar01(tc, dc, molefrac)) * tc * dc * pures[i]->get_R(molefrac);
            }
            // get<0>(pures)->get_Ar01()
            for (size_t i = 0; i < points.size(); i++) {


                // get corresponding vle (p-x) for isotherm
                auto idx = find_iso_idx(isotherms, points[i].T);

                // calculate reduced pressures 
                auto pr_liq = exp(points[i].x1 * log(pc[0]) + points[i].x2 * log(pc[1]));
                auto pr_vap = exp(points[i].y1 * log(pc[0]) + points[i].y2 * log(pc[1]));

                // experimental point reduced
                auto p_exp_red_liq = points[i].P * 1E6 / pr_liq;
                auto p_exp_red_vap = points[i].P * 1E6 / pr_vap;
                auto p_exp = points[i].P * 1E6;

                // Calculate reduced pressure for liquid sat
                auto p_red_liq = divide_vec(vle[idx].pliq, pr_liq);
                auto p_red_vap = divide_vec(vle[idx].pvap, pr_vap);

                auto distance_liq_1 = -1000.0;
                auto distance_liq_2 = -1000.0;
                auto distance_vap_1 = -1000.0;
                auto distance_vap_2 = -1000.0;

                // get xL0 and xV0 for the experimental pressure
                auto [xL0, xV0] = get_x_dev(model, points[i].T, p_exp, vle[idx].pliq, vle[idx].rhovecL, vle[idx].rhovecV);

                // set x satus to SUCCESS // later handle ? to do
                points[i].STAT_VLE_LIQ_X = CALC_STATUS::SUCCESS;
                points[i].STAT_VLE_VAP_X = CALC_STATUS::SUCCESS;

                points[i].deviation_vle_liq_x = get_vle_dev(points[i].STAT_VLE_LIQ_X, points[i].x1, xL0);
                points[i].deviation_vle_vap_x = get_vle_dev(points[i].STAT_VLE_VAP_X, points[i].y1, xV0);



                auto pL = 0.0; auto pV = 0.0;
                // Get bubble / dew pressure for the different phases 
                switch (points[i].VLE_P) {
                case(VLE_PHASE::BOTH):
                    tie(points[i].STAT_VLE_LIQ_P, pL) = get_p_sat(model, points[i].x1, points[i].T, vle[idx].xL_0, vle[idx].rhovecL, vle[idx].rhovecV);
                    tie(points[i].STAT_VLE_VAP_P, pV) = get_p_sat(model, points[i].y1, points[i].T, vle[idx].xV_0, vle[idx].rhovecV, vle[idx].rhovecL);
                    points[i].deviation_vle_liq_p = get_vle_dev(points[i].STAT_VLE_LIQ_P, p_exp, pL);
                    points[i].deviation_vle_vap_p = get_vle_dev(points[i].STAT_VLE_VAP_P, p_exp, pV);
                    break;
                case(VLE_PHASE::LIQ):
                    tie(points[i].STAT_VLE_LIQ_P, pL) = get_p_sat(model, points[i].x1, points[i].T, vle[idx].xL_0, vle[idx].rhovecL, vle[idx].rhovecV);
                    points[i].deviation_vle_liq_p = get_vle_dev(points[i].STAT_VLE_LIQ_P, p_exp, pL);
                    points[i].STAT_VLE_VAP_P = CALC_STATUS::FAILED;
                    break;
                case(VLE_PHASE::VAP):
                    tie(points[i].STAT_VLE_VAP_P, pV) = get_p_sat(model, points[i].y1, points[i].T, vle[idx].xL_0, vle[idx].rhovecL, vle[idx].rhovecV);
                    points[i].deviation_vle_vap_p = get_vle_dev(points[i].STAT_VLE_VAP_P, p_exp, pV);
                    points[i].STAT_VLE_LIQ_P = CALC_STATUS::FAILED;
                    break;
                default:
                    break;
                }
            }
        }

        // calculate deviation of pvt
        void calc_dev_pvt() {
            for (size_t i = 0; i < points.size(); i++) {
                auto dexp = points[i].D * 1E3;
                auto dev1 = get_dev(dexp, points[i].dense_l);
                auto dev2 = get_dev(dexp, points[i].dense_v);
                if (dev1.rel_abs < dev2.rel_abs) {
                    points[i].deviation = dev1;
                }
                else
                {
                    points[i].deviation = dev2;
                }
            }
        }

        void calc_dev_snd() {
            for (size_t i = 0; i < points.size(); i++) {
                auto dev1 = get_dev(points[i].SND, points[i].snd_l);
                auto dev2 = get_dev(points[i].SND, points[i].snd_v);
                if (dev1.rel_abs < dev2.rel_abs) {
                    points[i].deviation = dev1;
                    points[i].D = points[i].dense_l;
                }
                else
                {
                    points[i].deviation = dev2;
                    points[i].D = points[i].dense_v;
                }
            }
        }

        // calculate deviation of excess volume
        void calc_dev_ve() {
            for (size_t i = 0; i < points.size(); i++) {
                auto dexp = points[i].D;
                auto  dev1 = get_dev(dexp, points[i].ve_l);
                auto  dev2 = get_dev(dexp, points[i].ve_v);
                if (dev1.rel_abs < dev2.rel_abs) {
                    points[i].deviation = dev1;
                }
                else
                {
                    points[i].deviation = dev2;
                }
            }
        }

        void write_pvt_dev(string& file_out, bool& swaped) {
            ofstream myfile;
            myfile.open(file_out);
            for (size_t i = 2; i < head.size(); i++) {
                myfile << head[i] << endl;
            }
            for (size_t i = 0; i < points.size(); i++) {
                vector<double> x_;
                for (size_t j = 0; j < points[i].x.size(); j++) {
                    x_.push_back(points[i].x[j]);
                }
                if (swaped) { swap(x_[0], x_[1]); }
                string formatted_str = format("{0:>10.5f}{1:>10.5f}{2:>10.4f}{3:>10.4f}{4:>10.4f}{5:>8s}{6:>6d}{7:>12s}{8:>4d}", points[i].P, points[i].D, points[i].T, x_[0], points[i].deviation.rel, points[i].author, points[i].id, points[i].mix_name, points[i].nr);
                myfile << formatted_str << endl;
            }
            myfile.close();
        }


    };

    template <typename ...Users>
    decltype(auto)
        create_user_group(Users&&... users) {
        return make_tuple(forward<Users>(users)...);
    }


    struct stats_state {
        double tmin, tmax, xmin, xmax, pmin, pmax, AARD;
        int failed_pts, succ_pts;
        stats_state(double AARD, int succ_pts, int failed_pts, double tmin, double tmax, double pmin, double pmax, double xmin, double xmax) : tmin(tmin), tmax(tmax), xmin(xmin), xmax(xmax), pmin(pmin), pmax(pmax), AARD(AARD), succ_pts(succ_pts), failed_pts(failed_pts) {
        }

        auto get_stats_line() {
            string format_ = "{0:40s};{1:18d};{2:18d};";
            vector<double> values;

            int arg_count = 3;
            // Check temperature
            if (tmin == tmax) {
                format_ = format_ + "{" + to_string(arg_count) + ":10.2f};             ";
                arg_count += 1;
                values.push_back(tmin);
            }
            else {
                format_ = format_ + "{" + to_string(arg_count) + ":10.2f} - " + "{" + to_string(arg_count + 1) + ":<10.2f};";
                arg_count += 2;
                values.push_back(tmin);
                values.push_back(tmax);
            }

            if (pmin == pmax) {
                format_ = format_ + "{" + to_string(arg_count) + ":10.2f};             ";
                arg_count += 1;
                values.push_back(pmin);
            }
            else {
                format_ = format_ + "{" + to_string(arg_count) + ":10.2f} - " + "{" + to_string(arg_count + 1) + ":<10.2f};";
                arg_count += 2;
                values.push_back(pmin);
                values.push_back(pmax);
            }

            if (xmin == xmax) {
                format_ = format_ + "{" + to_string(arg_count) + ":10.2f};             ";
                arg_count += 1;
                values.push_back(xmin);
            }
            else {
                format_ = format_ + "{" + to_string(arg_count) + ":10.2f} - " + "{" + to_string(arg_count + 1) + ":<10.2f};";
                arg_count += 2;
                values.push_back(xmin);
                values.push_back(xmax);
            }
            values.push_back(AARD);
            format_ = format_ + "{" + to_string(arg_count) + ":10.2f}";
            string formatted_string;
            if (values.size() == 1) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0]); }
            if (values.size() == 2) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0], values[1]); }
            if (values.size() == 3) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0], values[1], values[2]); }
            if (values.size() == 4) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0], values[1], values[2], values[3]); }
            if (values.size() == 5) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0], values[1], values[2], values[3], values[4]); }
            if (values.size() == 6) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0], values[1], values[2], values[3], values[4], values[5]); }
            if (values.size() == 7) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0], values[1], values[2], values[3], values[4], values[5], values[6]); }
            if (values.size() == 8) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7]); }
            if (values.size() == 9) { formatted_string = format(format_, " ", succ_pts, failed_pts, values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8]); }

            return formatted_string;
        }
    };

    struct stats_state_vle {
        double tmin, tmax, x0min, x0max, y0min, y0max, pmin, pmax, AARD_x, AARD_y;
        int succ_pts;
        stats_state_vle(double AARD_x, double AARD_y, int succ_pts, double tmin, double tmax, double pmin, double pmax, double x0min, double x0max, double y0min, double y0max) :
            tmin(tmin), tmax(tmax), x0min(x0min), x0max(x0max), y0min(y0min), y0max(y0max), pmin(pmin), pmax(pmax), AARD_x(AARD_x), AARD_y(AARD_y), succ_pts(succ_pts) {};

        auto get_stats_line() {

            vector<string> format_;
            vector<double> values;

            auto build_format = [](vector<string>& format, vector<double>& values, vector<double> bounds, bool last = false) {

                if (bounds.size() > 1) {
                    if (bounds[0] == bounds[1]) {
                        if(!(last)){format.push_back("{0:10.2f};             ");}
                        else if(last) { format.push_back("{0:10.2f}              "); };
                        
                        values.push_back(bounds[0]);
                    }
                    else {
                        format.push_back("{0:10.2f} - ");
                        if (!(last)) { format.push_back("{0:<10.2f};"); }
                        else if (last) { format.push_back("{0:<10.2f}"); };
                        values.push_back(bounds[0]);
                        values.push_back(bounds[1]);
                    }
                }
                else
                {
                    values.push_back(bounds[0]);
                    if (!(isnan(bounds[0]))) {
                        if (!(last)) { format.push_back("{0:10.2f};"); }
                        else if (last) { format.push_back("{0:10.2f}"); };
                    }
                    else
                    {
                        if (!(last)) { format.push_back("{0:10s};"); }
                        else if (last) { format.push_back("{0:10s}"); };
                    }
                }
                return make_tuple(format, values);
            };

            tie(format_, values) = build_format(format_, values, { tmin, tmax });
            tie(format_, values) = build_format(format_, values, { pmin, pmax });
            tie(format_, values) = build_format(format_, values, { x0min, x0max });
            tie(format_, values) = build_format(format_, values, { y0min, y0max });
            tie(format_, values) = build_format(format_, values, { AARD_x });
            tie(format_, values) = build_format(format_, values, { AARD_y }, true);


            string formatted_string;
            formatted_string = format("{0:40s};", " ");
            formatted_string += format("{0:18d};", succ_pts);
            for (size_t i = 0; i < values.size(); i++) {
                if (!(isnan(values[i]))) {
                    formatted_string += format(format_[i], values[i]);
                }
                else {
                    formatted_string += format(format_[i], " - ");
                }
            }

            return formatted_string;
        }
    };

    class statistics
    {

    private:
        double AARD, AARD_x, AARD_p, AARD_ortho;
        data_set points_succ, points_fail;
        vector<data_set> author_points;

    public:

        auto get_overall_stats_hom(data_set& points_succ) {
            auto AARD = 0.0;
            auto succ_pts = 0;
            auto failed_pts = 0;
            for (auto p : points_succ.get_points()) {
                if (p.STATUS == CALC_STATUS::SUCCESS) {
                    AARD += p.deviation.rel_abs;
                    succ_pts += 1;
                }
                else
                {
                    failed_pts += 1;
                }

            }
            AARD = AARD / succ_pts;
            // BOUNDS FOR ALL AUTHORS
            // Get min / max temperature
            auto T = points_succ.get_temperatures();
            auto Tmin = *min_element(T.begin(), T.end());
            auto Tmax = *max_element(T.begin(), T.end());

            // Get min / max pressure
            auto p = points_succ.get_pressures();
            auto pmin = *min_element(p.begin(), p.end());
            auto pmax = *max_element(p.begin(), p.end());

            // Get min / max composition
            auto x = points_succ.get_x0();
            auto xmin = *min_element(x.begin(), x.end());
            auto xmax = *max_element(x.begin(), x.end());

            return stats_state(AARD, succ_pts, failed_pts, Tmin, Tmax, pmin, pmax, xmin, xmax);
        }

        auto get_overall_stats_vle(data_set& points) {

            auto AARD_x = 0.0;
            auto AARD_y = 0.0;
            auto succ_pts_liq = 0;
            auto succ_pts_vap = 0;
            // AARDs in pressure
            for (auto p : points.get_points()) {
                if (p.STAT_VLE_LIQ_P == CALC_STATUS::SUCCESS) {
                    AARD_x += p.deviation_vle_liq_p.rel_abs;
                    succ_pts_liq += 1;
                }
                if (p.STAT_VLE_VAP_P == CALC_STATUS::SUCCESS) {
                    AARD_y += p.deviation_vle_vap_p.rel_abs;
                    succ_pts_vap += 1;
                }
            }

            AARD_x = AARD_x / succ_pts_liq;
            AARD_y = AARD_y / succ_pts_vap;

            auto T = points.get_temperatures();
            auto Tmin = *min_element(T.begin(), T.end());
            auto Tmax = *max_element(T.begin(), T.end());

            // Get min / max pressure
            auto p = points.get_pressures();
            auto pmin = *min_element(p.begin(), p.end());
            auto pmax = *max_element(p.begin(), p.end());

            // Get min / max composition
            auto x0 = points.get_x0_vle();
            auto x0min = *min_element(x0.begin(), x0.end());
            auto x0max = *max_element(x0.begin(), x0.end());

            auto y0 = points.get_y0_vle();
            auto y0min = *min_element(y0.begin(), y0.end());
            auto y0max = *max_element(y0.begin(), y0.end());


            return stats_state_vle(AARD_x, AARD_y, succ_pts_liq + succ_pts_vap, Tmin, Tmax, pmin, pmax, x0min, x0max, y0min, y0max);
        }


        void create_stats(data_set& data, ofstream& stats) {


            auto prop = data.get_prop();
            auto points = data.get_points();

            
            stats << "#" << prop_t_str[prop] << endl;
            string formatted_str;
            switch (prop)
            {
            case(VLE):
                formatted_str = "";
                formatted_str = format("{0:40s};{1:18s};{2:26s};{3:26s};{4:26s};{5:26s};{6:26s};{7:26s}", "AUTHOR", "PTSCALC", "TBOUNDS", "PBOUNDS", "XBOUNDS", "YBOUNDS", "AARDY", "AARDX");
                stats << formatted_str << endl;
                for (auto id : data.get_biblios()) {
                    auto idx = findItems(data.get_ids(), id);
                    auto points_aut = data.get_points_by_idx(idx);
                    if (!(points_aut.empty())) {
                        data_set set;
                        set.insert_points(points_aut);
                        // overall statistic (only this one is create for vle)
                        auto reference = data.get_reference(id);
                        stats << format("{0:40s};;;;;;;", reference) << endl;
                        stats_state_vle state_vle = get_overall_stats_vle(set);
                        auto formatted_str = state_vle.get_stats_line();
                        stats << formatted_str << endl;
                    }
                }

                break;
            default:

                // PVT;SND ETC
                formatted_str = "";
                formatted_str = format("{0:40s};{1:18s};{2:18s};{3:26s};{4:26s};{5:26s};{6:26s}", "AUTHOR", "PTSCALC", "PTSCALCFAILED", "TBOUNDS", "PBOUNDS", "XBOUNDS", "AARD");
                stats << formatted_str << endl;

                for (auto id : data.get_biblios()) {
                    // get idx of the data points of the corrsponding id
                    auto idx = findItems(data.get_ids(), id);

                    // get points of author
                    auto points_aut = data.get_points_by_idx(idx);

                    data_set set;
                    set.insert_points(points_aut);

                    // overall statistic
                    auto reference = data.get_reference(id);
                    stats << format("{0:40s};;;;;;", reference) << endl;
                    stats_state state = get_overall_stats_hom(set);
                    auto formatted_str = state.get_stats_line();
                    stats << formatted_str << endl;

                    auto x_uniques = set.get_unique_x0();

                    // Structure for all different compsitions for an author
                    if (x_uniques.size() > 1) {
                        for (auto x_ : x_uniques) {
                            // get index of all points of that author and the compositon
                            auto idx1 = findItems(set.get_x0(), x_);
                            auto points_aut_x = set.get_points_by_idx(idx1);
                            data_set set_aut_x;
                            set_aut_x.insert_points(points_aut_x);
                            stats_state state = get_overall_stats_hom(set_aut_x);
                            auto formatted_str = state.get_stats_line();
                            stats << formatted_str << endl;
                            auto test = 0.0;
                        }
                    }
                }
                break;
            }



        }
    };

    // get ancillary equations form multi fluid model
    auto get_anc(compare_settings& s) {
        const auto mixture_model = build_multifluid_model(s.mix, s.coolprop_root, s.BIPcollection, s.depcollection);
        auto jj = nlohmann::json::parse(mixture_model.get_meta());
        auto ancillaries = jj.at("pures")[0].at("ANCILLARIES");
        teqp::MultiFluidVLEAncillaries anc(ancillaries);
        return anc;
    }

    // Create ideal fluid model (is always the on in the Coolprop GERG files)
    auto get_ideal_fluid_model(compare_settings& s) {
        auto pure_data = nlohmann::json::array();
        for (auto f : s.mix) {
            string path = s.coolprop_root + "/dev/fluids/" + f + ".json";
            auto j = convert_CoolProp_idealgas(path, 0);
            pure_data.push_back(j);
        }
        IdealHelmholtz mix_ideal(pure_data);
        return mix_ideal;
    }

    // Function get the elongated model
    auto make_elong(const nlohmann::json& spec) {
        using namespace teqp::cppinterface::adapter;
        return make_owned(Mie::MieElong(spec.at("model_path"), spec.at("components"), spec.at("combination")));
    }

    // Function to get pure fluid in std containter
    auto make_pures(const nlohmann::json& j) {

        // Get component names from mixture
        auto spec = j.at("model");
        vector<string> flds = spec.at("components");

        // create copys of the json and get pure fluid as single models 
        vector< nlohmann::json> pure_j;
        vector<unique_ptr<teqp::cppinterface::AbstractModel>> pures;

        for (auto i = 0; i < flds.size(); i++)
        {
            pure_j.push_back(j);
            pure_j[i]["model"].at("components") = { flds[i] };
            pures.emplace_back(((pure_j[i]["kind"] != "MieElong") ? make_model(pure_j[i]) : make_elong(pure_j[i]["model"])));
        }
        return pures;
    }

    inline auto get_compare_settings(string& path) {
        string filepath = filesystem::is_regular_file(path) ? path : path;
        nlohmann::json j = load_a_JSON_file(filepath);
        compare_settings set;
        auto spec = j.at("model");
        set.model_path = spec.at("model_path");
        set.srk_path = spec.at("srk_path");
        set.BIPcollection = spec.at("BIP");
        set.depcollection = spec.at("departure");
        set.coolprop_root = spec.at("coolprop_root");
        vector<string> mixture = spec.at("components");
        vector<string> dev_files = spec.at("dev_files");
        set.dev_files = dev_files;
        set.mix = mixture;
        set.combination = spec.at("combination");
        set.dev_dir = spec.at("dev_dir");
        return set;
    }


    template<typename MixType, typename PureType>
    auto call_tester(MixType& mix, PureType& pure) {
        auto T = 200.0;
        auto rho = 13000.0;
        auto molefrac = (Eigen::ArrayXd(2) << 0.5, 0.5).finished();
        auto res = mix->get_Ar00(T, rho, molefrac);
        auto molefrac1 = (Eigen::ArrayXd(1) << 1.0).finished();
        auto res2 = pure[0]->get_Ar00(T, rho, molefrac1);

        auto rhovecV0 = (Eigen::ArrayXd(2) << 1.0, 0.0).finished();
        auto rhovecL0 = (Eigen::ArrayXd(2) << 1.0, 0.0).finished();
        auto res_vle = mix->trace_VLE_isotherm_binary(T, rhovecL0, rhovecV0);


        return 0.0;
    }


    auto check_swap_pures(psrk& helper) {
        auto swap = false;
        for (size_t i = 0; i < helper.TC.size() - 1; i++) {
            if (helper.TC[i] < helper.TC[i + 1]) { // first fluids tc is lower than second one, swap!
                swap = true;
            }
        }
        return swap;
    }

    // Simply swap the fluids in the json
    auto swap_fluids(nlohmann::json& j, compare_settings& set, bool& swapped, psrk& anc) {
        if (swapped) {
            swap(set.mix[0], set.mix[1]);
            anc.swap_content();
            j["model"].at("components") = { set.mix[0],set.mix[1] }; // swap json data
        }
        return j;
    }

    auto create_out_name(string& path, string& components, string& name) {
        auto s = components.substr(0, components.find("_"));
        return path + "/" + s + "_" + name + ".plt";
    }

    auto create_statistics(std::string& set_path) {

        nlohmann::json jj;

        // Read json for fluid and model information
        try {
            jj = load_a_JSON_file(set_path);
        }
        catch ( ... ) {
            throw;
        }

        auto settings = get_compare_settings(set_path);
        auto helper = psrk(settings);
        auto flds_swap = check_swap_pures(helper);
        auto j = swap_fluids(jj, settings, flds_swap, helper);
        auto pures = make_pures(j);
        auto ancs = get_anc(settings);
        auto mixture = (j["kind"] != "MieElong") ? make_model(j) : make_elong(j.at("model"));
        auto mix_ideal = get_ideal_fluid_model(settings);


        string stats_out = settings.dev_dir + "/" + "stats.txt";
        ofstream stats;
        stats.open(stats_out);
        vector< vector<data_point>> data;
        for (size_t i = 0; i < settings.dev_files.size(); i++) {
            std::string t_name = settings.dev_files[i];
            std::string ser = ".DEV";
            std::string::size_type k = t_name.find(ser);
            t_name.erase(k,ser.length());
            stats << "@";
            stats << t_name << endl;
            data_set dat;
            string dev_path = settings.dev_dir + "/" + settings.dev_files[i];
            auto prop_name = dat.read_dev_file(dev_path, flds_swap);
            auto file_out = create_out_name(settings.dev_dir, settings.dev_files[i], prop_name);
            dat.get_binary_vle(mixture, pures, ancs, helper, settings.dev_dir);
            dat.perform_dev_calc(mixture, pures, mix_ideal, helper, flds_swap, file_out);
            statistics stat;
            stat.create_stats(dat, stats);
        }
        stats << '#' << endl;
        stats.close();

        return 0;
    }


}