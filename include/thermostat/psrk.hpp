#pragma once
#include "teqp/json_tools.hpp"
#include "thermostat/stat_types.hpp"


// Implementation of the Soave-Redlich-Kwong-EOS in pressure
namespace thermostat {

    class psrk {

    public:
        std::vector<double> TC, PC, DC, w , mw;
        psrk(compare_settings& s) {
            std::string filepath = std::filesystem::is_regular_file(s.srk_path) ? s.srk_path : s.srk_path;
            nlohmann::json j = teqp::load_a_JSON_file(filepath);
            for (size_t i = 0; i < s.mix.size(); i++) {
                TC.push_back(j.at("srk").at(s.mix[i]).at("tc"));
                PC.push_back(j.at("srk").at(s.mix[i]).at("pc"));
                DC.push_back(j.at("srk").at(s.mix[i]).at("dc"));
                w.push_back(j.at("srk").at(s.mix[i]).at("omega"));
                mw.push_back(j.at("srk").at(s.mix[i]).at("mw"));
            }
        }
        
        void swap_content() {
            for (size_t i = 0; i < TC.size()-1; i++) {
                std::swap(TC[i], TC[i+1]);
                std::swap(PC[i], PC[i + 1]);
                std::swap(DC[i], DC[i + 1]);
                std::swap(w[i], w[i + 1]);
                std::swap(mw[i], mw[i + 1]);
            }
            
        }
        
        auto get_mw(Eigen::ArrayXd& molefrac) {
            auto m = 0.0;
            for (size_t i = 0; i < molefrac.size(); i++)
            {
                m = m + molefrac(i) * mw[i];
            }
            return m;
        }

        std::vector<double> root_psrk(std::vector<double>& p) {

            auto pie = 3.14159265358979;
            std::vector<double> r(3);
            auto a = p[0];
            auto b = p[1];
            auto c = p[2];
            auto nr = 0;
            auto nn = 0.0;
            auto qq = 0.0;
            auto pp = 0.0;
            auto aa = 0.0;
            auto bb = 0.0;

            pp = (a * a - 3. * b) / 9.;
            qq = (2. * pow(a, 3.0) - 9. * a * b + 27. * c) / 54.0;
            if (pow(pp, 3.0) > pow(qq, 2.0)) {
                nn = acos(-qq / pow(pp, 1.5));
                nr = 3;
                r[0] = 2. * sqrt(pp) * cos(nn / 3.) - a / 3.;
                r[1] = 2. * sqrt(pp) * cos((nn + 2. * pie) / 3.) - a / 3.;
                r[2] = 2. * sqrt(pp) * cos((nn + 4. * pie) / 3.) - a / 3.;
            }
            else {
                nr = -1;
                aa = -signum(qq) * pow((abs(qq) + sqrt(pow(qq, 2.0) - pow(pp, 3.0))), 1.0 / 3.0);
                if (aa == 0) {
                    bb = 0.0;
                }
                else {
                    bb = pp / aa;
                }

                r[0] = (aa + bb) - a / 3.0;
                r[1] = 0.0;
                r[2] = 0.0;
            }

            // --creating output vector------------------------------------------------------------------ -
            std::vector<double> volume(3);
            if (nr < 0) {
                volume[0] = -1.0;
                volume[1] = r[0];
                volume[2] = 0.0;
            }


            else {
                volume[0] = 3.0;
                auto idx = max_element(std::begin(r), std::end(r));
                volume[1] = r[*idx];
                auto qq = 1.0;
                pp = 1.0;
                nn = 1.0;
                if (r[0] > 0) { qq = r[0]; }

                if (r[1] > 0) { pp = r[1]; };
                if (r[2] > 0) { nn = r[2]; };
                std::vector<double> r_vec = { qq, pp, nn };
                idx = min_element(std::begin(r_vec), std::end(r_vec));
                volume[2] = *idx;
            }

            return volume;
        }



        auto get_density(const double& T, const double& P, const Eigen::ArrayXd& molefrac, const double& REQ, phase& phase_est) {
            auto ncomp = molefrac.size();
            std::vector<double> aa(ncomp), bb(ncomp), cc(3), msrk(ncomp), alfa_lc(ncomp);
            std::vector<double> ai(ncomp), bi(ncomp), ci(ncomp);
            auto a1 = 0.0;
            auto c = 0.0;
            auto b = 0.0;
            for (auto i = 0; i < molefrac.size(); i++) {
                auto tr = T / TC[i];
                aa[i] = 0.42748 * REQ * REQ * TC[i] * TC[i] / PC[i];
                msrk[i] = 0.48 + 1.574 * w[i] - 0.176 * pow(w[i], 2.0);
                alfa_lc[i] = pow(1.0 + msrk[i] * (1.0 - sqrt(tr)), 2.0);
                ai[i] = aa[i] * alfa_lc[i];
                bi[i] = 0.08664 * REQ * TC[i] / PC[i];
                ci[i] = 0.40768 * REQ * TC[i] / PC[i] * (0.29441 - PC[i] * 1.0 / DC[i] / REQ / TC[i]);
                a1 = a1 + molefrac(i) * sqrt(ai[i]);
                b = b + molefrac(i) * bi[i];
                c = c + molefrac(i) * ci[i];
            }

            auto a = 0.0;
            if (ncomp == 1) {
                a = ai[0];
                b = bi[0];
                c = ci[0];
                cc[0] = -REQ * T / P;
                cc[1] = -b * b - REQ * T * b / P + a / P;
                cc[2] = -a * b / P;
            }
            else {
                a = pow(a1, 2.0);
                cc[0] = -REQ * T / P;
                cc[1] = -b * b - REQ * T * b / P + a / P; //rmix
                cc[2] = -a * b / P;
            }



            std::vector<double> volume = root_psrk(cc);

            auto psrk_val = 0.0;
            auto VSRK = 0.0;
            if (volume[0] < -10.0) {
                psrk_val = -1.0;
            }
            else {
                if (volume[0] < -0.0) {
                    VSRK = 1. / (volume[1]);
                }
                else {
                    if (phase_est == LIQUID) {
                        VSRK = 1. / (volume[2]);
                    }
                    else {
                        VSRK = 1. / (volume[1]);
                    }
                }
            }

            auto vf1 = 0.0;
            if (phase_est == LIQUID) {
                vf1 = 1.0 / VSRK;
                vf1 = vf1 - c;
                VSRK = 1 / vf1;
            }

            psrk_val = VSRK;
            return psrk_val;
        }

    };

}