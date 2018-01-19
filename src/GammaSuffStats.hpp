#include <vector>
#include "Eigen/Dense"
#include "Random.hpp"
#include "doctest.h"
using DMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using AAProfile = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using BMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class GammaSuffStats {
    double sum_logx{0};
    std::vector<double> sum_x;
    int sum_ind_total{0};
    std::vector<int> sum_ind;

    int Ncond{0}, Nsite{0}, Naa{20};
    std::vector<DMatrix>* fitness{nullptr};
    std::vector<BMatrix>* ind_conv{nullptr};

  public:
    GammaSuffStats() = default;

    GammaSuffStats(int Ncond, int Nsite, std::vector<DMatrix>* fitness,
                   std::vector<BMatrix>* ind_conv)
        : Ncond(Ncond), Nsite(Nsite), fitness(fitness), ind_conv(ind_conv) {}

    void collect() {
        sum_logx = 0;
        sum_ind_total = 0;
        sum_x.clear();
        sum_ind.clear();
        auto& x = *fitness;
        auto& ind = *ind_conv;
        for (int aa = 0; aa < Naa; aa++) {
            sum_x.push_back(0);
            sum_ind.push_back(0);
            for (int k = 0; k < Ncond; k++) {
                for (int i = 0; i < Nsite; i++) {
                    if (ind[k](i, aa)) {
                        sum_ind_total += 1;
                        sum_ind[aa] += 1;
                        sum_logx += log(x[k](i, aa));
                        sum_x[aa] += x[k](i, aa);
                    }
                }
            }
        }
    }

    double partial_density(double shape, AAProfile invrate) {
        double result =
            sum_ind_total * (shape * log(shape) - log(tgamma(shape))) + (shape - 1) * sum_logx;
        for (int aa = 0; aa < Naa; aa++) {
            result += -shape * (log(invrate[aa]) * sum_ind[aa] + sum_x[aa] / invrate[aa]);
        }
        return result;
    }

    double partial_density_invrate(double shape, AAProfile invrate) {
        double result = 0;
        for (int aa = 0; aa < Naa; aa++) {
            result += -shape * (log(invrate[aa]) * sum_ind[aa] + sum_x[aa] / invrate[aa]);
        }
        return result;
    }
};

/*
#### TESTS #########################################################################################
*/
TEST_CASE("Comparing gamma suff stat to logprob computed normally.") {
    auto InitUniformDirichlet = [](Eigen::VectorXd& v) {
        double tot = 0.;
        for (int i = 0; i < v.size(); i++) {
            v[i] = Random::sExpo();
            tot += v[i];
        }
        for (int i = 0; i < v.size(); i++) {
            v[i] /= tot;
        }
    };

    double fitness_shape;
    int Ncond{2}, Nsite{5}, Naa{20};
    AAProfile fitness_inv_rates;
    std::vector<DMatrix> fitness;

    fitness_shape = Random::sExpo();
    fitness_inv_rates = AAProfile(Naa);

    InitUniformDirichlet(fitness_inv_rates);

    fitness = std::vector<DMatrix>(Ncond, Eigen::MatrixXd(Nsite, Naa));
    for (auto& G_k : fitness) {
        for (int i = 0; i < Nsite; i++)
            for (int aa = 0; aa < Naa; aa++) {
                G_k(i, aa) = Random::Gamma(fitness_shape, fitness_shape / fitness_inv_rates[aa]);
            }
    }

    auto ind_conv = std::vector<BMatrix>(Ncond, BMatrix(Nsite, Naa));
    for (int k = 0; k < Ncond; k++)
        for (int i = 0; i < Nsite; i++)
            for (int aa = 0; aa < Naa; aa++) {
                ind_conv[k](i, aa) = (Random::Uniform() < 0.5);
            }

    GammaSuffStats gamma_suff_stats(Ncond, Nsite, &fitness, &ind_conv);

    auto partial_gamma_log_density = [](double alpha, double m, double x) {
        return -alpha * log(m) - alpha * x / m;
    };

    auto fitness_log_density = [&]() {
        double logprob = 0;
        for (int k = 0; k < Ncond; k++)
            for (int i = 0; i < Nsite; i++)
                for (int aa = 0; aa < Naa; aa++)
                    logprob += ind_conv[k](i, aa) * partial_gamma_log_density(fitness_shape,
                                                                              fitness_inv_rates[aa],
                                                                              fitness[k](i, aa));
        return logprob;
    };

    auto partial_gamma_log_density2 = [](double alpha, double m, double x) {
        double beta = alpha / m;
        return alpha * log(beta) - log(tgamma(alpha)) + (alpha - 1) * log(x) - beta * x;
    };

    auto fitness_log_density2 = [&]() {
        double logprob = 0;
        for (int k = 0; k < Ncond; k++)
            for (int i = 0; i < Nsite; i++)
                for (int aa = 0; aa < Naa; aa++)
                    logprob += ind_conv[k](i, aa) *
                               partial_gamma_log_density2(fitness_shape, fitness_inv_rates[aa],
                                                          fitness[k](i, aa));
        return logprob;
    };

    gamma_suff_stats.collect();
    CHECK(gamma_suff_stats.partial_density(fitness_shape, fitness_inv_rates) -
              fitness_log_density2() <
          10e-10);
    CHECK(gamma_suff_stats.partial_density_invrate(fitness_shape, fitness_inv_rates) -
              fitness_log_density() <
          10e-10);
}
