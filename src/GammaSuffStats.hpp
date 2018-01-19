#include <vector>
#include "Eigen/Dense"
#include "Random.hpp"
#include "doctest.h"
using DMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using AAProfile = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using BMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class GammaSuffStats {
    double sum_x{0}, sum_logx{0};
    int Ncond{0}, Nsite{0}, Naa{20};
    std::vector<DMatrix>* fitness{nullptr};

  public:
    GammaSuffStats() = default;

    GammaSuffStats(int Ncond, int Nsite, std::vector<DMatrix>* fitness)
        : Ncond(Ncond), Nsite(Nsite), fitness(fitness) {}

    void gather() {
        sum_logx = 0;
        sum_x = 0;
        auto& x = *fitness;
        for (int k = 0; k < Ncond; k++) {
            for (int i = 0; i < Nsite; i++) {
                for (int aa = 0; aa < Naa; aa++) {
                    sum_logx += log(x[k](i, aa));
                    sum_x += x[k](i, aa);
                }
            }
        }
    }

    double partial_density_shape(double shape, AAProfile invshape) {
        double result = Naa * Nsite * Ncond * (shape * log(shape) - log(tgamma(shape))) +
                        (shape - 1) * sum_logx;
        for (int aa = 0; aa < Naa; aa++) {
            result += -Ncond * Nsite * shape * (log(invshape[aa]) + sum_x / invshape[aa]);
        }
        return result;
    }

    double partial_density_invshape(double shape, AAProfile invshape) {
        double result = 0;  // this term was not dependent on invshape
        for (int aa = 0; aa < Naa; aa++) {
            result += -Ncond * Nsite * shape * (log(invshape[aa]) + sum_x / invshape[aa]);
        }
        return result;
    }
};

/*
  #### TESTS
  #########################################################################################
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

    GammaSuffStats gamma_suff_stats(Ncond, Nsite, &fitness);

    auto partial_gamma_log_density = [](double alpha, double m, double x) {
        double beta = alpha / m;
        return alpha * log(beta) - beta * x;
    };

    auto fitness_log_density = [&]() {
        double logprob =
            gamma_suff_stats.partial_density_invshape(fitness_shape, fitness_inv_rates);
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
        double logprob = -fitness_shape;
        for (int k = 0; k < Ncond; k++)
            for (int i = 0; i < Nsite; i++)
                for (int aa = 0; aa < Naa; aa++)
                    logprob += ind_conv[k](i, aa) * partial_gamma_log_density2(fitness_shape,
                                                                              fitness_inv_rates[aa],
                                                                              fitness[k](i, aa));
        return logprob;
    };

    gamma_suff_stats.gather();
    CHECK(gamma_suff_stats.partial_density_shape(fitness_shape, fitness_inv_rates) ==
          fitness_log_density2());
    CHECK(gamma_suff_stats.partial_density_invshape(fitness_shape, fitness_inv_rates) ==
          fitness_log_density());
}
