
#include "ufun.h"
double cal_energy(
    const Vec &phi_means,
    const Vec2D &chis,
    const Vec &Ls,
    const Vec &zs,
    const Vec &kappas,
    double v,
    Vec3D &omegas,
    Vec3D &phis,
    Vec2D &psi,
    Vec &Js,
    Vec &Lbetas,
    Vec &zetas,
    Vec &ps,
    Vec &local_energy);

// Main function implementing Gibbs dynamics
std::tuple<double, double, double, double, double, int>
Gibbs_dynamics(
    const Vec &phi_means,
    const Vec2D &chis,
    const Vec &Ls,
    const Vec &zs,
    const Vec &kappas,
    double v,
    Vec3D &omegas,
    Vec3D &phis,
    Vec2D &psi,
    int steps_inner,
    double acceptance_omega,
    double acceptance_J,
    double acceptance_Lbeta,
    double acceptance_zeta,
    Vec &Js,
    Vec &Lbetas,
    Vec &zetas,
    Vec &ps,
    bool flag_C = false,
    bool flag_zetas = false,
    bool flag_ps = false,
    double C = 1,
    double threshold_incomp = 1e-10,
    double threshold_omega = 1e-8,
    double threshold_J = 1e-6,
    double threshold_Lbeta = 1e-8,
    double threshold_zeta = 1e-7,
    bool flag_oneperiod = false)
{
    std::default_random_engine generator(0);

    double max_abs_incomp, max_omega_diff, max_J_diff, max_Lbeta_error, max_zeta_error;

    int num_comps = omegas.size();
    int num_beta = omegas[0].size();
    int num_coord = omegas[0][0].size();
    // print(num_comps);
    // print(num_beta);
    // print(num_coord);

    double chi_sum_sum = 0.0;
    for (int i = 0; i < num_comps; i++)
        for (int j = 0; j < num_comps; j++)
        {
            chi_sum_sum += chis[i][j];
        }
    // print("chi_sum_sum=");
    // print(chi_sum_sum);

    for (int istep = 0; istep < steps_inner; ++istep)
    {
        Vec2D k2s(num_beta, Vec(num_coord, 0.0));

        for (int i = 0; i < num_beta; ++i)
        {
            double samplingRate = 1.0 / (Lbetas[i] / num_coord);
            int numDataPoints = num_coord;
            std::vector<double> freq = calculateFrequencies(samplingRate, numDataPoints);

            for (int j = 0; j < num_coord; ++j)
            {
                k2s[i][j] = (2.0 * M_PI * freq[j]) * (2.0 * M_PI * freq[j]);
            }
        }
        // print("k2s=");
        // print(k2s);

        Vec2D inverse_k2s = k2s;
        for (int i = 0; i < num_beta; ++i)
        {
            inverse_k2s[i][0] = 0.0; // Avoid division by zero
            for (int j = 1; j < num_coord; ++j)
            {
                inverse_k2s[i][j] = 1.0 / k2s[i][j];
            }
        }
        // print("inverse_k2s=");
        // print(inverse_k2s);

        double rate_LPF = kappas[0] * 10;
        Vec2D LPFkernel(num_beta, Vec(num_coord, 0.0));
        for (int i = 0; i < num_beta; ++i)
        {
            for (int j = 0; j < num_coord; ++j)
            {
                LPFkernel[i][j] = 1 / (1 / acceptance_omega + rate_LPF * std::abs(k2s[i][j]));
            }
            LPFkernel[i][0] = 1 / (1 / acceptance_omega);
        }

        // print(LPFkernel);

        Vec Qs(num_comps, 0.0);
        // print(Qs);
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    Qs[itr_comp] += std::exp(-omegas[itr_comp][itr_beta][itr_coord] * Ls[itr_comp]) * Js[itr_beta];
                }
            }
            Qs[itr_comp] /= double(num_coord * num_beta);
        }
        // print("Qs=");
        // print(Qs);

        Vec2D incomp(num_beta, Vec(num_coord, 1.0));
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            double factor = phi_means[itr_comp] / Qs[itr_comp];
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    phis[itr_comp][itr_beta][itr_coord] = factor * std::exp(-omegas[itr_comp][itr_beta][itr_coord] * Ls[itr_comp]);
                    incomp[itr_beta][itr_coord] -= phis[itr_comp][itr_beta][itr_coord];
                }
            }
        }
        // print("phis=");
        // print(phis);
        // print("incomp=");
        // print(incomp);
        max_abs_incomp = 0;

        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                if (max_abs_incomp < std::abs(incomp[itr_beta][itr_coord]))
                {
                    max_abs_incomp = std::abs(incomp[itr_beta][itr_coord]);
                }
            }
        }
        // print(max_abs_incomp);

        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                psi[itr_beta][itr_coord] = 0.0;
            }
        }

        Vec2D rho(num_beta, Vec(num_coord, 0.0));
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
                {
                    rho[itr_beta][itr_coord] += zs[itr_comp] * phis[itr_comp][itr_beta][itr_coord];
                }
                rho[itr_beta][itr_coord] *= (4.0 * M_PI / v);
            }
            psi[itr_beta] = Convolution(rho[itr_beta], inverse_k2s[itr_beta]);
        }
        // print("psi=");
        // print(psi);

        Vec charge(num_beta, 0.0);
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    charge[itr_beta] += phis[itr_comp][itr_beta][itr_coord] * zs[itr_comp];
                }
            }
            charge[itr_beta] /= (num_coord);
        }
        // print("charge=");
        // print(charge);

        Vec3D chiphis = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {

                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    chiphis[itr_comp][itr_beta][itr_coord] = 0.0;
                    for (int i = 0; i < num_comps; i++)
                    {
                        chiphis[itr_comp][itr_beta][itr_coord] += chis[itr_comp][i] * phis[i][itr_beta][itr_coord];
                    }
                }
            }
        }
        // print("chiphis=");
        // print(chiphis);

        // kappaterm=-kappa_i*nabla^2 phi_i(beta,r)
        Vec3D kappaterm = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                kappaterm[itr_comp][itr_beta] = Convolution(phis[itr_comp][itr_beta], k2s[itr_beta]);
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    kappaterm[itr_comp][itr_beta][itr_coord] *= kappas[itr_comp];
                }
            }
        }
        // print("kappaterm=");
        // print(kappaterm);

        Vec3D omega_temp = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    omega_temp[itr_comp][itr_beta][itr_coord] =
                        chiphis[itr_comp][itr_beta][itr_coord] + kappaterm[itr_comp][itr_beta][itr_coord] + zs[itr_comp] * psi[itr_beta][itr_coord];
                    if (flag_zetas)
                    {
                        omega_temp[itr_comp][itr_beta][itr_coord] += zs[itr_comp] * zetas[itr_beta];
                    }
                    if (flag_C)
                    {
                        omega_temp[itr_comp][itr_beta][itr_coord] += 2 * C * zs[itr_comp] * charge[itr_beta];
                    }
                }
            }
        }
        // print("omega_temp=");
        // print(omega_temp);

        // laplace_incomp=-nabla^2 incomp(beta,r)
        Vec2D laplace_incomp = std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0));
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            laplace_incomp[itr_beta] = Convolution(incomp[itr_beta], k2s[itr_beta]);
        }

        // print("laplace_incomp=");
        // print(laplace_incomp);

        Vec3D omega_temp_shift = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    omega_temp_shift[itr_comp][itr_beta][itr_coord] = omega_temp[itr_comp][itr_beta][itr_coord] + kappas[itr_comp] * laplace_incomp[itr_beta][itr_coord];
                }
            }
        }

        // print("omega_temp_shift=");
        // print(omega_temp_shift);

        Vec2D xi(num_beta, Vec(num_coord, 0.0));
        // xi is -xi in the note
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                xi[itr_beta][itr_coord] += chi_sum_sum * incomp[itr_beta][itr_coord];
                for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
                {
                    xi[itr_beta][itr_coord] += omega_temp_shift[itr_comp][itr_beta][itr_coord] - omegas[itr_comp][itr_beta][itr_coord];
                }
                xi[itr_beta][itr_coord] /= num_comps;
            }
        }
        // print("xi=");
        // print(xi);

        Vec3D omega_new = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    omega_new[itr_comp][itr_beta][itr_coord] = omega_temp[itr_comp][itr_beta][itr_coord] - xi[itr_beta][itr_coord];
                }
            }
        }
        // print("omega_new=");
        // print(omega_new);

        // laplace_psi=-nabla^2 psi(beta,r)
        Vec2D laplace_psi(num_beta, Vec(num_coord, 0.0));
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            laplace_psi[itr_beta] = Convolution(psi[itr_beta], k2s[itr_beta]);
        }

        Vec local_energy(num_beta, 0.0); // this is -eta(beta)in the note
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
                {
                    local_energy[itr_beta] +=
                        (-phis[itr_comp][itr_beta][itr_coord] * omega_new[itr_comp][itr_beta][itr_coord] + 0.5 * phis[itr_comp][itr_beta][itr_coord] * chiphis[itr_comp][itr_beta][itr_coord] + 0.5 * phis[itr_comp][itr_beta][itr_coord] * kappaterm[itr_comp][itr_beta][itr_coord] + psi[itr_beta][itr_coord] * zs[itr_comp] * phis[itr_comp][itr_beta][itr_coord] - phis[itr_comp][itr_beta][itr_coord] / Ls[itr_comp]);

                    // if (flag_zetas)
                    // {
                    //     local_energy[itr_beta] += zetas[itr_beta] * zs[itr_comp] * phis[itr_comp][itr_beta][itr_coord];
                    // }
                }

                local_energy[itr_beta] += xi[itr_beta][itr_coord] * incomp[itr_beta][itr_coord] - v / 8.0 / M_PI * psi[itr_beta][itr_coord] * laplace_psi[itr_beta][itr_coord];
            }
            local_energy[itr_beta] /= num_coord;
            if (flag_zetas)
            {
                local_energy[itr_beta] += zetas[itr_beta] * charge[itr_beta];
            }
            if (flag_C)
            {
                local_energy[itr_beta] += C * charge[itr_beta] * charge[itr_beta];
            }
        }
        // print("local_energy");
        // print(local_energy);

        double local_energy_mean = 0.0;
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            local_energy_mean += local_energy[itr_beta] * Js[itr_beta];
        }
        local_energy_mean /= num_beta;
        // print("local_energy_mean");
        // print(local_energy_mean);

        Vec Jdiff(num_beta, 0.0);
        max_J_diff = 0.0;
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            Jdiff[itr_beta] = local_energy_mean - local_energy[itr_beta];
            max_J_diff = std::max(max_J_diff, std::abs(Jdiff[itr_beta]));
        }

        // print("Jdiff");
        // print(Jdiff);
        // print("max_J_diff");
        // print(max_J_diff);

        double Js_step_upperbound = 0.001;
        double Js_max_change = max_J_diff * acceptance_J;
        double additional_factor = Js_step_upperbound / std::max(Js_max_change, Js_step_upperbound);
        double Jsmean = 0.0;
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            Js[itr_beta] += additional_factor * acceptance_J * Jdiff[itr_beta];
            Jsmean += Js[itr_beta];
        }
        Jsmean /= num_beta;
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            Js[itr_beta] += (1.0 - Jsmean);
        }
        // print("Js=");
        // print(Js);

        int flag_bad_J = 0;
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            if (Js[itr_beta] < 1e-3) // replace the compartment by its neighbor
            {
                flag_bad_J = 1;

                // double amplitude = 5.0;
                // std::normal_distribution<double> distribution(0.0, amplitude);
                // for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
                // {
                //     for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                //     {
                //         omegas[itr_comp][itr_beta][itr_coord] = distribution(generator);
                //     }
                // }

                int itr_beta_copied = itr_beta - 1;
                if (itr_beta_copied < 0)
                {
                    itr_beta_copied += num_beta;
                }

                for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
                {
                    for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                    {
                        omegas[itr_comp][itr_beta][itr_coord] = omegas[itr_comp][itr_beta_copied][itr_coord];
                    }
                }

                Js[itr_beta] = Js[itr_beta_copied];
                Lbetas[itr_beta] = Lbetas[itr_beta_copied];
                zetas[itr_beta] = zetas[itr_beta_copied];
            }
        }
        if (flag_bad_J > 0)
        {
            print("reproduce a J(beta) and omega");
            // run with 0 acceptances
            std::tuple<double, double, double, double, double, int> errorsinside1 = Gibbs_dynamics(
                phi_means,
                chis,
                Ls,
                zs,
                kappas,
                v,
                omegas,
                phis,
                psi,
                30000,
                acceptance_omega,
                0.,
                0.,
                0.,
                Js,
                Lbetas,
                zetas,
                ps,
                flag_C,
                flag_zetas,
                flag_ps,
                C,
                threshold_incomp,
                threshold_omega,
                threshold_J,
                threshold_Lbeta,
                threshold_zeta);

            std::tuple<double, double, double, double, double, int> errorsinside2 = Gibbs_dynamics(
                phi_means,
                chis,
                Ls,
                zs,
                kappas,
                v,
                omegas,
                phis,
                psi,
                steps_inner - istep,
                acceptance_omega,
                acceptance_J,
                acceptance_Lbeta,
                acceptance_zeta,
                Js,
                Lbetas,
                zetas,
                ps,
                flag_C,
                flag_zetas,
                flag_ps,
                C,
                threshold_incomp,
                threshold_omega,
                threshold_J,
                threshold_Lbeta,
                threshold_zeta,
                flag_oneperiod);

            return errorsinside2;

            // break;
        }

        max_omega_diff = 0.0;
        Vec3D Domega = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
        Vec Domegamean(num_comps, 0.0);

        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    Domega[itr_comp][itr_beta][itr_coord] = omega_new[itr_comp][itr_beta][itr_coord] - omegas[itr_comp][itr_beta][itr_coord];
                    Domegamean[itr_comp] += Domega[itr_comp][itr_beta][itr_coord]; // why  not time Js[itr_beta] here? try this later.
                }
            }
            Domegamean[itr_comp] /= (num_beta * num_coord);
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    Domega[itr_comp][itr_beta][itr_coord] -= Domegamean[itr_comp];
                    max_omega_diff = std::max(std::abs(Domega[itr_comp][itr_beta][itr_coord]), max_omega_diff);
                }
            }
        }
        // print("Domega=");
        // print(Domega);
        // print("Domegamean");
        // print(Domegamean);
        // print("max_omega_diff=");
        // print(max_omega_diff);

        // calculate Domega_real using LPFkernel
        Vec3D Domega_real = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
        Vec Domegamean_real(num_comps, 0.0);
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                Domega_real[itr_comp][itr_beta] = Convolution(Domega[itr_comp][itr_beta], LPFkernel[itr_beta]);
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    Domegamean_real[itr_comp] += Domega_real[itr_comp][itr_beta][itr_coord];
                }
            }
            Domegamean_real[itr_comp] /= (num_beta * num_coord);
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    Domega_real[itr_comp][itr_beta][itr_coord] -= Domegamean_real[itr_comp];
                }
            }
        }
        // print("Domega_real=");
        // print(Domega_real);
        // print("Domegamean_real");
        // print(Domegamean_real);

        // update omega
        Vec omegamean(num_comps, 0.0);
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    omegas[itr_comp][itr_beta][itr_coord] += additional_factor * Domega_real[itr_comp][itr_beta][itr_coord];
                    omegamean[itr_comp] += omegas[itr_comp][itr_beta][itr_coord];
                }
            }
            omegamean[itr_comp] /= (num_beta * num_coord);

            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    omegas[itr_comp][itr_beta][itr_coord] -= omegamean[itr_comp];
                }
            }
        }
        // print("omegas=");
        // print(omegas);

        // Here we apply omega(x)=omega(x)+omega(L-x))/2
        Vec3D omegas_sym = omegas;
        omegamean = Vec(num_comps, 0.0);
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    omegas_sym[itr_comp][itr_beta][itr_coord] = (omegas[itr_comp][itr_beta][itr_coord] + omegas[itr_comp][itr_beta][num_coord - itr_coord - 1]) * 0.5;
                    omegamean[itr_comp] += omegas_sym[itr_comp][itr_beta][itr_coord];
                }
            }
            omegamean[itr_comp] /= (num_beta * num_coord);
        }
        for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    omegas[itr_comp][itr_beta][itr_coord] = omegas_sym[itr_comp][itr_beta][itr_coord] - omegamean[itr_comp];
                }
            }
        }

        if (flag_ps)
        {
            Vec2D qs(num_comps, Vec(num_beta, 0.0));
            Qs = Vec(num_comps, 0.0);
            // print(Qs);
            for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
            {
                for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
                {
                    for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                    {
                        qs[itr_comp][itr_beta] += std::exp(-omegas[itr_comp][itr_beta][itr_coord] * Ls[itr_comp]);
                    }
                    qs[itr_comp][itr_beta] /= num_coord;
                    Qs[itr_comp] += qs[itr_comp][itr_beta] * Js[itr_beta];
                }
                Qs[itr_comp] /= double(num_beta);
            }
            // print("qs=");
            // print(qs);
            // print("Qs=");
            // print(Qs);

            Vec2D hatphi(num_comps, Vec(num_beta, 0.0));
            for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
            {
                for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
                {
                    hatphi[itr_comp][itr_beta] = phi_means[itr_comp] * qs[itr_comp][itr_beta] / Qs[itr_comp];
                }
            }
            Vec rhosum(num_beta, 0.0);
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
                {
                    rhosum[itr_beta] += zs[itr_comp] * hatphi[itr_comp][itr_beta];
                }
            }

            // print("rhosum");
            // print(rhosum);

            double pzsum = 0.0;
            for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
            {
                pzsum += std::abs(zs[itr_comp]) * ps[itr_comp];
            }
            // print("pzsum=");
            // print(pzsum);

            Vec2D domega(num_comps, Vec(num_beta, 0.0));
            for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
            {
                for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
                {
                    double tmp = 1.0 - rhosum[itr_beta] / hatphi[itr_comp][itr_beta] * ps[itr_comp] * zs[itr_comp] / std::abs(zs[itr_comp]) / pzsum;

                    if (std::abs(zs[itr_comp]) > 1e-10)
                    {
                        if (tmp > 0)
                        {
                            domega[itr_comp][itr_beta] = std::log(tmp) / (-Ls[itr_comp]);
                        }
                        else
                        {
                        }
                    }
                }
            }
            // print("domega=");
            // print(domega);

            omegamean = Vec(num_comps, 0.0);
            for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
            {
                for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
                {
                    for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                    {
                        omegas[itr_comp][itr_beta][itr_coord] += domega[itr_comp][itr_beta];
                        omegamean[itr_comp] += omegas[itr_comp][itr_beta][itr_coord];
                    }
                }
                omegamean[itr_comp] /= (num_beta * num_coord);
                for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
                {
                    for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                    {
                        omegas[itr_comp][itr_beta][itr_coord] -= omegamean[itr_comp];
                    }
                }
            }
            // print("omegas=");
            // print(omegas);
        }

        if (flag_zetas)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                zetas[itr_beta] += additional_factor * acceptance_zeta * charge[itr_beta] * Js[itr_beta];
            }
        }
        // print("zetas");
        // print(zetas);

        // note that we may need the updated phis, psi, etc.. to calculate rate_Lbeta
        Vec rate_Lbeta(num_beta, 0.0);
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
            {
                for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
                {
                    rate_Lbeta[itr_beta] += 0.5 * phis[itr_comp][itr_beta][itr_coord] * kappaterm[itr_comp][itr_beta][itr_coord];
                }
            }
            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                rate_Lbeta[itr_beta] += -v / 8.0 / M_PI * psi[itr_beta][itr_coord] * laplace_psi[itr_beta][itr_coord];
            }
            rate_Lbeta[itr_beta] /= (num_coord);
            rate_Lbeta[itr_beta] *= Js[itr_beta] * (-2.0 / Lbetas[itr_beta]);
        }
        // print("rate_Lbeta");
        // print(rate_Lbeta);
        max_Lbeta_error = 0.0;
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            Lbetas[itr_beta] -= additional_factor * acceptance_Lbeta * rate_Lbeta[itr_beta];
            max_Lbeta_error = std::max(max_Lbeta_error, std::abs(acceptance_Lbeta * rate_Lbeta[itr_beta] / Lbetas[itr_beta]));
        }
        // print("max_Lbeta_error");
        // print(max_Lbeta_error);
        // print("additional factor");
        // print(additional_factor);
        // print("Lbetas");
        // print(Lbetas);

        max_zeta_error = 0.0;
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            max_zeta_error = std::max(max_zeta_error, std::abs(charge[itr_beta]));
        }
        // print("max_zeta_error");
        // print(max_zeta_error);

        //     double max_Lbeta_error = *std::max_element(local_width.begin(), local_width.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });

        if ((istep + 1) % 10000 == 0)
        {
            std::cout << std::setprecision(8);
            std::cout << "Step: " << std::setw(5) << istep
                      << ",    max_abs_incomp: " << std::setw(13) << max_abs_incomp
                      << ",    max_omega_diff: " << std::setw(13) << max_omega_diff
                      << ",    max_J_diff: " << std::setw(13) << max_J_diff
                      << ",    max_zeta_error: " << std::setw(13) << max_zeta_error
                      << ",    max_Lbeta_error: " << std::setw(13) << max_Lbeta_error
                      << ",    charges:";
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                std::cout << std::setw(13) << charge[itr_beta] << ' ';
            }
            std::cout << ",    Lbetas: ";
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                std::cout << std::setw(13) << Lbetas[itr_beta] << ' ';
            }
            std::cout << ",    Js: ";
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                std::cout << std::setw(13) << Js[itr_beta] << ' ';
            }
            std::cout << ",    zetas: ";
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                std::cout << std::setw(13) << zetas[itr_beta] << ' ';
            }
            std::cout << ",    additional_factor: ";
            std::cout << std::setw(13) << additional_factor;
            // print("pre_local_energy");
            // print(local_energy);
            double total_energy = cal_energy(
                phi_means,
                chis,
                Ls,
                zs,
                kappas,
                v,
                omegas,
                phis,
                psi,
                Js,
                Lbetas,
                zetas,
                ps,
                local_energy);
            print("local_energy");
            print(local_energy);
            std::cout << ",    total_energy: ";
            std::cout << std::setw(13) << total_energy;
            std::cout << std::endl;
        }

        ////// if there are more than one peaks, we rescale to one period.
        if ((istep + 1) % 10000 == 0 && istep > 10000 && flag_oneperiod)
        {
            for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
            {
                // print("size of omega [itr_beta][num_comps - 1]");
                // int size = (omegas[num_comps - 1][itr_beta]).size();
                // print(size);
                int period = find_period(omegas[num_comps - 1][itr_beta], num_coord);
                // print("period=");
                // print(period);

                if (period > 1)
                {
                    print("period=");
                    print(period);
                    print("rescale the profile");
                    for (int icomp = 0; icomp < num_comps; icomp++)
                    {
                        rescale(omegas[icomp][itr_beta], period, num_coord);
                    }
                    Lbetas[itr_beta] /= period;
                }
            }
        }
        else
        {
            if (max_abs_incomp < threshold_incomp && max_omega_diff < threshold_omega && max_J_diff < threshold_J && max_Lbeta_error < threshold_Lbeta && max_zeta_error < threshold_zeta)
            {
                return std::make_tuple(max_abs_incomp, max_omega_diff, max_J_diff, max_Lbeta_error, max_zeta_error, 1); // the last value is 1 if converged, else is 0
            }
        }
    }
    return std::make_tuple(max_abs_incomp, max_omega_diff, max_J_diff, max_Lbeta_error, max_zeta_error, 0);
}

double cal_energy(
    const Vec &phi_means,
    const Vec2D &chis,
    const Vec &Ls,
    const Vec &zs,
    const Vec &kappas,
    double v,
    Vec3D &omegas,
    Vec3D &phis,
    Vec2D &psi,
    Vec &Js,
    Vec &Lbetas,
    Vec &zetas,
    Vec &ps,
    Vec &local_energy)
{
    int num_comps = omegas.size();
    int num_beta = omegas[0].size();
    int num_coord = omegas[0][0].size();
    double chi_sum_sum = 0.0;
    for (int i = 0; i < num_comps; i++)
        for (int j = 0; j < num_comps; j++)
        {
            chi_sum_sum += chis[i][j];
        }

    Vec2D k2s(num_beta, Vec(num_coord, 0.0));

    for (int i = 0; i < num_beta; ++i)
    {
        double samplingRate = 1.0 / (Lbetas[i] / num_coord);
        int numDataPoints = num_coord;
        std::vector<double> freq = calculateFrequencies(samplingRate, numDataPoints);

        for (int j = 0; j < num_coord; ++j)
        {
            k2s[i][j] = (2.0 * M_PI * freq[j]) * (2.0 * M_PI * freq[j]);
        }
    }
    // print("k2s=");
    // print(k2s);

    Vec2D inverse_k2s = k2s;
    for (int i = 0; i < num_beta; ++i)
    {
        inverse_k2s[i][0] = 0.0; // Avoid division by zero
        for (int j = 1; j < num_coord; ++j)
        {
            inverse_k2s[i][j] = 1.0 / k2s[i][j];
        }
    }

    Vec2D incomp(num_beta, Vec(num_coord, 1.0));
    for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
    {

        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                incomp[itr_beta][itr_coord] -= phis[itr_comp][itr_beta][itr_coord];
            }
        }
    }

    double total_energy = 0.0;

    Vec3D chiphis = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
    for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
    {
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {

            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                chiphis[itr_comp][itr_beta][itr_coord] = 0.0;
                for (int i = 0; i < num_comps; i++)
                {
                    chiphis[itr_comp][itr_beta][itr_coord] += chis[itr_comp][i] * phis[i][itr_beta][itr_coord];
                }
            }
        }
    }

    // kappaterm=-kappa_i*nabla^2 phi_i(beta,r)
    Vec3D kappaterm = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
    for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
    {
        for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
        {
            kappaterm[itr_comp][itr_beta] = Convolution(phis[itr_comp][itr_beta], k2s[itr_beta]);
            for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
            {
                kappaterm[itr_comp][itr_beta][itr_coord] *= kappas[itr_comp];
            }
        }
    }

    // laplace_psi=-nabla^2 psi(beta,r)
    Vec2D laplace_psi(num_beta, Vec(num_coord, 0.0));
    for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
    {
        laplace_psi[itr_beta] = Convolution(psi[itr_beta], k2s[itr_beta]);
    }

    local_energy = Vec(num_beta, 0.0); // this is -eta(beta)in the note
    for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
    {
        for (int itr_coord = 0; itr_coord < num_coord; ++itr_coord)
        {
            for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
            {
                local_energy[itr_beta] +=
                    phis[itr_comp][itr_beta][itr_coord] * std::log(phis[itr_comp][itr_beta][itr_coord]) / Ls[itr_comp] + 0.5 * phis[itr_comp][itr_beta][itr_coord] * chiphis[itr_comp][itr_beta][itr_coord] + 0.5 * phis[itr_comp][itr_beta][itr_coord] * kappaterm[itr_comp][itr_beta][itr_coord] + psi[itr_beta][itr_coord] * zs[itr_comp] * phis[itr_comp][itr_beta][itr_coord];
            }

            local_energy[itr_beta] += -v / 8.0 / M_PI * psi[itr_beta][itr_coord] * laplace_psi[itr_beta][itr_coord];
        }
        local_energy[itr_beta] /= num_coord;
    }
    // print("local_energy");
    // print(local_energy);

    double local_energy_mean = 0.0;
    for (int itr_beta = 0; itr_beta < num_beta; ++itr_beta)
    {
        local_energy_mean += local_energy[itr_beta] * Js[itr_beta];
    }
    local_energy_mean /= num_beta;

    total_energy = local_energy_mean;

    for (int itr_comp = 0; itr_comp < num_comps; ++itr_comp)
    {
        total_energy -= phi_means[itr_comp] * std::log(phi_means[itr_comp]) / Ls[itr_comp];
    }

    return total_energy;
}
