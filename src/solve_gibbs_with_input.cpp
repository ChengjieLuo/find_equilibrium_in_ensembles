#include "gibbs_dynamics.h"

// Example usage of the Convolution function
int main(int argc, char *argv[])
{

    if (argc < 3)
    {
        print("./solve_gibbs_with_input.out {num_comps} {num_beta} {num_coord} {phi_means_str} {chis_str} {Ls_str} {zs_str} {kappas_str} {v} {omega_type} {omega_value} {steps_inner1} {steps_inner2} {acceptance_omega} {acceptance_J} {acceptance_Lbeta} {acceptance_zeta} {Js_type} {Js_value} {Lbeta_type} {Lbeta_value} {zetas_type} {zetas_value} {flag_C} {C} {ps_type} {ps_value} {flag_zetas} {flag_ps} {savedata_pre} {flag_save_separate} {threshold_incomp} {threshold_omega} {threshold_J} {threshold_Lbeta} {threshold_zeta}");
        print("##############Inputs##############");
        print("num_comps: 1 int, number of components. eg. 5");
        print("num_beta: 1 int, number of compartments. eg. 2");
        print("num_coord: 1 int, number of coordinates. eg. 128");
        print("phi_mean_str: num_comps double, mean value of volume fractions for num_comps species, the total value should be 1.0. eg. 0.1 0.1 0.1 0.1 0.6");
        print("chis_str: num_comps*num_comps double, interaction matrix.  eg. \n \
        0.0 -5.0 0.0 0.0 0.0 \n \
        -5.0 0.0 0.0 0.0 0.0\n \
        0.0 0.0 0.0 0.0 0.0\n \
        0.0 0.0 0.0 0.0 0.0\n \
        0.0 0.0 0.0 0.0 0.0");
        print("Ls_str: num_comps double, numbers of monomers. eg. 10 10 1 1 1");
        print("zs_str: num_comps double, charge densities. eg. 0.8 ,-0.2, -1.0 ,1.0, 0.0");
        print("kappas_str: num_comps double, interfacial parameters, which is square of interfacial width. eg. 10 10 10 10 0");
        print("v: double, volume of each monomer. eg. 100.0");
        print("omega_type: one string, type of omega input, 'file' or 'random' or 'random_coord_independent'. eg. random");
        print("omega_value: one string or one double,\n \t if omega_type is 'file', then the omega_value is the filename which includes omegas, there should be num_comps*num_beta*num_coord values in the file;\n \t if type is 'random', then the omega_value is the amplitude of a Gaussian random distribution with mean 0, all the omega values are randomly produced;\n \t  if type is 'random_coord_independent', then the omega_value is the amplitude of a Gaussian random distribution with mean 0, but omega(i,beta,r) are the same for all position r.\n \t eg. 5.0");
        print("steps_inner1: one int, number of steps in the first loop. In this loop, to be stable, only omega will be updated. eg. 10001");
        print("steps_inner2: one int, number of steps in the second loop. In this loop, all fields will be updated. eg. 100001");
        print("acceptance_omega: one double, acceptance rate of omegas. eg. 0.001");
        print("acceptance_J:one double, acceptance rate of Js. eg. 0.001");
        print("acceptance_Lbeta:one double, acceptance rate of Lbeta. eg. 10");
        print("acceptance_zeta:one double, acceptance rate of zeta. eg. 10");

        print("Js_type: one string, type of Js input, 'file' or 'random'. eg. random");
        print("Js_value: one string or one double,\n \t if Js_type is 'file', then the Js_value is the filename which includes Js, there should be num_beta values in the file;\n \t if type is 'random', then the omega_value is the amplitude of a Gaussian random distribution with mean 1.0, all the Js values are randomly produced;\n \t eg. 0.0");

        print("Lbeta_type: one string, type of Lbeta input, 'file' or 'equal'. eg. equal");
        print("Lbeta_value: one string or one double,\n \t if Lbeta_type is 'file', then the Lbeta_value is the filename which includes Lbetas, there should be num_beta values in the file;\n \t if type is 'equal', then the Lbeta_value is the Lbeta value for all compartments;\n \t eg. 30.0");
        

        print("zetas_type: one string, type of zetas input, 'file' or 'equal'. eg. equal");
        print("zetas_value: one string or one double,\n \t if zetas_type is 'file', then the zetas_value is the filename which includes zetas, there should be num_beta values in the file;\n \t if type is 'equal', then the zetas_value is the zeta value for all compartments;\n \t eg. 0.0");

        print("flag_C: one string. flag indicating whether to use parameter C to constrain the total charges, 'true' or 'false'. eg. false ");
        print("C: one double, value of the parameter C. eg. 100");

        print("ps_type: one string, type of ps input, 'file' or 'equal'. eg. equal");
        print("ps_value: one string or one double,\n \t if ps_type is 'file', then the ps_value is the filename which includes ps, there should be num_beta values in the file;\n \t if type is 'equal', then the ps_value is the ps value for all compartments;\n \t eg. 0.0");

        print("flag_zetas: one string. flag indicating whether updating zetas, 'true' or 'false'. eg. false ");

        print("flag_ps: one string. flag indicating whether updating ps, 'true' or 'false'. eg. false ");

        print("savedata_pre: one string, pre of the file to be saved, the data will be saved in file savedata_pre_alldata.txt. The data structure is \n\
        Lbetas(num_beta*double)  Js(num_beta*double) total_energy(1*double) errors(5*double, [max_abs_incomp, max_omega_diff, max_J_diff, max_Lbeta_error, max_zeta_error])\n \
        phis(num_comps*num_beta*num_coord: num_comps*num_beta rows and each row includes num_coord values)\n \
        omegas(num_comps*num_beta*num_coord: num_comps*num_beta rows and each row includes num_coord values)\n \
        psi(num_beta*num_coord: num_beta rows and each row includes num_coord values)\n \
         eg. ./data/result");

        print("flag_save_separate: one string, flag indicating whether to save the final data in different files, 'true' or 'false'. If true, _Lbetas.txt, _Js.txt, _phis.txt, _omegas.txt, _psi.txt, _zetas.txt. eg. true");


        print("########### the following parameters can be omitted#########");
        print("threshold_incomp: one double, threshold of max_abs_incomp, eg. 1e-10");
        print("threshold_omega: one double, threshold of max_omega_diff, eg. 1e-8");
        print("threshold_J: one double, threshold of max_J_diff, eg. 1e-6");
        print("threshold_Lbeta: one double, threshold of max_Lbeta_error, eg. 1e-8");
        print("threshold_zeta: one double, threshold of max_zeta_error, eg. 1e-7");

        print("---------- not successful ---------");
        print("PLEASE use the format mentioned above to successfully run the simulation.");

        return 0;
    }

    int arg_index = 1;

    int num_comps = std::stoi(argv[arg_index++]);
    int num_beta = std::stoi(argv[arg_index++]);
    int num_coord = std::stoi(argv[arg_index++]);

    Vec phi_means(num_comps, 0.0);
    for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
    {
        phi_means[itr_comp] = std::stod(argv[arg_index++]);
    }

    Vec2D chis(num_comps, Vec(num_comps, 0.0));
    // double chi = std::stod(argv[arg_index++]);
    // chis[0][1] = chi;
    // chis[1][0] = chi;
    // if (chis[0][1] < 0)
    // {
    //     for (int i = 0; i < num_comps; i++)
    //     {
    //         for (int j = 0; j < num_comps; j++)
    //         {
    //             chis[i][j] -= chi;
    //         }
    //     }
    // }

    for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
    {
        for (int itr_comp2 = 0; itr_comp2 < num_comps; itr_comp2++)
        {
            chis[itr_comp][itr_comp2] = std::stod(argv[arg_index++]);
        }
    }

    Vec Ls(num_comps, 0.0);
    for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
    {
        Ls[itr_comp] = std::stod(argv[arg_index++]);
    }

    Vec zs(num_comps, 0.0);
    for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
    {
        zs[itr_comp] = std::stod(argv[arg_index++]);
    }

    Vec kappas(num_comps, 0.0);
    for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
    {
        kappas[itr_comp] = std::stod(argv[arg_index++]);
    }

    double v = std::stod(argv[arg_index++]);

    Vec3D omegas = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));

    std::string flag_omega_type = argv[arg_index++];
    if (flag_omega_type == "file")
    {
        std::string filename = argv[arg_index++];
        std::ifstream ifile(filename, std::ios::in);
        if (!ifile.is_open())
        {
            std::cerr << "There was a problem opening the input file!\n";
            exit(1); // exit or do additional error checking
        }
        for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                {
                    ifile >> omegas[itr_comp][itr_beta][itr_coord];
                }
            }
        }
        ifile.close();
    }
    else if (flag_omega_type == "random")
    {
        double amplitude = std::stod(argv[arg_index++]);
        std::default_random_engine generator(0);
        std::normal_distribution<double> distribution(0.0, amplitude);
        for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                {
                    omegas[itr_comp][itr_beta][itr_coord] = distribution(generator);
                }
            }
        }
    }
    else if (flag_omega_type == "random_coord_independent")
    {
        double amplitude = std::stod(argv[arg_index++]);
        std::default_random_engine generator(0);
        std::normal_distribution<double> distribution(0.0, amplitude);
        for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                double tmp = distribution(generator);
                for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                {
                    omegas[itr_comp][itr_beta][itr_coord] = tmp;
                }
            }
        }
    }

    Vec3D phis = std::vector<std::vector<std::vector<double>>>(num_comps, std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0)));
    // print("phis=");
    // print(phis);

    Vec2D psi = std::vector<std::vector<double>>(num_beta, std::vector<double>(num_coord, 0.0));
    // print("psi=");
    // print(psi);

    int steps_inner1 = std::stoi(argv[arg_index++]);
    int steps_inner2 = std::stoi(argv[arg_index++]);

    double acceptance_omega = std::stod(argv[arg_index++]);
    double acceptance_J = std::stod(argv[arg_index++]);
    double acceptance_Lbeta = std::stod(argv[arg_index++]);
    double acceptance_zeta = std::stod(argv[arg_index++]);

    Vec Js = std::vector<double>(num_beta, 1.0);
    std::string flag_J_type = argv[arg_index++];
    if (flag_J_type == "file")
    {
        std::string filename = argv[arg_index++];
        // print(filename);
        std::ifstream ifile(filename, std::ios::in);
        if (!ifile.is_open())
        {
            std::cerr << "There was a problem opening the input file!\n";
            exit(1); // exit or do additional error checking
        }
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            ifile >> Js[itr_beta];
        }

        ifile.close();
    }
    else if (flag_J_type == "random")
    {
        double amplitude = std::stod(argv[arg_index++]);
        std::default_random_engine generator(5);
        std::normal_distribution<double> distribution(0.0, amplitude);
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            Js[itr_beta] += distribution(generator);
        }
    }
    else
    {
        print("error in initialize Js.");
        return -1;
    }

    Vec Lbetas = std::vector<double>(num_beta, 10.0);
    std::string flag_L_type = argv[arg_index++];
    if (flag_L_type == "file")
    {
        std::string filename = argv[arg_index++];
        std::ifstream ifile(filename, std::ios::in);
        if (!ifile.is_open())
        {
            std::cerr << "There was a problem opening the input file!\n";
            exit(1); // exit or do additional error checking
        }
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            ifile >> Lbetas[itr_beta];
        }
        ifile.close();
    }
    else if (flag_L_type == "equal")
    {
        double valueL = std::stod(argv[arg_index++]);
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {

            Lbetas[itr_beta] = valueL;
        }
    }
    else
    {
        print("error in initialize Lbetas.");
        return -1;
    }

    Vec zetas = std::vector<double>(num_beta, 0.0);
    std::string flag_zeta_type = argv[arg_index++];
    if (flag_zeta_type == "file")
    {
        std::string filename = argv[arg_index++];
        std::ifstream ifile(filename, std::ios::in);
        if (!ifile.is_open())
        {
            std::cerr << "There was a problem opening the input file!\n";
            exit(1); // exit or do additional error checking
        }
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            ifile >> zetas[itr_beta];
        }
        ifile.close();
    }
    else if (flag_zeta_type == "equal")
    {
        double valuezeta = std::stod(argv[arg_index++]);
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            zetas[itr_beta] = valuezeta;
        }
    }
    else
    {
        print("error in initialize zetas.");
        return -1;
    }

    bool flag_C = false;
    // print(argv[arg_index]);
    std::string flag_C_string = argv[arg_index++];
    if (flag_C_string == "true" || flag_C_string == "1")
    {
        flag_C = true;
    }
    else if (flag_C_string == "false" || flag_C_string == "0")
    {
        flag_C = false;
    }
    else
    {
        print("error in initialize flag_C.");
        return -1;
    }

    double C = std::stod(argv[arg_index++]);

    Vec ps = std::vector<double>(num_comps, 1.0);
    std::string flag_ps_type = argv[arg_index++];
    if (flag_ps_type == "file")
    {
        std::string filename = argv[arg_index++];
        std::ifstream ifile(filename, std::ios::in);
        if (!ifile.is_open())
        {
            std::cerr << "There was a problem opening the input file!\n";
            exit(1); // exit or do additional error checking
        }
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            ifile >> ps[itr_beta];
        }
        ifile.close();
    }
    else if (flag_ps_type == "equal")
    {
        double valueps = std::stod(argv[arg_index++]);
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            ps[itr_beta] = valueps;
        }
    }
    else
    {
        print("error in initialize ps.");
        return -1;
    }

    bool flag_zetas = false;
    std::string flag_zeta_string = argv[arg_index++];
    if (flag_zeta_string == "true" || flag_zeta_string == "1")
    {
        flag_zetas = true;
    }
    else if (flag_zeta_string == "false" || flag_zeta_string == "0")
    {
        flag_zetas = false;
    }
    else
    {
        print("error in initialize flag_zeta.");
        return -1;
    }

    bool flag_ps = false;
    std::string flag_ps_string = argv[arg_index++];
    if (flag_ps_string == "true" || flag_ps_string == "1")
    {
        flag_ps = true;
    }
    else if (flag_ps_string == "false" || flag_ps_string == "0")
    {
        flag_ps = false;
    }
    else
    {
        print("error in initialize flag_ps.");
        return -1;
    }
    std::string filenamepre = argv[arg_index++];
    std::string flag_save_separate_string = argv[arg_index++];

    double threshold_incomp = 1e-10;
    double threshold_omega = 1e-8;
    double threshold_J = 1e-6;
    double threshold_Lbeta = 1e-8;
    double threshold_zeta = 1e-7;
    print(flag_save_separate_string);
    // print(argc);
    // print(arg_index);

    if (argc > arg_index + 1)
    {
        threshold_incomp = std::stod(argv[arg_index++]);
        threshold_omega = std::stod(argv[arg_index++]);
        threshold_J = std::stod(argv[arg_index++]);
        threshold_Lbeta = std::stod(argv[arg_index++]);
        threshold_zeta = std::stod(argv[arg_index++]);
    }

    print("phi_means=");
    print(phi_means);
    print("chis=");
    print(chis);
    print("Ls=");
    print(Ls);
    print("zs=");
    print(zs);
    print("kappas=");
    print(kappas);
    print("v=");
    print(v);
    // print("omegas=");
    // print(omegas);
    print("acceptance_omega");
    print(acceptance_omega);
    print("acceptance_J");
    print(acceptance_J);
    print("acceptance_Lbeta");
    print(acceptance_Lbeta);
    print("acceptance_zeta");
    print(acceptance_zeta);
    print("Js=");
    print(Js);
    print("Lbetas=");
    print(Lbetas);
    print("zetas=");
    print(zetas);
    print("flag_C=");
    print(flag_C);
    print("C=");
    print(C);
    print("ps=");
    print(ps);
    print("flag_zetas=");
    print(flag_zetas);
    print("flag_ps=");
    print(flag_ps);
    print("threshold_incomp=");
    print(threshold_incomp);
    print("threshold_omega=");
    print(threshold_omega);
    print("threshold_J=");
    print(threshold_J);
    print("threshold_Lbeta=");
    print(threshold_Lbeta);
    print("threshold_zeta=");
    print(threshold_zeta);
    print("flag_save_separate_string");
    print(flag_save_separate_string);

    // all acceptances are 0 except acceptance_omega
    std::tuple<double, double, double, double, double> errors1 = Gibbs_dynamics(
        phi_means,
        chis,
        Ls,
        zs,
        kappas,
        v,
        omegas,
        phis,
        psi,
        steps_inner1,
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

    std::tuple<double, double, double, double, double> errors = Gibbs_dynamics(
        phi_means,
        chis,
        Ls,
        zs,
        kappas,
        v,
        omegas,
        phis,
        psi,
        steps_inner2,
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
        threshold_zeta);

    Vec local_energy(num_beta, 0.0);
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

    // write results to file
    // std::string filenamepre = argv[arg_index++];
    std::ofstream myfile(filenamepre + "_alldata.txt", std::ofstream::out);
    if (myfile.is_open())
    {
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            myfile << Lbetas[itr_beta] << ' ';
        }
        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            myfile << Js[itr_beta] << ' ';
        }

        myfile << total_energy << ' ';
        myfile << std::get<0>(errors) << ' ';
        myfile << std::get<1>(errors) << ' ';
        myfile << std::get<2>(errors) << ' ';
        myfile << std::get<3>(errors) << ' ';
        myfile << std::get<4>(errors) << ' ';
        myfile << std::endl;
        for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                {
                    myfile << phis[itr_comp][itr_beta][itr_coord] << ' ';
                }
                myfile << std::endl;
            }
        }

        for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                {
                    myfile << omegas[itr_comp][itr_beta][itr_coord] << ' ';
                }
                myfile << std::endl;
            }
        }

        for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
        {
            for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
            {
                myfile << psi[itr_beta][itr_coord] << ' ';
            }
            myfile << std::endl;
        }
        myfile.close();
    }
    else
    {
        std::cout << "unable to open file.";
    }

    // save to different files
    if (flag_save_separate_string == "true" || flag_save_separate_string == "1")
    {
        std::ofstream myfile1;
        myfile1.open(filenamepre + "_Lbetas.txt", std::ofstream::out);
        if (myfile1.is_open())
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                myfile1 << std::setprecision(14);
                myfile1 << Lbetas[itr_beta] << ' ';
            }
            myfile1.close();
        }
        else
        {
            print("cannot open " + filenamepre + "_Lbetas.txt");
            return -1;
        }

        myfile1.open(filenamepre + "_Js.txt", std::ofstream::out);
        if (myfile1.is_open())
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                myfile1 << std::setprecision(14);
                myfile1 << Js[itr_beta] << ' ';
            }
            myfile1.close();
        }
        else
        {
            print("cannot open " + filenamepre + "_Js.txt");
            return -1;
        }

        myfile1.open(filenamepre + "_phis.txt", std::ofstream::out);
        if (myfile1.is_open())
        {
            for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
            {
                for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
                {
                    for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                    {
                        myfile1 << std::setprecision(14);
                        myfile1 << phis[itr_comp][itr_beta][itr_coord] << ' ';
                    }
                    myfile1 << std::endl;
                }
            }

            myfile1.close();
        }

        myfile1.open(filenamepre + "_omegas.txt", std::ofstream::out);
        if (myfile1.is_open())
        {
            for (int itr_comp = 0; itr_comp < num_comps; itr_comp++)
            {
                for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
                {
                    for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                    {
                        myfile1 << std::setprecision(14);
                        myfile1 << omegas[itr_comp][itr_beta][itr_coord] << ' ';
                    }
                    myfile1 << std::endl;
                }
            }

            myfile1.close();
        }

        myfile1.open(filenamepre + "_psi.txt", std::ofstream::out);
        if (myfile1.is_open())
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                for (int itr_coord = 0; itr_coord < num_coord; itr_coord++)
                {
                    myfile1 << std::setprecision(14);
                    myfile1 << psi[itr_beta][itr_coord] << ' ';
                }
                myfile1 << std::endl;
            }
            myfile1.close();
        }
        myfile1.open(filenamepre + "_zetas.txt", std::ofstream::out);
        if (myfile1.is_open())
        {
            for (int itr_beta = 0; itr_beta < num_beta; itr_beta++)
            {
                myfile1 << std::setprecision(14);
                myfile1 << zetas[itr_beta] << ' ';
            }
            myfile1.close();
        }
    }
    return 0;
}