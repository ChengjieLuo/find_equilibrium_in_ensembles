#include "gibbs_dynamics.h"

// Example usage of the Convolution function
int main(int argc, char *argv[])
{

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
    double chi = std::stod(argv[arg_index++]);
    chis[0][1] = chi;
    chis[1][0] = chi;
    if (chis[0][1] < 0)
    {
        for (int i = 0; i < num_comps; i++)
        {
            for (int j = 0; j < num_comps; j++)
            {
                chis[i][j] -= chi;
            }
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
                double tmp=distribution(generator);
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
    print(argc);
    print(arg_index);

    if(argc>arg_index+1)
    {
        threshold_incomp=std::stod(argv[arg_index++]);
        threshold_omega =std::stod(argv[arg_index++]);
        threshold_J =std::stod(argv[arg_index++]);
        threshold_Lbeta =std::stod(argv[arg_index++]);
        threshold_zeta =std::stod(argv[arg_index++]);
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
    print("omegas=");
    print(omegas);
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
    std::tuple<double, double, double, double, double> errors1=Gibbs_dynamics(
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

    std::tuple<double, double, double, double, double> errors=Gibbs_dynamics(
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
        myfile<<total_energy<<' ';
        myfile<<std::get<0>(errors)<<' ';
        myfile<<std::get<1>(errors)<<' ';
        myfile<<std::get<2>(errors)<<' ';
        myfile<<std::get<3>(errors)<<' ';
        myfile<<std::get<4>(errors)<<' ';
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
                    myfile1 << zetas[itr_beta]<< ' ';
            }
            myfile1.close();
        }
    }
    return 0;
}