#include "AMRSimulation.hpp"

int AMRSimulation::load_deck(std::string &deck_address, pt::ptree &deck) {

    try {
        pt::read_json(deck_address, deck);
    } catch(std::exception& err) {
        cout << "unable to open input deck" << endl;
        return 1;
    }
    return 0;
}

int AMRSimulation::get_box_t_params(pt::ptree &deck) {

    std::string project_name = deck.get<std::string>("sim_name", "no_name_found");
    x_min = deck.get<double>("xmin", 0.0), x_max = deck.get<double>("xmax", 0.0);
    p_min = deck.get<double>("pmin", -1.0), p_max = deck.get<double>("pmax",1.0);
    int bcs_int = deck.get<int>("bcs",0);
    
    if (bcs_int < 0 || bcs_int >= last_bc) {
        cout << "Invalid boundary condition provided in input deck. Valid BCs are: " << endl;
        cout << "0: periodic, 1: open" << endl;
        return 1;
    }
    
    bcs = static_cast<BoundaryConditions> (bcs_int);

    Lx = x_max - x_min;

    int which_quad = deck.get<int>("quadrature",1);
    quad = static_cast<Quadrature>(which_quad);

    int relativistic_int = deck.get<int>("relativistic", 0);
    relativistic = (relativistic_int > 0);
    num_steps = deck.get<int>("num_steps", 10);//atoi(argv[19]);//120;
    n_steps_remesh = deck.get<int>("remesh_period", 1); //atoi(argv[20]);
    n_steps_diag = deck.get<int>("diag_period", 1); //atoi(argv[21]);
    dt = deck.get<double> ("dt", 0.5); //atof(argv[22]);//0.5;

    return 0;
}

ElectricField* AMRSimulation::make_field_return_ptr(pt::ptree &deck) {
    
    ElectricField* calculate_e;

    double greens_epsilon = deck.get<double>("greens_epsilon",0.2);//atof(argv[12]);//0.2;
    int use_treecode = deck.get<int>("use_treecode", 0); //atoi(argv[13]);
    double beta = deck.get<double>("beta", -1.0); //atof(argv[14]);
    double mac = deck.get<double>("mac", -1.0); //atof(argv[15]);
    int degree = deck.get<int>("degree", -1); //atoi(argv[16]);
    int max_source = deck.get<int>("max_source", 2000); //atoi(argv[17]);
    int max_target = deck.get<int>("max_target", 2000); //atoi(argv[18]);
    
    if (use_treecode > 0) {
        if (bcs!=periodic_bcs) { // open_bcs
            if (0 <= beta && beta <= 1.0) {
                calculate_e = new E_MQ_Treecode_openbcs(greens_epsilon, beta);
            } else {
                int verbosity = 0;
                calculate_e = new E_MQ_Treecode_openbcs(greens_epsilon, 
                                                mac, degree, 
                                                max_source, max_target, 
                                                verbosity);
            }
        } else {
            if (0 <= beta && beta <= 1.0) {
                calculate_e = new E_MQ_Treecode(Lx, greens_epsilon, beta);
            } else {
                int verbosity = 0;
                calculate_e = new E_MQ_Treecode(Lx, greens_epsilon, 
                                                mac, degree, 
                                                max_source, max_target, 
                                                verbosity);
            }
        }
    } else {
        if (bcs==periodic_bcs) {
            calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
        } else { //open bcs
            calculate_e = new E_MQ_DirectSum_openbcs(greens_epsilon);
        }
    }

    return calculate_e;
}


void AMRSimulation::make_external_field(pt::ptree &deck) {
    try {
        pt::ptree &external_field_deck = deck.get_child("external_field");

        int use_external_int = external_field_deck.get<int>("use_external_field", 0);
        if (use_external_int > 0) {
            use_external_field = true;
            int external_type = external_field_deck.get<int>("external_type", 0);
            double Am = external_field_deck.get<double>("Am");
            double k = external_field_deck.get<double>("k");
            double omega = external_field_deck.get<double>("omega");
            switch (external_type)
            {
                case 1: {
                    calculate_e_external = new ExternalLogistic(Am, k, omega);
                }
                break;
            
                default: // 0
                {
                    calculate_e_external = new ExternalSine(Am, k, omega);
                }
                break;
            }
        }
        else {
            use_external_field = false;
            calculate_e_external = nullptr;
        }
    } catch(std::exception& e) {
        cout << "External field parameters not found.  External field disabled." << endl;
        use_external_field = false;
        calculate_e_external = nullptr;
    }
}

distribution* AMRSimulation::make_f0_return_ptr(pt::ptree &species_deck_portion) {
    pt::ptree deck = species_deck_portion;

    double kx = 2.0 * M_PI / Lx * deck.get<double>("normalized_wavenumber",1.0);
    double amp = deck.get<double>("amp", 0.0);//0.5;
    double pth = deck.get<double>("pth", 1.0);//atof(argv[9]);//1.0;
    double pstr = deck.get<double>("pstr", 0.0); //atof(argv[10]);
    int ics_type = deck.get<int>("ics_type", 1);//atoi(argv[6]);
    distribution* f0;
    switch (ics_type)
    {
        case 1: // weak Landau Damping
            f0 = new F0_LD(pth, pstr, kx, amp);
            break;
        case 2: // strong Landau Damping
            f0 = new F0_LD(pth, pstr, kx, amp);
            break;
        case 3: // 'strong' two-stream
            f0 = new F0_strong_two_stream(pth, kx, amp);
            break;
        case 4: // 'colder' two-stream
            f0 = new F0_colder_two_stream(pth, pstr, kx, amp);
            break;
        case 5: // Friedman beam problem 
        {
            double Tstar = pth * pth;
            f0 = new F0_Friedman_beam(amp, Tstar, x_max);
        }
            break;
        case 6: // maxwell-Jutner
        {
            double mu;
            try {
                mu = deck.get<double>("mu");
            }
            catch (pt::ptree_error& mu_error) {
                cout << "relativistic inverse temperature mu not found in deck" << endl;
                cout << "setting mu to default, 1.0" << endl;
                mu = 1.0;
            }
            f0 = new Maxwell_Jutner(mu, pstr, kx, amp);
        }
            break;
        case 7: // relativistic two-stream(double p_th, double p_str, double k, double amp)
        {
            double mu;
            try {
                mu = deck.get<double>("mu");
            }
            catch (pt::ptree_error& mu_error) {
                cout << "relativistic inverse temperature mu not found in deck" << endl;
                cout << "setting mu to default, 1.0" << endl;
                mu = 1.0;
            }
            f0 = new Relativistic_Two_Stream(pstr, mu, kx, amp);
        }
            break;
        case 8: // relativistic traveling wave
            double wave_beta;
            try {
                wave_beta = deck.get<double>("wave_beta");
            }
            catch (pt::ptree_error& wave_beta_error) {
                cout << "wave beta not found in deck" << endl;
                cout << "setting wave_beta to default, 0.8" << endl;
                wave_beta = 0.8;
            }
            f0 = new Relativistic_Wave(amp, wave_beta, pth);
            break;
        default:
            cout << "Using default (Landau damping) initial conditions" << endl;
            f0 = new F0_LD(pth, pstr, kx, amp);
            break;
    }
    return f0;
}

AMRStructure* AMRSimulation::make_species_return_ptr(pt::ptree &species_deck_portion, distribution* f0) {

    cout << "Trying to load species" << endl;

    pt::ptree deck = species_deck_portion;
    // get parameters
    std::string sp_name = deck.get_child("name").get_value<std::string>();
    double kx = 2.0 * M_PI / Lx * deck.get<double>("normalized_wavenumber",1.0);
    double amp = deck.get<double>("amp", 0.0);//0.5;
    double pth = deck.get<double>("pth", 1.0);//atof(argv[9]);//1.0;
    double pstr = deck.get<double>("pstr", 0.0); //atof(argv[10]);
    int sim_type = deck.get<int>("sim_type", 1);//atoi(argv[6]);
    double q = deck.get<double>("q", -1.0);
    double m = deck.get<double>("m", 1.0);
    if (sim_type==5) { q = 1.0; }
    int initial_height = deck.get<int>("initial_height",6);//atoi(argv[11]);//6;
    int p_height = deck.get<int>("p_height",0);
    int max_height = deck.get<int>("max_height", initial_height);
    bool do_adaptively_refine = deck.get<bool> ("adaptively_refine", false);//false;
    std::vector<double> amr_epsilons; 
    try {
        for (pt::ptree::value_type &eps : deck.get_child("amr_epsilons")) {
            amr_epsilons.push_back(eps.second.get_value<double>() );
        }
    } catch(std::exception& err) {
        cout << "Unable to find amr refinement values.  Disabling amr." << endl;
        amr_epsilons = std::vector<double>();

        do_adaptively_refine = false;
    }


    AMRStructure *species = new AMRStructure{sim_dir, sp_name,
                f0, q, m,
                initial_height, p_height,max_height,
                x_min, x_max, p_min, p_max, bcs,
                quad, 
                do_adaptively_refine, amr_epsilons};
    // 
    return species;
}

void AMRSimulation::get_qms() {
    for(auto &species : species_list) {
        species_qs.push_back(species->q);
        species_qms.push_back(species->qm);
        species_ms.push_back(species->q / species->qm);
    }
}

    // distribution* f0;
    // double q = -1.0, m = 1.0;

    /*
    int which_quad = deck.get<int>("quadrature",1);
    quad = static_cast<Quadrature>(which_quad);
    // if (which_quad == 1) {
    //     quad = trap;
    // } else {
    //     quad = simpsons;
    // }

    // int initial_height = deck.get<int>("initial_height",6);//atoi(argv[11]);//6;
    // int v_height = deck.get<int>("v_height",0);
    // int max_height = deck.get<int>("max_height", initial_height);
    greens_epsilon = deck.get<double>("greens_epsilon",0.2);//atof(argv[12]);//0.2;
    use_treecode = deck.get<int>("use_treecode", 0); //atoi(argv[13]);
    beta = deck.get<double>("beta", -1.0); //atof(argv[14]);
    mac = deck.get<double>("mac", -1.0); //atof(argv[15]);
    degree = deck.get<int>("degree", -1); //atoi(argv[16]);
    max_source = deck.get<int>("max_source", 2000); //atoi(argv[17]);
    max_target = deck.get<int>("max_target", 2000); //atoi(argv[18]);

    // int nxsqrt = pow(2, initial_height + 1) + 1;
    // int nx = nxsqrt * nxsqrt;

    ElectricField* calculate_e;
    if (use_treecode > 0) {
        if (bcs!=periodic_bcs) { // open_bcs
            if (0 <= beta && beta <= 1.0) {
                calculate_e = new E_MQ_Treecode_openbcs(greens_epsilon, beta);
            } else {
                int verbosity = 0;
                calculate_e = new E_MQ_Treecode_openbcs(greens_epsilon, 
                                                mac, degree, 
                                                max_source, max_target, 
                                                verbosity);
            }
        } else {
            if (0 <= beta && beta <= 1.0) {
                calculate_e = new E_MQ_Treecode(Lx, greens_epsilon, beta);
            } else {
                int verbosity = 0;
                calculate_e = new E_MQ_Treecode(Lx, greens_epsilon, 
                                                mac, degree, 
                                                max_source, max_target, 
                                                verbosity);
            }
        }
    } else {
        if (bcs==periodic_bcs) {
            calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
        } else { //open bcs
            calculate_e = new E_MQ_DirectSum_openbcs(greens_epsilon);
        }
    }

    int num_steps = deck.get<int>("num_steps", 10);//atoi(argv[19]);//120;
    int n_steps_remesh = deck.get<int>("remesh_period", 1); //atoi(argv[20]);
    int n_steps_diag = deck.get<int>("diag_period", 1); //atoi(argv[21]);
    double dt = deck.get<double> ("dt", 0.5); //atof(argv[22]);//0.5;
    bool do_adaptively_refine = deck.get<bool> ("adaptively_refine", false);//false;


    std::vector<double> amr_epsilons; 
    try {
        for (pt::ptree::value_type &eps : deck.get_child("amr_epsilons")) {
            amr_epsilons.push_back(eps.second.get_value<double>() );
        }
    } catch(std::exception& err) {
        cout << "Unable to find amr refinement values.  Disabling amr." << endl;
        amr_epsilons = std::vector<double>();

        do_adaptively_refine = false;
    }

    cout << "============================" << endl;
    cout << "Running a FARRSIGHT simulation" << endl;
    cout << "sim dir: " << sim_dir << endl;
    cout << "deck found in: " << input_deck << endl;
    cout << x_min << " <= x <= " << x_max << endl;
    cout << v_min << " <= v <= " << v_max << endl;
    switch (bcs) {
        case (open_bcs) : cout << "Using open boundary conditions" << endl;
            break;
        default : // periodic
            cout << "Using periodic boundary conditions" << endl;
            break;
    }
    cout << "k = " << kx << ", amp = " << amp << ", vth = " << vth << ", vstr = " << vstr <<  endl;
    cout << "height " << initial_height << ", v height " << v_height << endl;
    switch (quad) {
        case (simpsons) : cout << "Using Simpson's rule" << endl;
            break;
        default : cout << "Using trap rule" << endl;
            break;
    }
    cout << "green's epsilon = " << greens_epsilon << endl;
    cout << "Taking " << num_steps << " steps with dt = " << dt << endl;
    cout << "Remesh every " << n_steps_remesh << " step(s), diagnostic dump every " << n_steps_diag << " step(s)" << endl;
    cout << "use treecode flag " << use_treecode << endl;
    if (use_treecode > 0) { 
        if (0 <= beta && beta <= 1.0) {
            cout << "Using treecode with beta " << beta << endl;
        } else {
            cout << "Using treecode with mac " << mac << " and degree " << degree << endl;
        }
    } else {
        cout << "using direct sum" << endl;
    }
    if (do_adaptively_refine) {
        cout << "Adaptively refining, to height at most " << max_height << endl;
        cout << "Refinement epsilon(s) : ";
        std::copy(amr_epsilons.begin(), amr_epsilons.end(), std::ostream_iterator<double>(cout, " "));
        cout << endl;
    } else {
        cout <<"Not adaptively refining." << endl;
    }
    cout << "============================" << endl;

    auto sim_start = high_resolution_clock::now();
    


    // int v_height = 2;
    // initial_height = 0;

    // AMRStructure amr{sim_dir, f0, q, m,
    //             initial_height, v_height,max_height,
    //             x_min, x_max, v_min, v_max, bcs,
    //             calculate_e, quad, num_steps, dt,
    //             do_adaptively_refine, amr_epsilons};


    return 0;
}
*/
