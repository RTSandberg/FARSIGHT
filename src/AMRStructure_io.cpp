#include "AMRStructure.hpp"

int AMRStructure::write_particles_to_file() {
    std::ofstream x_file;
    std::ofstream v_file;
    std::ofstream f_file;
    std::ofstream qw_file;
    std::ofstream e_file;

    x_file.open(sim_dir + "simulation_output/xs/xs_" + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    v_file.open(sim_dir + "simulation_output/vs/vs_" + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    f_file.open(sim_dir + "simulation_output/fs/fs_" + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    qw_file.open(sim_dir + "simulation_output/qws/qws_" + std::to_string(iter_num), std::ios::out | std::ios::binary); 
    e_file.open(sim_dir + "simulation_output/es/es_" + std::to_string(iter_num), std::ios::out | std::ios::binary); 

    // std::cout << "#xs " << xs.size() << std::endl;
    // std::cout << "#vs " << vs.size() << std::endl;
    // std::cout << "#fs " << fs.size() << std::endl;
    // std::cout << "#qw " << q_ws.size() << std::endl;
    // std::cout << "#es " << es.size() << std::endl;

    if (!x_file | !v_file | !f_file | !qw_file | !e_file ) {
        cout << "Unable to open step " << iter_num << " particle data files" << endl;
        return 1;
    }

    // while copying seems to be one of the best ways, I can't make it work right now.
    // copy(particles.begin(), particles.end(), std::ostreambuf_iterator<char>(out_file));

    // for (int ii = 0; ii < particles.size(); ++ii){
    assert(xs.size() == vs.size() && vs.size() == fs.size() && fs.size() == q_ws.size() && q_ws.size() == es.size());
    for (int ii = 0; ii < xs.size(); ++ii) {
        // double x = particles[ii].get_x();
        // double v = particles[ii].get_v();
        // double f = particles[ii].get_f();
        // double qw = particles[ii].get_qw();
        // double e = particles[ii].get_e();
        double x = xs[ii];
        double v = vs[ii];
        double f = fs[ii];
        double qw = q_ws[ii];
        double e = es[ii];
        x_file.write((char *) &x, sizeof(double));
        v_file.write((char *) &v, sizeof(double));
        f_file.write((char *) &f, sizeof(double));
        qw_file.write((char *) &qw, sizeof(double));
        e_file.write((char *) &e, sizeof(double));
    }

    if (!x_file.good() | !v_file.good() | !f_file.good() | !qw_file.good() | !e_file.good()) {
        cout << "Error occurred writing step " << iter_num << " particle data files." << endl;
        return 1;
    }
    x_file.close();
    v_file.close();
    f_file.close();
    qw_file.close();
    e_file.close();
    // cout << "Successfully wrote step " << iter_num << " particle data files" << endl;

    return 0;
}

int AMRStructure::write_panels_to_file() {
    std::ofstream panel_file;

    panel_file.open(sim_dir + "simulation_output/panels/leaf_point_inds_" + std::to_string(iter_num), ios::out | ios::binary);
    
    if (!panel_file) {
        cout << "Unable to open step " << iter_num << " panel data files" << endl;
        return 1;
    }

    for (int ii=0; ii < panels.size(); ++ii ) {
        const int *inds = panels[ii].point_inds;
        panel_file.write( (char*) inds, 9*sizeof(int));
    }

    panel_file.close();
    if (!panel_file.good() ) {
        cout << "Error writing step " << iter_num << " panel data files" << endl;
        return 1;
    }

    // cout << "Successfully wrote step " << iter_num << " panel data files" << endl;

    return 0;
}

int AMRStructure::write_to_file() {
    write_particles_to_file();
    write_panels_to_file();
    return 0;
}



void AMRStructure::print_amr() {
    printf("This mesh structure has %lu panels.\n", panels.size() );
    for (const Panel& panel : panels) {
        panel.print_panel();
    }
    // printf("\n structure particle data is \n");
    // for (int ii = 0; ii < particles.size(); ++ii) {
    //     printf("particle %i ", ii);
    //     particles[ii].print_particle();
    // }
}
// AMRStructure_io.cpp
std::ostream& operator<<(std::ostream& os, const AMRStructure& amr) {
    os << "=================" << endl;
    os << "This is an AMR structure" << endl;
    os << "=================" << endl;
    os << "Computational domain: (x,v) in [" << amr.x_min << ", " << amr.x_max << "]x[" << amr.v_min << ", " << amr.v_max << "]" << endl; 
    os << "Species charge: " << amr.q << ", species mass: " << amr.q/amr.qm << endl;
    os << "-----------------" << endl;
    os << "Initial height: " << amr.initial_height << ", height: " << amr.height << ", max height: " << amr.max_height << endl;
    if (amr.do_adaptively_refine) { os << "This structure is adaptively refined" << endl;}
    else { os << "This structure is NOT adaptively refined" << endl;}

    if (amr.need_further_refinement) { os << "Structure is not done being refined" << endl;}
    else { os << "Structure has reached full refinement" << endl; }


    os << "---------" << endl;
    // os << "time stepping parameters" << endl << "--------" << endl;
    os << "num steps: " <<  amr.iter_num << ", dt: " << amr.dt << endl;
    os << "---------" << endl;

    // os << "Field parameters " << endl << "---------" << endl;
    os << "Field softening parameter (aka Green's function epsilon): " << amr.greens_epsilon << endl;
    os << "---------" << endl;

    os << "Panels: " << endl << "==============" << endl;
    std::copy(amr.panels.begin(), amr.panels.end(), std::ostream_iterator<Panel>(os));
    os <<  "==============" << endl;

    os << "Point data" << endl << "==============" << endl;
    os << "xs: ";
    std::copy(amr.xs.begin(), amr.xs.end(), std::ostream_iterator<double>(os, " "));
    os << endl << "vs: ";
    std::copy(amr.vs.begin(), amr.vs.end(), std::ostream_iterator<double>(os, " "));
    os << endl << "fs: ";
    std::copy(amr.fs.begin(), amr.fs.end(), std::ostream_iterator<double>(os, " "));
    os << endl << "qw: ";
    std::copy(amr.q_ws.begin(), amr.q_ws.end(), std::ostream_iterator<double>(os, " "));
    os << endl << "es: ";
    std::copy(amr.es.begin(), amr.es.end(), std::ostream_iterator<double>(os, " "));

    return os;
}

void AMRStructure::print_panel_points() {
    auto panel_it = panels.begin();

    for (int ii = 0; ii < panels.size(); ii++) {
        cout << "==============" << endl;
        cout << "Panel " << ii << " points" << endl;
        cout << "xs: ";
        for (int jj = 0; jj < 9; jj++) {
            cout << xs[panel_it->point_inds[jj]] << ", ";
        }
        cout << endl << "vs: ";
        for (int jj = 0; jj < 9; jj++) {
            cout << vs[panel_it->point_inds[jj]] << ", ";
        }
        cout << endl << "fs: ";
        for (int jj = 0; jj < 9; jj++) {
            cout << fs[panel_it->point_inds[jj]] << ", ";
        }
        cout << endl << "qw: ";
        for (int jj = 0; jj < 9; jj++) {
            cout << q_ws[panel_it->point_inds[jj]] << ", ";
        }
        cout << endl;
        panel_it++;
    }
}
