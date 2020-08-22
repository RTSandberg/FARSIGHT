#include "AMRStructure.hpp"


int AMRStructure::create_prerefined_mesh() {
    // printf("setting initial mesh of height %i.\n", initial_height);
    if (initial_height < 1) {
        throw std::invalid_argument("height must be greater than 1");
    }
    double dx = (x_max - x_min) / 4;
    double dv = (v_max - v_min) / 4;
    std::vector<double> xs_init, vs_init;
    for(int ii = 0; ii < 5; ++ii) {
        xs_init.push_back(ii * dx);
        vs_init.push_back(v_min + ii * dv);
    }
    panels.clear();
    xs.clear();
    vs.clear();
    fs.clear();
    xs.reserve(25);
    vs.reserve(25);
    for (int ii = 0; ii < 5; ii += 2) {
        for (int jj = 0; jj < 5; jj+=2) { 
            xs.push_back(xs_init[ii]); 
            vs.push_back(vs_init[jj]);
        }
    }
    for (int ii = 0; ii < 2; ++ii) {
        xs.push_back(xs_init[2*ii]); xs.push_back(xs_init[2*ii]);
        vs.push_back(vs_init[1]); vs.push_back(vs_init[3]);
        for (int jj = 0; jj < 5; ++jj) {
            xs.push_back(xs_init[1 + 2*ii]);
            vs.push_back(vs_init[jj]);
        }
    }
    for (int jj = 1; jj < 5; jj+=2) {
        xs.push_back(xs_init[4]);
        vs.push_back(vs_init[jj]);
    }
    fs = std::vector<double>(25, 1.0);
    //   2   15      5    22    8 (2 by periodic bcs)
    
    //   10  14[2]   17   21[4] 24 (10 by pbcs)
    
    //   1   13      4    20     7 (1 by periodic bcs)
    
    //   9   12[1]   16   19[3] 23 (9 by pbcs)
    
    //   0   11     3    18     6 (0 by periodic bcs)

    panels.push_back(Panel{});
    panels[0].is_left_bdry = true;
    panels[0].is_right_bdry = true;
    panels.push_back(Panel{1, 1, 0, 0, 3, 2, 3, -2});
    panels[1].is_left_bdry = true;
    panels[1].set_point_inds(0,9,1,11,12,13,3,16,4);
    panels[1].needs_refinement = true;
    panels.push_back(Panel{2, 1, 0, 1, 4, -2, 4, 1});
    panels[2].is_left_bdry = true;
    panels[2].set_point_inds(1,10,2,13,14,15,4,17,5);
    panels[2].needs_refinement = true;
    panels.push_back(Panel{3, 1, 0, 2, 1, 4, 1, -2});
    panels[3].is_right_bdry = true;
    panels[3].set_point_inds(3,16,4,18,19,20,6,23,7);
    panels[3].needs_refinement = true;
    panels.push_back(Panel{4, 1, 0, 3, 2,-2,2,3});
    panels[4].is_right_bdry = true;
    panels[4].set_point_inds(4,17,5,20,21,22,7,24,8);
    panels[4].needs_refinement = true;

    panels[0].set_child_inds_start(1);
    minimum_unrefined_index = 1;

    // call refine
    for (int level = 1; level < initial_height; ++level) {
        int num_panels_pre_refine = panels.size();

        for (auto panel_it = panels.begin() + minimum_unrefined_index; panel_it != panels.end(); ++panel_it) {
            panel_it->needs_refinement = true;
        }
        refine_panels( [] (double x, double v) {return 1.0;} , false);
        minimum_unrefined_index = num_panels_pre_refine;
    }
    is_initial_mesh_set = true;

    return 0;
}

void AMRStructure::refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine) {
    std::vector <double> new_xs;
    std::vector <double> new_vs;
    std::vector <double> new_fs;
    std::vector <int> prospective_leaf_inds;
    // int new_vert_ind = particles.size();
    int new_vert_ind = xs.size();
    int num_panels_before_this_iter = panels.size();
    


    for (int jj = minimum_unrefined_index; jj < num_panels_before_this_iter; ++jj) {
        Panel* panel= &(panels[jj]);
        if (panel->needs_refinement ) {
            // printf("refining panel %i\n", panel->get_panel_ind() );
            std::vector<double> panel_xs;
            std::vector<double> panel_vs;
            double dx, dv;

            const int* panel_points = panel->point_inds;
            for (int ii = 0; ii < 9; ++ii) {
                int point_ind = panel_points[ii];
                // panel_xs.push_back(particles[vertex_ind].get_x() );
                // panel_vs.push_back(particles[vertex_ind].get_v() );
                panel_xs.push_back(xs[point_ind]);
                panel_vs.push_back(vs[point_ind]);
            }
            dx = panel_xs[3] - panel_xs[0];
            dv = panel_vs[1] - panel_vs[0];
            double sub_dx = 0.5 * dx;
            double sub_dv = 0.5 * dv;

            int num_new_panels = panels.size();
            // double x_left = panel_xs[0];
            // double x_mid = x_left + .5 * dx;
            // double x_right = panel_xs[2];
            // double v_bottom = panel_vs[0];
            // double v_mid = v_bottom + .5 * dv;
            // double v_top = panel_vs[1];
            double subpanel_xs[5], subpanel_vs[5];
            for (int ii = 0; ii < 5; ii ++) {
                subpanel_xs[ii] = panel_xs[0] + sub_dx * ii;
                subpanel_vs[ii] = panel_vs[0] + sub_dv * ii;
            }

            // int left_vert_ind, bottom_vert_ind, mid_vert_ind, top_vert_ind, right_vert_ind;
            int point_9_ind, point_10_ind, point_11_ind, point_15_ind, 
                point_18_ind, point_22_ind, point_23_ind, point_24_ind;
            int child_0_bottom_nbr_ind = -1;
            int child_0_left_nbr_ind = -1;
            int child_1_left_nbr_ind = -1;
            int child_1_top_nbr_ind = -1;
            int child_2_bottom_nbr_ind = -1;
            int child_2_right_nbr_ind = -1;
            int child_3_right_nbr_ind = -1;
            int child_3_top_nbr_ind = -1;

            // generate new vertices

            Panel* panel_parent;
            // check left neighbor
            if (panel->left_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_left = &(panels[panel_parent->left_nbr_ind]);
                if (!parent_left->is_refined) {
                    parent_left->needs_refinement = true;
                    need_further_refinement = true;
                }
                point_9_ind = new_vert_ind;
                // point_10_ind = new_vert_ind++;
                point_10_ind = point_9_ind + 1;
                new_vert_ind += 2;
                new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                new_vs.push_back(subpanel_vs[1]); new_vs.push_back(subpanel_vs[3]);
            } else {
                Panel* panel_left = &(panels[panel->left_nbr_ind]);
                if (! panel_left->is_refined) {
                    point_9_ind = new_vert_ind++;
                    point_10_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                    new_vs.push_back(subpanel_vs[1]); new_vs.push_back(subpanel_vs[3]);
                }
                else {
                    child_0_left_nbr_ind = panel_left->child_inds_start +2;
                    Panel* child_0_left_nbr = &(panels[child_0_left_nbr_ind]);
                    child_0_left_nbr->right_nbr_ind = num_new_panels;
                    child_1_left_nbr_ind = panel_left->child_inds_start + 3;
                    Panel* child_1_left_nbr = &(panels[child_1_left_nbr_ind]);
                    child_1_left_nbr->right_nbr_ind = num_new_panels + 1;
                    if (panel->is_left_bdry) {
                        point_9_ind = new_vert_ind++;
                        point_10_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[0]); new_xs.push_back(subpanel_xs[0]);
                        new_vs.push_back(subpanel_vs[1]); new_vs.push_back(subpanel_vs[3]);
                    } else {
                        point_9_ind = child_0_left_nbr->point_inds[7];
                        point_10_ind = child_1_left_nbr->point_inds[7];
                    }
                }
            }
            // check bottom neighbor
            if (panel->bottom_nbr_ind == -2) {
                child_0_bottom_nbr_ind = -2;
                child_2_bottom_nbr_ind = -2;
                point_11_ind = new_vert_ind++;
                point_18_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_vs.push_back(subpanel_vs[0]); new_vs.push_back(subpanel_vs[0]);
            } else if (panel->bottom_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                
                Panel* parent_bottom = &(panels[panel_parent->bottom_nbr_ind]);
                if (!parent_bottom->is_refined ) {
                    parent_bottom->needs_refinement = true;
                    need_further_refinement = true;
                }
                point_11_ind = new_vert_ind++;
                point_18_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_vs.push_back(subpanel_vs[0]); new_vs.push_back(subpanel_vs[0]);
            } else {
                Panel* panel_bottom = &(panels[panel->bottom_nbr_ind]);
                if (! panel_bottom->is_refined ) {
                    point_11_ind = new_vert_ind++;
                    point_18_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                    new_vs.push_back(subpanel_vs[0]); new_vs.push_back(subpanel_vs[0]);
                }
                else {
                    child_0_bottom_nbr_ind = panel_bottom->child_inds_start + 1;
                    Panel* child_0_bottom_nbr = &(panels[child_0_bottom_nbr_ind]);
                    child_0_bottom_nbr->top_nbr_ind = num_new_panels;
                    child_2_bottom_nbr_ind = panel_bottom->child_inds_start +3;
                    Panel* child_2_bottom_nbr = &(panels[child_2_bottom_nbr_ind]);
                    child_2_bottom_nbr->top_nbr_ind = num_new_panels + 2;

                    point_11_ind = child_0_bottom_nbr->point_inds[5];
                    point_18_ind = child_2_bottom_nbr->point_inds[5];
                }
            }

            // check top neighbor
            if (panel->top_nbr_ind == -2) {
                child_1_top_nbr_ind = -2;
                child_3_top_nbr_ind = -2;
                point_15_ind = new_vert_ind++;
                point_22_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_vs.push_back(subpanel_vs[4]); new_vs.push_back(subpanel_vs[4]);
            } else if (panel->top_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                
                Panel* parent_top = &(panels[panel_parent->top_nbr_ind]);
                if (!parent_top->is_refined ) {
                    parent_top->needs_refinement = true;
                    need_further_refinement = true;
                }
                point_15_ind = new_vert_ind++;
                point_22_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                new_vs.push_back(subpanel_vs[4]); new_vs.push_back(subpanel_vs[4]);
            } else {
                Panel* panel_top = &(panels[panel->top_nbr_ind]);
                if (! panel_top->is_refined ) {
                    point_15_ind = new_vert_ind++;
                    point_22_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[1]); new_xs.push_back(subpanel_xs[3]);
                    new_vs.push_back(subpanel_vs[4]); new_vs.push_back(subpanel_vs[4]);
                }
                else {
                    child_1_top_nbr_ind = panel_top->child_inds_start;
                    Panel* child_1_top_nbr = &(panels[child_1_top_nbr_ind]);
                    child_1_top_nbr->bottom_nbr_ind = num_new_panels + 1;
                    child_3_top_nbr_ind = panel_top->child_inds_start + 2;
                    Panel* child_3_top_nbr = &(panels[child_3_top_nbr_ind]);
                    child_3_top_nbr->bottom_nbr_ind = num_new_panels + 3;

                    point_15_ind = child_1_top_nbr->point_inds[3];
                    point_22_ind = child_3_top_nbr->point_inds[3];
                }
            }

            // check right neighbor
            if (panel->right_nbr_ind == -1) {
                panel_parent = &(panels[panel->parent_ind]);
                Panel* parent_right = &(panels[panel_parent->right_nbr_ind]);
                if (!parent_right->is_refined ) {
                    parent_right->needs_refinement = true;
                    need_further_refinement = true;
                }
                point_23_ind = new_vert_ind++;
                point_24_ind = new_vert_ind++;
                new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                new_vs.push_back(subpanel_vs[1]); new_vs.push_back(subpanel_vs[3]);
            } else {
                Panel* panel_right = &(panels[panel->right_nbr_ind]);
                if (! panel_right->is_refined) {
                    point_23_ind = new_vert_ind++;
                    point_24_ind = new_vert_ind++;
                    new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                    new_vs.push_back(subpanel_vs[1]); new_vs.push_back(subpanel_vs[3]);
                }
                else {
                    child_2_right_nbr_ind = panel_right->child_inds_start;
                    Panel* child_2_right_nbr = &(panels[child_2_right_nbr_ind]);
                    child_2_right_nbr->left_nbr_ind = num_new_panels+2;
                    child_3_right_nbr_ind = panel_right->child_inds_start + 1;
                    Panel* child_3_right_nbr = &(panels[child_3_right_nbr_ind]);
                    child_3_right_nbr->left_nbr_ind = num_new_panels + 3;

                    if (panel->is_right_bdry) {
                        point_23_ind = new_vert_ind++;
                        point_24_ind = new_vert_ind++;
                        new_xs.push_back(subpanel_xs[4]); new_xs.push_back(subpanel_xs[4]);
                        new_vs.push_back(subpanel_vs[1]); new_vs.push_back(subpanel_vs[3]);
                    } else {
                        point_23_ind = child_2_right_nbr->point_inds[1];
                        point_24_ind = child_3_right_nbr->point_inds[1];
                    }
                }
            } // end check right neighbor

            // add interior points
            int point_12_ind = new_vert_ind;
            int point_16_ind = point_12_ind + 3;
            int point_19_ind = point_16_ind + 2;
            for (int ii = 0; ii < 3; ii++) {
                new_xs.push_back(subpanel_xs[1]);
                new_vs.push_back(subpanel_vs[1+ii]);
            }
            for (int ii = 0; ii < 2; ii++) {
                new_xs.push_back(subpanel_xs[2]);
                new_vs.push_back(subpanel_vs[1+2*ii]);
            }
            for (int ii = 0; ii < 3; ii++) {
                new_xs.push_back(subpanel_xs[3]);
                new_vs.push_back(subpanel_vs[1+ii]);
            }
            new_vert_ind += 8;

            // generate new panels
            // add these to list of prospective_panel_indices
            if (do_adaptive_refine) {
                for (int ii = num_new_panels; ii < num_new_panels + 4; ++ii) {
                    prospective_leaf_inds.push_back(ii);
                }
            }
            // panel->child_inds_start = num_new_panels;
            panel->set_child_inds_start(num_new_panels);
            // printf("post refinement, panel looks like:\n");
            // panel->print_panel();
            int child_level = panel->level + 1;
            int panel_ind = panel->panel_ind;
            int* point_inds = panel->point_inds;
            // for (int ii = 0; ii < ; ii++) {
            //     panel_vertex_inds[ii] = panel->point_inds[ii];
            // }
            panels.push_back(Panel {num_new_panels, child_level, panel_ind, 0, 
                    point_inds[0], point_9_ind, point_inds[1],
                    point_11_ind, point_12_ind, point_12_ind + 1,
                    point_inds[3], point_16_ind, point_inds[4],
                    child_0_left_nbr_ind, num_new_panels + 1, 
                    num_new_panels + 2, child_0_bottom_nbr_ind,
                    panel->is_left_bdry, false});
            panels.push_back(Panel {num_new_panels+1, child_level, panel_ind, 1, 
                    point_inds[1], point_10_ind, point_inds[2],
                    point_12_ind+1, point_12_ind+2, point_15_ind,
                    point_inds[4], point_16_ind+1, point_inds[5],
                    child_1_left_nbr_ind, child_1_top_nbr_ind,
                    num_new_panels + 3, num_new_panels,
                    panel->is_left_bdry, false});
            panels.push_back(Panel {num_new_panels+2, child_level, panel_ind, 2, 
                    point_inds[3], point_16_ind, point_inds[4],
                    point_18_ind, point_19_ind, point_19_ind+1,
                    point_inds[6], point_23_ind, point_inds[7],
                    num_new_panels, num_new_panels+3, 
                    child_2_right_nbr_ind, child_2_bottom_nbr_ind,
                    false, panel->is_right_bdry});
            panels.push_back(Panel {num_new_panels+3, child_level, panel_ind, 3,
                    point_inds[4], point_16_ind+1, point_inds[5],
                    point_19_ind+1, point_19_ind+2, point_22_ind,
                    point_inds[7], point_24_ind, point_inds[8],
                    num_new_panels+1, child_3_top_nbr_ind,
                    child_3_right_nbr_ind, num_new_panels+2,
                    false, panel->is_right_bdry});

        } // end if panel is flagged
    } //end for loop through panels

    // set fs
    new_fs.reserve(new_xs.size() );
    for (int ii = 0; ii < new_xs.size(); ++ii) {
        new_fs.push_back( f(new_xs.at(ii), new_vs.at(ii)) );
    }

    for (int ii = 0; ii < new_xs.size(); ++ii) {
        xs.push_back(new_xs[ii]); vs.push_back(new_vs[ii]); fs.push_back(new_fs[ii]);
        // particles.push_back(Particle(new_xs.at(ii), new_vs.at(ii), new_fs.at(ii), 0.0));
    }

    // test 
    if (do_adaptive_refine) { 
        for (int ii = 0; ii < prospective_leaf_inds.size(); ++ii) {
            test_panel(prospective_leaf_inds.at(ii));
        }
    }
}

void AMRStructure::generate_mesh(std::function<double (double,double)> f, 
                                 bool do_adaptive_refine, bool is_initial_step) 
{
    bool verbose=false;
    create_prerefined_mesh();

    if (is_initial_step) {
        for (int ii = 0; ii < xs.size(); ii++) {
            fs[ii] = (*f0)(xs[ii],vs[ii]);
        }
    } else {
        int nx_points = 2*npanels_x + 1;
        int nv_points = 2*npanels_v + 1;

        auto start = high_resolution_clock::now();
        interpolate_to_initial_xvs(fs,xs,vs, nx_points, nv_points,verbose);
        auto stop = high_resolution_clock::now();
        add_time(interp_time, duration_cast<duration<double>>(stop - start) );
    }

    int num_panels_pre_refine = panels.size();

    for (int ii = minimum_unrefined_index; ii < num_panels_pre_refine; ++ii) {
        test_panel(ii);
    }
    while (need_further_refinement) {
        need_further_refinement = false;
        refine_panels(f, do_adaptive_refine);
    }

    set_leaves_weights();

    Q0 = 0;
    for (int ii = 0; ii < q_ws.size(); ii++) {
        Q0 += q_ws[ii];
    }
}

void AMRStructure::test_panel(int panel_ind) {
    double panel_fs[9];
    auto panel_it = panels.begin() + panel_ind;
    for (int ii = 0; ii < 9; ++ii) {
        // panel_fs[ii] = particles.at(panel_it->get_vertex_ind(ii)).get_f();
        panel_fs[ii] = fs[panel_it->point_inds[ii]];
    }
    double max_f = panel_fs[0];
    double min_f = panel_fs[0];
    for (int ii = 1; ii < 9; ++ii) {
        double fii = panel_fs[ii];
        if (max_f < fii) { max_f = fii; }
        if (min_f > fii) { min_f = fii; }
    }

    if (panel_it->level < max_height & max_f - min_f > 100) { 
        panel_it->needs_refinement = true; 
        need_further_refinement = true;
    }
}


void AMRStructure::set_leaves_weights() {
    leaf_inds = std::vector<int> ();
    // for (int ii = 0; ii < particles.size(); ++ii) {
    //     // particles[ii].q_w = 0;
    // }
    q_ws = std::vector<double> (xs.size());
    recursively_set_leaves_weights(0);

    // for (int ii = 0; ii < particles.size(); ++ii) {
    //     particles[ii].q_w *= particles[ii].f;
    // }
    for (int ii = 0; ii < xs.size(); ++ii) {
        q_ws[ii] *= fs[ii];
    }

}

void AMRStructure::recursively_set_leaves_weights(int panel_ind) {
    auto panel_it = panels.begin() + panel_ind;
    if (panel_it->is_refined) {
        int child_start = panel_it->child_inds_start;
        for (int ii = 0; ii < 4; ii++) {
            recursively_set_leaves_weights(child_start + ii);
        }
    }
    else {
        leaf_inds.push_back(panel_ind);
        double dx = xs[panel_it->point_inds[3]] - xs[panel_it->point_inds[0]];
        // double v0 = particles[panel_it->vertex_inds[0]].v;
        // double v1 = particles[panel_it->vertex_inds[1]].v;
        double v0 = vs[panel_it->point_inds[0]];
        double v1 = vs[panel_it->point_inds[1]];
        double dv = v1 - v0;
        switch (quad) {
            case simpsons : {
                double qdxdv9 = q*dx * dv / 9.0;
                double weights[9] = {1.0,4.0,1.0, 4.0, 16.0,4.0,1.0,4.0,1.0};
                for (int ii = 0; ii < 9; ii++) {
                    q_ws[panel_it->point_inds[ii]] += qdxdv9 * weights[ii];
                }
                break;
            }
            default : {// trap 
                double qdxdv4 = q*dx * dv / 4.0;
                // cout << "area factor " << qdxdv4 << endl;
                double weights[9] = {1.0,2.0,1.0,2.0,4.0,2.0,1.0, 2.0, 1.0};
                for (int ii = 0; ii < 9; ii++) {
                    q_ws[panel_it->point_inds[ii]] += qdxdv4 * weights[ii];
                    // q_ws[panel_it->point_inds[ii]] += weights[ii];
                }
                break;
            }
        }
    }
}

void AMRStructure::remesh() {
    // create copy of current panels and particle data
    do_adaptively_refine = false;

    auto start = high_resolution_clock::now();
    old_panels = std::vector<Panel> (); old_panels.reserve(panels.size() );
    old_xs = std::vector<double> (); old_xs.reserve(xs.size());
    old_vs = std::vector<double> (); old_vs.reserve(xs.size());
    old_fs = std::vector<double> (); old_fs.reserve(xs.size());


    for (const auto& panel : panels) {
        old_panels.push_back(Panel (panel));
    }
    for (int ii = 0; ii < xs.size(); ++ii ) {
        old_xs.push_back(xs[ii]);
        old_vs.push_back(vs[ii]);
        old_fs.push_back(fs[ii]);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // cout << "Old data copy time " << duration.count() << " microseconds." << endl << endl;

    bool is_initial_step = false;
    generate_mesh([&] (double x, double v) { return interpolate_from_mesh(x,v,true);} , do_adaptively_refine, is_initial_step);
    init_e();
}
