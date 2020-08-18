#include "AMRStructure.hpp"


int AMRStructure::find_leaf_containing_point_from_neighbor(double& tx, double& tv, int leaf_ind, std::set<int>& history, bool verbose) {
    Panel* panel = &(old_panels[leaf_ind]);
    double x_bl = old_xs[panel->point_inds[0]]; double v_bl = old_vs[panel->point_inds[0]];
    double x_tl = old_xs[panel->point_inds[2]]; double v_tl = old_vs[panel->point_inds[2]];
    double x_br = old_xs[panel->point_inds[6]]; double v_br = old_vs[panel->point_inds[6]];
    double x_tr = old_xs[panel->point_inds[8]]; double v_tr = old_vs[panel->point_inds[8]];
    if (verbose) {
        cout << "(x,v)_bl = (" << x_bl << ", " << v_bl << ")" << endl;
        cout << "(x,v)_tl = (" << x_tl << ", " << v_tl << ")" << endl;
        cout << "(x,v)_br = (" << x_br << ", " << v_br << ")" << endl;
        cout << "(x,v)_tr = (" << x_tr << ", " << v_tr << ")" << endl;
    }

    bool ineq_right = (x_tr - x_br) * (tv - v_br) > (v_tr - v_br) * (tx - (x_br));
    bool ineq_left = (x_tl - x_bl) * (tv - v_bl) <= (v_tl - v_bl) * (tx - x_bl);
    bool ineq_top = (x_tr - x_tl) * (tv - v_tl) < (v_tr - v_tl) * (tx - x_tl);
    bool ineq_bottom = (x_br - x_bl) * (tv - v_bl) >= (v_br - v_bl) * (tx - x_bl);
    int new_leaf_ind = leaf_ind;
    if (verbose) {
        cout << "(tx,tv) = (" << tx << ", " << tv << ")" << endl;
        cout << "panel " << leaf_ind << endl;
        cout << "ineq_right " << ineq_right << endl;
        cout << "(x_tr - x_tl) * (tv - v_tl)" << (x_tr - x_tl) * (tv - v_tl) << endl;
        cout << "(v_tr - v_tl) * (tx - x_tl)" << (v_tr - v_tl) * (tx - x_tl) << endl;
        cout << "ineq_top " << ineq_top << endl;
        cout << "ineq_bottom " << ineq_bottom << endl;
        cout << "ineq_left " << ineq_left << endl;
    }
    if (! ineq_right) {
        if (panel->is_right_bdry) { tx -= Lx; }
        if (verbose) {
            cout << "tx " << tx << endl;
        }

        if (panel->right_nbr_ind < 0) {
            new_leaf_ind = old_panels[panel->parent_ind].right_nbr_ind;
            if (verbose) {
                cout << "in parent right" << endl;
                cout << "next leaf test " << new_leaf_ind << endl;
            }
        } else {
            Panel* panel_right = &old_panels[panel->right_nbr_ind];
            if (panel_right->is_refined) { 
                new_leaf_ind = panel_right->child_inds_start; 
                if (verbose) {
                    cout << "in panel right children" << endl;
                    cout << "next leaf test " << new_leaf_ind << endl;
                }
            }
            else { 
                new_leaf_ind = panel->right_nbr_ind; 
                if (verbose) {
                    cout << "in panel right" << endl;
                    cout << "next leaf test " << new_leaf_ind << endl;
                }
            }
        }
    } else {
        if (!ineq_top && panel->top_nbr_ind != -2) {
            if (panel->top_nbr_ind == -1) {
                new_leaf_ind = old_panels[panel->parent_ind].top_nbr_ind;
                if (verbose) {
                    cout << "in parent top" << endl;
                    cout << "next leaf test " << new_leaf_ind << endl;
                }
            } else {
                Panel* panel_top = &old_panels[panel->top_nbr_ind];
                if (panel_top->is_refined) {
                    new_leaf_ind = panel_top->child_inds_start;
                    if (verbose) {
                        cout << "in top children" << endl;
                        cout << "next leaf test " << new_leaf_ind << endl;
                    }
                } else {
                    new_leaf_ind = panel->top_nbr_ind;
                    if (verbose) {
                        cout << "in panel top" << endl;
                        cout << "next leaf test " << new_leaf_ind << endl;
                    }
                }
            }
        } // end if checking top boundary 
        else {
            if (!ineq_bottom && panel->bottom_nbr_ind != -2) {
                if (panel->bottom_nbr_ind == -1) {
                    new_leaf_ind = old_panels[panel->parent_ind].bottom_nbr_ind;
                    if (verbose) {
                        cout << "in parent bottom" << endl;
                        cout << "next leaf test " << new_leaf_ind << endl;
                    }
                } else {
                    Panel* panel_bottom = &old_panels[panel->bottom_nbr_ind];
                    if (panel_bottom -> is_refined) {
                        new_leaf_ind = panel_bottom->child_inds_start + 1;
                        if (verbose) {
                            cout << "in parent bottom children" << endl;
                            cout << "next leaf test " << new_leaf_ind << endl;
                        }
                    } else {
                        new_leaf_ind = panel->bottom_nbr_ind;
                        if (verbose) {
                            cout << "in parent bottom" << endl;
                            cout << "next leaf test " << new_leaf_ind << endl;
                        }
                    }
                }
            } // end if checking bottom boundary 
            else {
                if (! ineq_left) {
                    if (panel->is_left_bdry) { 
                        tx += Lx; 
                        if (verbose) {
                            cout << "tx, " << tx << endl;
                        }
                    }
                    if (panel->left_nbr_ind < 0) {
                        new_leaf_ind = old_panels[panel->parent_ind].left_nbr_ind;
                        if (verbose) {
                            cout << "in parent left" << endl;
                            cout << "next leaf test " << new_leaf_ind << endl;
                        }
                    } else {
                        Panel* panel_left = &old_panels[panel->left_nbr_ind];
                        if (panel_left->is_refined) {
                            new_leaf_ind = panel_left->child_inds_start+2;
                            if (verbose) {
                                cout << "in right children" << endl;
                                cout << "next leaf test " << new_leaf_ind << endl;
                            }
                        }
                        else {
                            new_leaf_ind = panel->left_nbr_ind;
                            if (verbose) {
                                cout << "in panel left" << endl;
                                cout << "next leaf test " << new_leaf_ind << endl;
                            }
                        }
                    }
                } // end if checking left boundary
            } // end else for left
        } //end else for bottom/left
    } // end else for top/bottom/left
    if (verbose) {
        cout << "History" << endl;
        for (std::set<int>::iterator it = history.begin(); it != history.end(); ++it) {
            cout << *it << " ";
        }
        cout << endl;
    }
    if (history.find(new_leaf_ind) == history.end()) {
        history.emplace(new_leaf_ind);
        if (verbose) {
            cout << "running find_leaf again for (x,v)=(" << tx << ", " << tv << ")" << endl;
        }
        new_leaf_ind = find_leaf_containing_point_from_neighbor(tx,tv,new_leaf_ind, history, verbose);
    }
    return new_leaf_ind;
}

int AMRStructure::find_leaf_containing_xv_recursively(double  &x, const double &v, int panel_ind, bool verbose) {
    int leaf_ind;
    int subpanel_ind;
    int child_inds_start;
    // double x_temp = x;
    
    Panel* panel = &(old_panels[panel_ind]);
    child_inds_start = panel->child_inds_start;

    if (verbose) {
        cout << "In panel " << panel_ind << endl;
    }
    if (! panel->is_refined ) {
        leaf_ind = panel_ind;
        if (verbose) {
            cout << "leaf panel!" << endl;
        }
    } else {
        double x_bl = old_xs[panel->point_inds[0]]; double v_bl = old_vs[panel->point_inds[0]];
        double x_ml = old_xs[panel->point_inds[1]]; double v_ml = old_vs[panel->point_inds[1]];
        double x_tl = old_xs[panel->point_inds[2]]; double v_tl = old_vs[panel->point_inds[2]];
        double x_bm = old_xs[panel->point_inds[3]]; double v_bm = old_vs[panel->point_inds[3]];
        double x_mm = old_xs[panel->point_inds[4]]; double v_mm = old_vs[panel->point_inds[4]];
        double x_tm = old_xs[panel->point_inds[5]]; double v_tm = old_vs[panel->point_inds[5]];
        double x_br = old_xs[panel->point_inds[6]]; double v_br = old_vs[panel->point_inds[6]];
        double x_mr = old_xs[panel->point_inds[7]]; double v_mr = old_vs[panel->point_inds[7]];
        double x_tr = old_xs[panel->point_inds[8]]; double v_tr = old_vs[panel->point_inds[8]];
        // left_shear_ineq = (x_tl - x_bl) * (v - v_bl) > (v_tl - v_bl) * (x_temp - x_bl)

        bool ineq_1_bottom = (x_mm - x_ml) * (v - v_ml) >= (v_mm - v_ml) * (x - x_ml);
        bool ineq_1_right = (x_tm - x_mm) * (v - v_mm) >= (v_tm - v_mm) * (x - x_mm);
        bool ineq_3_bottom = (x_mr - x_mm) * (v - v_mm) >= (v_mr - v_mm) * (x - x_mm);

        if (verbose) {
            std::cout << "ineq_1_bottom = " << ineq_1_bottom << endl;
            std::cout << "ineq_1_right = " << ineq_1_right << endl;
            std::cout << "ineq_3_bottom = " << ineq_3_bottom << endl;
        }

        if (ineq_1_bottom && ineq_1_right) {

            bool ineq_1_top = (x_tm - x_tl) * (v - v_tl) <= (v_tm - v_tl) * (x - x_tl);
            Panel* child_1 = &old_panels[child_inds_start+1];
            int child_1_top_nbr_ind = child_1->top_nbr_ind;
            if (ineq_1_top ||  child_1_top_nbr_ind < 0) {
                bool ineq_1_left = (x_tl - x_ml) * (v - v_ml) <= (v_tl - v_ml) * (x - x_ml);
                int child_1_left_nbr_ind = child_1->left_nbr_ind;
                if (ineq_1_left || child_1_left_nbr_ind < 0) {
                    subpanel_ind = child_inds_start + 1;
                    if (verbose) {
                        cout << "in child 1, panel " << subpanel_ind << endl;
                    }
                }
                else {
                    subpanel_ind = child_1_left_nbr_ind;
                    if (panel->is_left_bdry) {
                        x += Lx;
                    }
                }
            } else {
                subpanel_ind = child_1_top_nbr_ind;
            }
            
        } else if (ineq_3_bottom && !ineq_1_right)
        {
            bool ineq_3_top = (x_tr - x_tm) * (v - v_tm) <= (v_tr - v_tm) * (x - x_tm);
            Panel* child_3 = &old_panels[child_inds_start+3];
            int child_3_top_nbr_ind = child_3->top_nbr_ind;
            if (ineq_3_top || child_3_top_nbr_ind < 0) {
                bool ineq_3_right = (x_tr - x_mr) * (v - v_mr) <= (v_tr - v_mr) * (x - x_mr);
                int child_3_right_nbr_ind = child_3->right_nbr_ind;
                if (!ineq_3_right || child_3_right_nbr_ind < 0) {
                    subpanel_ind = child_inds_start + 3;
                    if (verbose) {
                        cout << "in child 3, panel " << subpanel_ind << endl;
                    }
                }
                else {
                    subpanel_ind = child_3_right_nbr_ind;
                    if (panel->is_right_bdry) {
                        x -= Lx;
                    }
                }
            } else {
                subpanel_ind = child_3_top_nbr_ind;
            }
        } else {
            bool ineq_0_right = (x_mm - x_bm) * (v - v_bm) >= (v_mm - v_bm) * (x - x_bm);
            if (verbose) {
                std::cout << "ineq_0_right = " << ineq_0_right << endl;
            }
            if (ineq_0_right) {
                bool ineq_0_bottom = (x_bm - x_bl) * (v - v_bl) >= (v_bm - v_bl) * (x - x_bl);
                Panel* child_0 = &old_panels[child_inds_start];
                int child_0_bottom_nbr_ind = child_0->bottom_nbr_ind;
                if (ineq_0_bottom || child_0_bottom_nbr_ind < 0) {
                    bool ineq_0_left = (x_ml - x_bl) * (v - v_bl) <= (v_ml - v_bl) * (x - x_bl);
                    int child_0_left_nbr_ind = child_0->left_nbr_ind;
                    if (ineq_0_left || child_0_left_nbr_ind < 0) {
                        
                        subpanel_ind = child_inds_start;
                        if (verbose) {
                            cout << "in child 0, panel " << subpanel_ind << endl;
                        }
                    }
                    else {
                        subpanel_ind = child_0_left_nbr_ind;
                        if (panel->is_left_bdry) { x += Lx; }
                    }
                } else {
                    subpanel_ind = child_0_bottom_nbr_ind;
                }
            } else
            {
                bool ineq_2_bottom = (x_br - x_bm) * (v - v_bm) >= (v_br - v_bm) * (x - x_bm);
                Panel* child_2 = &old_panels[child_inds_start+2];
                int child_2_bottom_nbr_ind = child_2->bottom_nbr_ind;
                if (ineq_2_bottom || child_2_bottom_nbr_ind < 0) {
                    bool ineq_2_right = (x_mr - x_br) * (v - v_br) <= (v_mr - v_br) * (x - x_br);
                    int child_2_right_nbr_ind = child_2->right_nbr_ind;
                    if (! ineq_2_right || child_2_right_nbr_ind < 0) {
                        subpanel_ind = child_inds_start + 2;
                        if (verbose) {
                            cout << "in child 2, panel " << subpanel_ind << endl;
                        }
                    }
                    else {
                        subpanel_ind = child_2_right_nbr_ind;
                        if (panel->is_right_bdry) {
                            x -= Lx;
                        }
                    }
                } else {
                    subpanel_ind = child_2_bottom_nbr_ind;
                }
            }
            
        }

        leaf_ind = find_leaf_containing_xv_recursively(x, v, subpanel_ind, verbose);

    }
    return leaf_ind;
}


void AMRStructure::shift_xs(std::vector<double>& shifted_xs, const std::vector<double>& xs, const std::vector<double>& vs) {
    bool verbose = false;
    // shifted_xs = std::vector<double> (xs.size() );
    double x_bl, x_tl, x_br, x_tr;
    double v_bl, v_tl, v_br, v_tr;
    x_bl = this->old_xs[0]; v_bl = this->old_vs[0];
    x_tl = this->old_xs[2]; v_tl = this->old_vs[2];
    x_br = this->old_xs[6]; v_br = this->old_vs[6];
    x_tr = this->old_xs[8]; v_tr = this->old_vs[8];
    if (verbose) {
        cout << "(x,v)_bl (" << x_bl << ", " << v_bl << ")" << endl;
        cout << "(x,v)_tl (" << x_tl << ", " << v_tl << ")" << endl;
        cout << "(x,v)_br (" << x_br << ", " << v_br << ")" << endl;
        cout << "(x,v)_tr (" << x_tr << ", " << v_tr << ")" << endl;
    }

    for (int ii = 0; ii < xs.size(); ++ii) {
        double v=vs[ii];
        double x_temp = xs[ii];
        bool ineq_00_left = (x_tl - x_bl) * (v - v_bl) <= (v_tl - v_bl) * (x_temp - x_bl);

        if (verbose) {
            cout << "(x,v)= (" << x_temp << ", " << v << ")" << endl;
            cout << "ineq_00_left, " << ineq_00_left << endl; 
        }

        while (not ineq_00_left) {
            x_temp += Lx;
            ineq_00_left = (x_tl - x_bl) * (v - v_bl) <= (v_tl - v_bl) * (x_temp - x_bl);
            if(verbose) {
                cout << "post shift x= (" << x_temp << ", ineq_00_left, " << ineq_00_left << endl; 
            }
        }
        bool ineq_00_right = (x_tr - x_br) * (v - v_br) > (v_tr - v_br) * (x_temp - (x_br));
        if (verbose) {
            cout << "ineq_00_right, " << ineq_00_left << endl; 
        }
        while (not ineq_00_right) {
            x_temp -= Lx;
            ineq_00_right = (x_tr - x_br) * (v - v_br) > (v_tr - v_br) * (x_temp - (x_br));
            if(verbose) {
                cout << "post shift x= (" << x_temp << ", ineq_00_right, " << ineq_00_left << endl; 
            }
        }
        shifted_xs[ii] = x_temp;
    }
}

void AMRStructure::interpolate_to_initial_xvs(
    std::vector<double>& fs, std::vector<double>& xs, std::vector<double>& vs, 
    int nx, int nv, bool verbose) 
{
    bool print_profile = false;

    auto start = high_resolution_clock::now();
    // bool verbose = false;
    std::vector<double> shifted_xs(xs.size() );
    shift_xs(shifted_xs, xs, vs);
    if (verbose) {
        cout << "shift xs" << endl;
        std::copy(shifted_xs.begin(), shifted_xs.end(), std::ostream_iterator<double>(cout, ", "));
        cout << endl;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    if (print_profile) {
        cout << "shift time " << duration.count() << " microseconds" << endl;
    }

    // have to sort points:
    start = high_resolution_clock::now();

    std::vector<int> sort_indices(xs.size());
    for (int ii = 0; ii < xs.size(); ii++ ) { sort_indices[ii] = ii; }
    bool do_sort = false;
    if (do_sort) {
        std::sort(sort_indices.begin(), sort_indices.end(),
            [&] (int a, int b) 
            { 
                if (vs[a] != vs[b]) { return vs[a] < vs[b]; }
                else {
                    return shifted_xs[a] < shifted_xs[b];
                }
            });
    }
    std::vector<double> sortxs(shifted_xs.size()), sortvs(vs.size());
    for (int ii = 0; ii < xs.size(); ii++) {
        sortxs[ii] = shifted_xs[sort_indices[ii]];
        sortvs[ii] = vs[sort_indices[ii]];
    }
    if (verbose) {
        cout << "sorted xs: " << endl;
        std::copy(sortxs.begin(), sortxs.end(), std::ostream_iterator<double>(cout, ", "));
        cout << endl << "sorted vs: " << endl;
        std::copy(sortvs.begin(), sortvs.end(), std::ostream_iterator<double>(cout, ", "));
        cout << endl;
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    if (print_profile) {
        cout << "sort time " << duration.count() << " microseconds" << endl;
    }
    
    start = high_resolution_clock::now();
    std::vector<double> sortfs(xs.size() );

    std::vector<int> leaf_panel_of_points(xs.size() );
    std::vector<std::vector<int> > point_in_leaf_panels_by_inds(old_panels.size() );

    int leaf_ind = find_leaf_containing_xv_recursively(sortxs[0], sortvs[0], 0, verbose);
    std::vector<int> first_column_leaf_inds(nv);

    if (verbose) {
        std::cout << "nx x nv= " << nx << " x " << nv << endl;
        cout << "xs size: " << xs.size() << endl; 
    }

    for (int ii =0; ii < nv; ++ii) {
        int point_ind = ii * nv;
        std::set<int> history;
        history.emplace(leaf_ind);
        // cout << "testing point " << point_ind << ", x= " << sortxs[point_ind] << ", v= " << sortvs[point_ind] << endl;
        leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortvs[point_ind], leaf_ind, history, verbose);
        first_column_leaf_inds[ii] = leaf_ind;
        point_in_leaf_panels_by_inds[leaf_ind].push_back(point_ind);
        leaf_panel_of_points[sort_indices[point_ind]] = leaf_ind;


        // cout << "point " << sort_indices[point_ind] << " in panel " << leaf_ind << " (sorted ind " << point_ind << ")" << endl;
    }

#pragma omp parallel
{
    // int num_threads = omp_get_num_threads();
    // if (omp_get_thread_num() == 0) {
    //     cout << "Number of threads from calc_E: " << num_threads << endl;
    // }
    #pragma omp for
    for (int ii = 0; ii < nv; ++ii) {
        int point_ind = ii * nv;
        int leaf_ind_c = first_column_leaf_inds[ii];
        for (int jj = 1; jj < nx; ++jj) {
            point_ind++;
            // cout << "Testing point " << point_ind << ", (x,v)= (" << sortxs[point_ind] << ", " << sortvs[point_ind] << ")" << endl;
            std::set<int> history;
            history.emplace(leaf_ind_c);
            leaf_ind_c = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortvs[point_ind], leaf_ind_c, history, verbose);
            // point_in_leaf_panels_by_inds[leaf_ind_c].push_back(point_ind);
            leaf_panel_of_points[sort_indices[point_ind]] = leaf_ind_c;

            // cout << "sorted ind " << point_ind << ", ind " << sort_indices[point_ind] << " is in panel " << leaf_ind << endl;
        }
    }    
} // end omp parallel

    for (int ii = 0; ii < leaf_panel_of_points.size(); ++ii) {
        point_in_leaf_panels_by_inds[leaf_panel_of_points[ii]].emplace_back(ii);
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    if (print_profile) {
        cout << "find leaf time " << duration.count() << " microseconds" << endl;
    }
    
    start = high_resolution_clock::now();
    for (int panel_ind = 0; panel_ind < old_panels.size(); panel_ind++) {
        if (point_in_leaf_panels_by_inds[panel_ind].size() > 0) {
            interpolate_from_panel_to_points(sortfs,sortxs,sortvs,point_in_leaf_panels_by_inds[panel_ind], panel_ind);
        }
    }
    for (int ii = 0; ii < fs.size(); ii++) {
        fs[sort_indices[ii]] = sortfs[ii];
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    if (print_profile) {
        cout << "evaluate interpolant time " << duration.count() << " microseconds" << endl;
    }
}

void AMRStructure::interpolate_from_panel_to_points(
    std::vector<double>& values, std::vector<double>& xs, std::vector<double>& vs,
    std::vector<int>& point_inds, int panel_ind) 
{
    bool verbose = false;
    Panel* panel = &(old_panels[panel_ind]);
    const int* panel_point_inds = panel->point_inds;
    double panel_xs[9], panel_vs[9], panel_fs[9];

    if (verbose) {
        cout << "Getting panel points" << endl;
    }
    for (int ii = 0; ii < 9; ++ii) {
        // Particle* part = &particles[vertex_inds[ii]];
        int pind = panel_point_inds[ii];
        // panel_xs[ii] = part->get_x();
        // panel_vs[ii] = part->get_v();
        // panel_fs[ii] = part->get_f();
        panel_xs[ii] = old_xs[pind];
        panel_vs[ii] = old_vs[pind];
        panel_fs[ii] = old_fs[pind];
    }

    if (verbose) {
        std::cout << "Interpolating from " << std::endl;
        for (int ii = 0; ii < 9; ++ii) {
            std::cout << "(x,v,f)=(" << panel_xs[ii] << ", " << panel_vs[ii] << ", " << panel_fs[ii] << ")" << std::endl;
        }
    
        std::cout << "onto";
        for (int ii = 0; ii < point_inds.size(); ii++) {
            int pind = point_inds[ii];
            std::cout << "(" << xs[pind] << ", " << vs[pind] << ")\n";
        }
    }

    // double dx, dv, 
    double panel_dx[9], panel_dv[9];
    // dx = x - panel_xs[4];
    // dv = v - panel_vs[4];
    std::vector<double> dxs(point_inds.size()), dvs(point_inds.size());
    for (int ii = 0; ii < point_inds.size(); ++ii) {
        int pind = point_inds[ii];
        dxs[ii] = xs[pind] - panel_xs[4];
        dvs[ii] = vs[pind] - panel_vs[4];
    }
    for (int ii = 0; ii < 9; ii ++) {
        panel_dx[ii] = panel_xs[ii] - panel_xs[4];
        panel_dv[ii] = panel_vs[ii] - panel_vs[4];
    }

    if (verbose) {
        std::cout << "test point distance from midpoint:" << std::endl;
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            std::cout << ii <<": dx=" << dxs[ii] <<", dv=" << dvs[ii] << std::endl;
        } 
        std::cout << "panel vertex distances from midpoint:" << std::endl;
        for (int ii = 0; ii < 9; ii++ ) {
            std::cout << ii << ": " << panel_dx[ii] << ", " << panel_dv[ii] << std::endl;
        }
    }

    Eigen::Matrix<double,9,9> A;
    for (int ii = 0; ii < 9; ++ii) {
        A(ii,0) = 1; A(ii,1) = panel_dx[ii];
        A(ii,2) = panel_dx[ii] * panel_dv[ii];
        A(ii,3) = panel_dv[ii];
        A(ii,4) = panel_dx[ii] * panel_dx[ii];
        A(ii,5) =  panel_dx[ii] * panel_dx[ii] * panel_dv[ii];
        A(ii,6) = panel_dx[ii] * panel_dx[ii] * panel_dv[ii] * panel_dv[ii];
        A(ii,7) = panel_dx[ii] * panel_dv[ii] * panel_dv[ii];
        A(ii,8) = panel_dv[ii] * panel_dv[ii];
    }
    Eigen::Map<Eigen::Matrix<double,9,1>> b(panel_fs);
    Eigen::Matrix<double,9,1> c = A.lu().solve(b);


    if (verbose) {
        std::cout << "Here is the matrix A:\n" << A << std::endl;
        std::cout << "Here is the f vector b:\n" << b << std::endl;
        std::cout << "Here is the coefficient vector c:\n" << c << std::endl;
    }
    Eigen::Matrix<double, Dynamic, Dynamic> Dx(point_inds.size(),9);
    for (int ii = 0; ii < point_inds.size(); ++ii) {
        double dxsq = dxs[ii] * dxs[ii];
        double dvsq = dvs[ii] * dvs[ii];
        Dx(ii,0) = 1; Dx(ii,1) = dxs[ii];
        Dx(ii,2) = dxs[ii] * dvs[ii];
        Dx(ii,3) = dvs[ii];
        Dx(ii,4) = dxsq;
        Dx(ii,5) =  dxsq * dvs[ii];
        Dx(ii,6) = dxsq * dvsq;
        Dx(ii,7) = dxs[ii] * dvsq;
        Dx(ii,8) = dvsq;
    }
    Eigen::Matrix<double, Dynamic,1> interp_vals = Dx * c;
    for (int ii = 0; ii < point_inds.size(); ++ii) {
        values[point_inds[ii]] = interp_vals(ii);
    }

    // return c(0) + c(1)*dx + c(2) * dx*dv + c(3) * dv +
    //         c(4) * dx*dx + c(5) * dx*dx*dv + c(6) * dx*dx*dv*dv +
    //         c(7) * dx*dv*dv + c(8) * dv*dv;

}

double AMRStructure::interpolate_from_panel(double x, double v, int panel_ind, bool verbose) {
    
    Panel* panel = &(old_panels[panel_ind]);
    const int* point_inds = panel->point_inds;
    double panel_xs[9], panel_vs[9], panel_fs[9];

    for (int ii = 0; ii < 9; ++ii) {
        // Particle* part = &particles[vertex_inds[ii]];
        int vind = point_inds[ii];
        // panel_xs[ii] = part->get_x();
        // panel_vs[ii] = part->get_v();
        // panel_fs[ii] = part->get_f();
        panel_xs[ii] = old_xs[vind];
        panel_vs[ii] = old_vs[vind];
        panel_fs[ii] = old_fs[vind];
    }

    if (verbose) {
        std::cout << "Interpolating from " << std::endl;
        for (int ii = 0; ii < 9; ++ii) {
            std::cout << "(x,v,f)=(" << panel_xs[ii] << ", " << panel_vs[ii] << ", " << panel_fs[ii] << ")" << std::endl;
        }
    
        std::cout << "onto (" << x << ", " << v << ")\n";
    }

    double dx, dv, panel_dx[9], panel_dv[9];
    dx = x - panel_xs[4];
    dv = v - panel_vs[4];
    for (int ii = 0; ii < 9; ii ++) {
        panel_dx[ii] = panel_xs[ii] - panel_xs[4];
        panel_dv[ii] = panel_vs[ii] - panel_vs[4];
    }

    if (verbose) {
        std::cout << "test point distance from midpoint: dx=" << dx <<", dv=" << dv << std::endl;
        std::cout << "panel vertex distances from midpoint:" << std::endl;
        for (int ii = 0; ii < 9; ii++ ) {
            std::cout << ii << ": " << panel_dx[ii] << ", " << panel_dv[ii] << std::endl;
        }
    }

    Eigen::Matrix<double,9,9> A;
    for (int ii = 0; ii < 9; ++ii) {
        A(ii,0) = 1; A(ii,1) = panel_dx[ii];
        A(ii,2) = panel_dx[ii] * panel_dv[ii];
        A(ii,3) = panel_dv[ii];
        A(ii,4) = panel_dx[ii] * panel_dx[ii];
        A(ii,5) =  panel_dx[ii] * panel_dx[ii] * panel_dv[ii];
        A(ii,6) = panel_dx[ii] * panel_dx[ii] * panel_dv[ii] * panel_dv[ii];
        A(ii,7) = panel_dx[ii] * panel_dv[ii] * panel_dv[ii];
        A(ii,8) = panel_dv[ii] * panel_dv[ii];
    }
    Eigen::Map<Eigen::Matrix<double,9,1>> b(panel_fs);
    Eigen::Matrix<double,9,1> c = A.lu().solve(b);


    if (verbose) {
        std::cout << "Here is the matrix A:\n" << A << std::endl;
        std::cout << "Here is the f vector b:\n" << b << std::endl;
        std::cout << "Here is the coefficient vector c:\n" << c << std::endl;
    }

    return c(0) + c(1)*dx + c(2) * dx*dv + c(3) * dv +
            c(4) * dx*dx + c(5) * dx*dx*dv + c(6) * dx*dx*dv*dv +
            c(7) * dx*dv*dv + c(8) * dv*dv;
}

double AMRStructure::interpolate_from_mesh(double x, double v, bool verbose) {
    int leaf_containing = find_leaf_containing_xv_recursively(x,v,0, verbose);
    return interpolate_from_panel(x,v,leaf_containing, verbose);
}

void AMRStructure::interpolate_from_mesh(std::vector<double>& values, std::vector<double>& xs, std::vector<double>& vs, bool verbose) {
    std::vector<double> shifted_xs;
    shift_xs(shifted_xs, xs, vs);

    std::vector<int> leaves;
    std::vector<std::vector<int> > point_in_leaf_panels_by_inds(old_panels.size() );
    for (int ii = 0; ii < xs.size(); ++ii) {
        int leaf_ind = find_leaf_containing_xv_recursively(shifted_xs[ii], vs[ii], 0, verbose);
        leaves.push_back(leaf_ind);
        point_in_leaf_panels_by_inds[leaf_ind].push_back(ii);
    }

    for (int panel_ind = 0; panel_ind < old_panels.size(); panel_ind++) {
        if (point_in_leaf_panels_by_inds[panel_ind].size() > 0) {
            interpolate_from_panel_to_points(values,shifted_xs,vs,point_in_leaf_panels_by_inds[panel_ind], panel_ind);
        }
    }

    // cout << endl << "interpolated fs: ";
    // std::copy(values.begin(), values.end(), std::ostream_iterator<double>(cout, " "));
}

void AMRStructure::interpolate_from_mesh_slow(std::vector<double>& values, std::vector<double>& xs, std::vector<double>& vs, bool verbose) {
    std::vector<double> shifted_xs;
    shift_xs(shifted_xs, xs, vs);

    std::vector<int> leaves;
    // std::vector<double> values;
    values = std::vector<double> ();
    for (int ii = 0; ii < xs.size(); ++ii) {
        leaves.push_back(find_leaf_containing_xv_recursively(shifted_xs[ii], vs[ii], 0, verbose));
        values.push_back(interpolate_from_panel(shifted_xs[ii], vs[ii], leaves[ii], verbose));
    }
    // cout << endl << "interpolated fs: ";
    // std::copy(values.begin(), values.end(), std::ostream_iterator<double>(cout, " "));
}

// double AMRStructure::InterpolateDistribution::operator() (double x, double v, bool verbose) {
//     // int leaf_containing = ::find_leaf_containing_xv_recursively(x,v,0, verbose);
//     // return ::interpolate_from_panel(x,v,leaf_containing, verbose);
//     this->find_leaf_containing_xv_recursively(x,v,0, verbose);
//     return 1.0;
// }