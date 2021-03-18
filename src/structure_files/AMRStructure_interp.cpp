#include "AMRStructure.hpp"

// #ifndef DEBUG
// #define DEBUG
// #define DEBUG_L2
// #define DEBUG_L3

/* Changes for refine_p
find leav recursively : if a panel is only refined in p, then child2 = child0, child child3 = child1
*/

int AMRStructure::find_leaf_containing_point_from_neighbor(double& tx, double& tp, bool& beyond_boundary, int leaf_ind, std::set<int>& history, bool verbose) {


    // trouble with interpolation and amr.  What?
#ifdef DEBUG
// if (iter_num >= 9) {
//     verbose = true;
// }
// if (fabs(tx + 0.0314) < 0.0003 && fabs(tp - 0.0122) < 0.0003) {
//     verbose = true;
// }
// verbose=true;
#endif
    // verbose = true;
    // end trouble shooting verbosity change

    if (leaf_ind == 0) {
        history.emplace(leaf_ind);
        cout << "If you see this then you are sending in a leaf ind 0 somewhere you hoped not to." << endl;
        return leaf_ind;
    } else {
        Panel* panel = &(old_panels[leaf_ind]);
        double x_bl = old_xs[panel->point_inds[0]]; double p_bl = old_ps[panel->point_inds[0]];
        double x_tl = old_xs[panel->point_inds[2]]; double p_tl = old_ps[panel->point_inds[2]];
        double x_mid = old_xs[panel->point_inds[4]];
        double x_br = old_xs[panel->point_inds[6]]; double p_br = old_ps[panel->point_inds[6]];
        double x_tr = old_xs[panel->point_inds[8]]; double p_tr = old_ps[panel->point_inds[8]];
        if (verbose) {
            cout << "(x,p)_bl = (" << x_bl << ", " << p_bl << ")" << endl;
            cout << "(x,p)_tl = (" << x_tl << ", " << p_tl << ")" << endl;
            cout << "(x,p)_br = (" << x_br << ", " << p_br << ")" << endl;
            cout << "(x,p)_tr = (" << x_tr << ", " << p_tr << ")" << endl;
        }
        // need to correct periodic distance
        if (tx - x_mid >= Lx/2 && bcs == periodic_bcs) { 
            tx -= Lx; 
            if (verbose) {
                cout << "shifting across boundary, tx= " << tx << endl;
            }
        }
        if (tx - x_mid < -Lx/2 && bcs == periodic_bcs) { 
            tx += Lx; 
            if (verbose) {
                cout << "shifting across boundary, tx= " << tx << endl;
            }
        }


        bool ineq_right = (x_tr - x_br) * (tp - p_br) > (p_tr - p_br) * (tx - (x_br));
        bool ineq_left = (x_tl - x_bl) * (tp - p_bl) <= (p_tl - p_bl) * (tx - x_bl);
        bool ineq_top = (x_tr - x_tl) * (tp - p_tl) < (p_tr - p_tl) * (tx - x_tl);
        bool ineq_bottom = (x_br - x_bl) * (tp - p_bl) >= (p_br - p_bl) * (tx - x_bl);
        int new_leaf_ind = leaf_ind;
        if (verbose) {
            cout << "(tx,tp) = (" << tx << ", " << tp << ")" << endl;
            cout << "testing leaf panel " << leaf_ind << " for containment" << endl;
            cout << "ineq_right " << ineq_right << endl;
            cout << "(x_tr - x_tl) * (tp - p_tl)" << (x_tr - x_tl) * (tp - p_tl) << endl;
            cout << "(p_tr - p_tl) * (tx - x_tl)" << (p_tr - p_tl) * (tx - x_tl) << endl;
            cout << "ineq_top " << ineq_top << endl;
            cout << "ineq_bottom " << ineq_bottom << endl;
            cout << "ineq_left " << ineq_left << endl;
        }
        if (! ineq_right && panel->right_nbr_ind != -2) {
            if (panel->is_right_bdry && bcs == periodic_bcs) { 
                tx -= Lx; 
                if (verbose) {
                    cout << "shifting across boundary, tx= " << tx << endl;
                }
            }

            if (panel->right_nbr_ind == -1) {
                new_leaf_ind = old_panels[panel->parent_ind].right_nbr_ind;
                if (verbose) {
                    cout << "in parent right" << endl;
                    cout << "next leaf test " << new_leaf_ind << endl;
                }
            } else {
                Panel* panel_right = &old_panels[panel->right_nbr_ind];
                if (panel_right->is_refined_xp) { 
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
                    if (panel_top->is_refined_xp) {
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
                        if (panel_bottom -> is_refined_xp) {
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
                    if (! ineq_left && panel->left_nbr_ind != -2) {
                        if (panel->is_left_bdry && bcs == periodic_bcs) { 
                            tx += Lx; 
                            if (verbose) {
                                cout << "shifting across boundary, now tx= " << tx << endl;
                            }
                        }
                        if (panel->left_nbr_ind == -1) {
#ifdef DEBUG
cout << "panel left: " << panel->left_nbr_ind << endl;
cout <<"length of panels_list " << old_panels.size() << endl;
#endif
                            new_leaf_ind = old_panels[panel->parent_ind].left_nbr_ind;
                            if (verbose) {
                                cout << "in parent left" << endl;
                                cout << "next leaf test " << new_leaf_ind << endl;
                            }
                        } else {
                            Panel* panel_left = &old_panels[panel->left_nbr_ind];
                            if (panel_left->is_refined_xp) {
                                new_leaf_ind = panel_left->child_inds_start+2;
                                if (verbose) {
                                    cout << "in left children" << endl;
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
                cout << "running find_leaf again for (x,p)=(" << tx << ", " << tp << ")" << endl;
            }
            new_leaf_ind = find_leaf_containing_point_from_neighbor(tx,tp,beyond_boundary, new_leaf_ind, history, verbose);
        } else if (verbose)
        {
            cout << "done searching." << endl;
        }
        if (!allow_boundary_extrapolation) {
            bool boundary_extrapolating_right = !ineq_right && panel->right_nbr_ind==-2;
            bool boundary_extrapolating_left = !ineq_left && panel->left_nbr_ind==-2;
            bool boundary_extrapolating_top = !ineq_top && panel->top_nbr_ind==-2;
            bool boundary_extrapolating_bottom = !ineq_bottom && panel->bottom_nbr_ind==-2;
            if (boundary_extrapolating_left || boundary_extrapolating_right || 
                boundary_extrapolating_top || boundary_extrapolating_bottom) {
                if (verbose) {
                    cout << "Not allowing boundary extrapolation and this point is flagged as beyond the boundary" << endl;
                }
                beyond_boundary = true;
            }
        }
        if (verbose) {
            cout << "Point (" << tx << ", " << tp << ") is in panel " << new_leaf_ind << endl;
        }
        return new_leaf_ind;
    }
}

int AMRStructure::find_leaf_containing_xp_recursively(double  &x, const double &p, bool& beyond_boundary, int panel_ind, bool verbose) {
    int leaf_ind;
    int subpanel_ind;
    int child_inds_start;

    #ifdef DEBUG
    verbose = true;
    #endif
    // double x_temp = x;
    // if ( fabs(x +1.57) < 0.5 && fabs(p +4.125) < 0.5) { verbose = true; }
    // else {verbose = false; }

    //trouble shooting
    
    Panel* panel = &(old_panels[panel_ind]);
    child_inds_start = panel->child_inds_start;

    if (verbose) {
        cout << "In panel " << panel_ind << endl;
        cout << *panel << endl;
        cout << "testing (x,p)=(" << x << ", " << p << ")" << endl;
    }
    if (! (panel->is_refined_xp || panel->is_refined_p ) ) {
        leaf_ind = panel_ind;
        if (verbose) {
            cout << "leaf panel!" << endl;
        }

        if (!allow_boundary_extrapolation) {
            double x_bl = old_xs[panel->point_inds[0]]; double p_bl = old_ps[panel->point_inds[0]];
            double x_tl = old_xs[panel->point_inds[2]]; double p_tl = old_ps[panel->point_inds[2]];
            double x_br = old_xs[panel->point_inds[6]]; double p_br = old_ps[panel->point_inds[6]];
            double x_tr = old_xs[panel->point_inds[8]]; double p_tr = old_ps[panel->point_inds[8]];
            bool ineq_right = (x_tr - x_br) * (p - p_br) > (p_tr - p_br) * (x - (x_br));
            bool ineq_left = (x_tl - x_bl) * (p - p_bl) <= (p_tl - p_bl) * (x - x_bl);
            bool ineq_top = (x_tr - x_tl) * (p - p_tl) < (p_tr - p_tl) * (x - x_tl);
            bool ineq_bottom = (x_br - x_bl) * (p - p_bl) >= (p_br - p_bl) * (x - x_bl);
            bool boundary_extrapolating_right = !ineq_right && panel->right_nbr_ind==-2;
            bool boundary_extrapolating_left = !ineq_left && panel->left_nbr_ind==-2;
            bool boundary_extrapolating_top = !ineq_top && panel->top_nbr_ind==-2;
            bool boundary_extrapolating_bottom = !ineq_bottom && panel->bottom_nbr_ind==-2;
            if (boundary_extrapolating_left || boundary_extrapolating_right || 
                boundary_extrapolating_top || boundary_extrapolating_bottom) {
                if (verbose) {
                    cout << "Not allowing boundary extrapolation and this point is flagged for beyond the boundary" << endl;
                }
                beyond_boundary = true;
            }
        }
    } else {
        double x_bl = old_xs[panel->point_inds[0]]; double p_bl = old_ps[panel->point_inds[0]];
        double x_ml = old_xs[panel->point_inds[1]]; double p_ml = old_ps[panel->point_inds[1]];
        double x_tl = old_xs[panel->point_inds[2]]; double p_tl = old_ps[panel->point_inds[2]];
        double x_bm = old_xs[panel->point_inds[3]]; double p_bm = old_ps[panel->point_inds[3]];
        double x_mm = old_xs[panel->point_inds[4]]; double p_mm = old_ps[panel->point_inds[4]];
        double x_tm = old_xs[panel->point_inds[5]]; double p_tm = old_ps[panel->point_inds[5]];
        double x_br = old_xs[panel->point_inds[6]]; double p_br = old_ps[panel->point_inds[6]];
        double x_mr = old_xs[panel->point_inds[7]]; double p_mr = old_ps[panel->point_inds[7]];
        double x_tr = old_xs[panel->point_inds[8]]; double p_tr = old_ps[panel->point_inds[8]];
        // left_shear_ineq = (x_tl - x_bl) * (p - p_bl) > (p_tl - p_bl) * (x_temp - x_bl)

        bool ineq_1_bottom = (x_mm - x_ml) * (p - p_ml) >= (p_mm - p_ml) * (x - x_ml);
        bool ineq_1_right = (x_tm - x_mm) * (p - p_mm) >= (p_tm - p_mm) * (x - x_mm);
        bool ineq_3_bottom = (x_mr - x_mm) * (p - p_mm) >= (p_mr - p_mm) * (x - x_mm);

        if (verbose) {
            std::cout << "ineq_1_bottom = " << ineq_1_bottom << endl;
            std::cout << "ineq_1_right = " << ineq_1_right << endl;
            std::cout << "ineq_3_bottom = " << ineq_3_bottom << endl;
        }

        if (ineq_1_bottom && ineq_1_right) {

            bool ineq_1_top = (x_tm - x_tl) * (p - p_tl) <= (p_tm - p_tl) * (x - x_tl);
            Panel* child_1 = &old_panels[child_inds_start+1];
            int child_1_top_nbr_ind = child_1->top_nbr_ind;
            if (ineq_1_top ||  child_1_top_nbr_ind < 0) {
                bool ineq_1_left = (x_tl - x_ml) * (p - p_ml) <= (p_tl - p_ml) * (x - x_ml);
                int child_1_left_nbr_ind = child_1->left_nbr_ind;
                if (ineq_1_left || child_1_left_nbr_ind < 0) {
                    subpanel_ind = child_inds_start + 1;
                    if (verbose) {
                        cout << "in child 1, panel " << subpanel_ind << endl;
                    }
                }
                else {
                    subpanel_ind = child_1_left_nbr_ind;
                    if (panel->is_left_bdry and bcs==periodic_bcs) {
                        x += Lx;
                    }
                }
            } else {
                subpanel_ind = child_1_top_nbr_ind;
            }
            
        } else if (ineq_3_bottom && !ineq_1_right)
        {
            bool ineq_3_top = (x_tr - x_tm) * (p - p_tm) <= (p_tr - p_tm) * (x - x_tm);
            Panel* child_3;
            if (panel->is_refined_p) {
                child_3 = &old_panels[child_inds_start +1];
            } else { // panel is refined in xv
                child_3 = &old_panels[child_inds_start+3];
            }
            int child_3_top_nbr_ind = child_3->top_nbr_ind;
            if (ineq_3_top || child_3_top_nbr_ind < 0) {
                bool ineq_3_right = (x_tr - x_mr) * (p - p_mr) <= (p_tr - p_mr) * (x - x_mr);
                int child_3_right_nbr_ind = child_3->right_nbr_ind;
                if (!ineq_3_right || child_3_right_nbr_ind < 0) {
                    if (panel->is_refined_p) {
                        subpanel_ind = child_inds_start + 1;
                    }
                    else {
                        subpanel_ind = child_inds_start + 3;
                    }
                    if (verbose) {
                        cout << "in child 3, panel " << subpanel_ind << endl;
                    }
                }
                else {
                    subpanel_ind = child_3_right_nbr_ind;
                    if (panel->is_right_bdry && bcs==periodic_bcs) {
                        x -= Lx;
                    }
                }
            } else {
                subpanel_ind = child_3_top_nbr_ind;
            }
        } else {
            bool ineq_0_right = (x_mm - x_bm) * (p - p_bm) >= (p_mm - p_bm) * (x - x_bm);
            if (verbose) {
                std::cout << "ineq_0_right = " << ineq_0_right << endl;
            }
            if (ineq_0_right) {
                bool ineq_0_bottom = (x_bm - x_bl) * (p - p_bl) >= (p_bm - p_bl) * (x - x_bl);
                Panel* child_0 = &old_panels[child_inds_start];
                int child_0_bottom_nbr_ind = child_0->bottom_nbr_ind;
                if (ineq_0_bottom || child_0_bottom_nbr_ind < 0) {
                    bool ineq_0_left = (x_ml - x_bl) * (p - p_bl) <= (p_ml - p_bl) * (x - x_bl);
                    int child_0_left_nbr_ind = child_0->left_nbr_ind;
                    if (ineq_0_left || child_0_left_nbr_ind < 0) {
                        
                        subpanel_ind = child_inds_start;
                        if (verbose) {
                            cout << "in child 0, panel " << subpanel_ind << endl;
                        }
                    }
                    else {
                        subpanel_ind = child_0_left_nbr_ind;
                        if (panel->is_left_bdry && bcs==periodic_bcs) { x += Lx; }
                    }
                } else {
                    subpanel_ind = child_0_bottom_nbr_ind;
                }
            } else
            {
                bool ineq_2_bottom = (x_br - x_bm) * (p - p_bm) >= (p_br - p_bm) * (x - x_bm);
                Panel* child_2;
                if (panel->is_refined_p) {
                     child_2 = &old_panels[child_inds_start];
                } else { // panel is refined in x and v
                     child_2 = &old_panels[child_inds_start+2];
                }
                int child_2_bottom_nbr_ind = child_2->bottom_nbr_ind;
                if (ineq_2_bottom || child_2_bottom_nbr_ind < 0) {
                    bool ineq_2_right = (x_mr - x_br) * (p - p_br) <= (p_mr - p_br) * (x - x_br);
                    int child_2_right_nbr_ind = child_2->right_nbr_ind;
                    if (! ineq_2_right || child_2_right_nbr_ind < 0) {
                        if (panel->is_refined_p) {
                            subpanel_ind = child_inds_start;
                        } else {
                            subpanel_ind = child_inds_start + 2;
                        }
                        if (verbose) {
                            cout << "in child 2, panel " << subpanel_ind << endl;
                        }
                    }
                    else {
                        subpanel_ind = child_2_right_nbr_ind;
                        if (panel->is_right_bdry && bcs == periodic_bcs) {
                            x -= Lx;
                        }
                    }
                } else {
                    subpanel_ind = child_2_bottom_nbr_ind;
                }
            }
            
        }

        leaf_ind = find_leaf_containing_xp_recursively(x, p, beyond_boundary, subpanel_ind, verbose);

    }
    return leaf_ind;
}


void AMRStructure::shift_xs(std::vector<double>& shifted_xs, const std::vector<double>& xs, const std::vector<double>& ps) {
    bool verbose = false;

    // if (iter_num >58) { verbose = true; }
    // shifted_xs = std::vector<double> (xs.size() );
    double x_bl, x_tl, x_br, x_tr;
    double p_bl, p_tl, p_br, p_tr;
    x_bl = this->old_xs[0]; p_bl = this->old_ps[0];
    x_tl = this->old_xs[2]; p_tl = this->old_ps[2];
    x_br = this->old_xs[6]; p_br = this->old_ps[6];
    x_tr = this->old_xs[8]; p_tr = this->old_ps[8];
    if (verbose) {
        cout << "(x,p)_bl (" << x_bl << ", " << p_bl << ")" << endl;
        cout << "(x,p)_tl (" << x_tl << ", " << p_tl << ")" << endl;
        cout << "(x,p)_br (" << x_br << ", " << p_br << ")" << endl;
        cout << "(x,p)_tr (" << x_tr << ", " << p_tr << ")" << endl;
    }

    for (int ii = 0; ii < xs.size(); ++ii) {
        double p=ps[ii];
        double x_temp = xs[ii];
        // trouble shooting in amr
        // if ( (fabs(x_temp) < 0.1 || fabs(x_temp -12.5664) < 0.1) && fabs(p -1.125) < 0.1) { verbose = true; }
        // else {verbose = false; }
        //end troubleshoot
        bool ineq_00_left = (x_tl - x_bl) * (p - p_bl) <= (p_tl - p_bl) * (x_temp - x_bl);


        if (verbose) {
            cout << "point " << ii << ": (x,p)= (" << x_temp << ", " << p << ")" << endl;
            cout << "ineq_00_left, " << ineq_00_left << endl; 
        }
        int counter = 0;
        while (not ineq_00_left) {
            x_temp += Lx;
            ineq_00_left = (x_tl - x_bl) * (p - p_bl) <= (p_tl - p_bl) * (x_temp - x_bl);
            if(verbose) {
                cout << "post shift x= (" << x_temp << ", ineq_00_left, " << ineq_00_left << endl; 
            }
            counter++;
            if (counter > 10) {
                throw std::runtime_error("too many shifts!");
            }
        }
        bool ineq_00_right = (x_tr - x_br) * (p - p_br) > (p_tr - p_br) * (x_temp - (x_br));
        if (verbose) {
            cout << "ineq_00_right, " << ineq_00_left << endl; 
        }
        counter = 0;
        while (not ineq_00_right) {
            x_temp -= Lx;
            ineq_00_right = (x_tr - x_br) * (p - p_br) > (p_tr - p_br) * (x_temp - (x_br));
            if(verbose) {
                cout << "post shift x= (" << x_temp << ", ineq_00_right, " << ineq_00_right << endl; 
            }
            counter++;
            if (counter > 10) {
                throw std::runtime_error("too many shifts!");
            }
        }
        shifted_xs[ii] = x_temp;
    }
}

void AMRStructure::interpolate_to_initial_xps(
    std::vector<double>& fs, std::vector<double>& xs, std::vector<double>& ps, 
    int nx, int np, bool verbose) 
{
    // bool print_profile = false;

    auto start = high_resolution_clock::now();
    #ifdef DEBUG
    // verbose = true;
    #endif
    std::vector<double> shifted_xs(xs.size() );
    if (bcs == periodic_bcs) {
#ifdef DEBUG
cout << "shifting " << endl;
#endif
        shift_xs(shifted_xs, xs, ps);
        if (verbose) {
            cout << "shift xs" << endl;
            std::copy(shifted_xs.begin(), shifted_xs.end(), std::ostream_iterator<double>(cout, ", "));
            cout << endl;
        }
    }
    else { // open bcs
        shifted_xs = xs;
    }
    auto stop = high_resolution_clock::now();
    // auto duration = duration_cast<microseconds>(stop - start);
    // if (print_profile) {
    //     cout << "shift time " << duration.count() << " microseconds" << endl;
    // }

    // have to sort points:
    start = high_resolution_clock::now();

#ifdef DEBUG
cout << "sorting " << endl;
#endif
    std::vector<int> sort_indices(xs.size());
    for (int ii = 0; ii < xs.size(); ii++ ) { sort_indices[ii] = ii; }
    // std::iota(sort_indices.begin(), sort_indices.end(), 0);
    bool do_sort = true;
    double sort_threshold = initial_dp / 10.0;
    if (do_sort) {
        std::sort(sort_indices.begin(), sort_indices.end(),
            [&] (int a, int b) 
            { 
                if (fabs(ps[a] - ps[b]) >= sort_threshold) { return ps[a] < ps[b]; }
                else {
                    return shifted_xs[a] < shifted_xs[b];
                }
            });
        // std::cout << "sorted indices " << endl;
        // std::copy(sort_indices.begin(), sort_indices.end(), std::ostream_iterator<int>(cout, " "));
        // std::cout << endl;
        // cout << "xs.size " << xs.size() << endl;
        // cout << "sort_indices.size " << sort_indices.size() << endl;
        // cout << "shifted_xs.size " << shifted_xs.size() << endl;
    }
#ifdef DEBUG
cout << "Done sorting" << endl;
#endif
    std::vector<double> sortxs(shifted_xs.size()), sortps(ps.size());
    for (int ii = 0; ii < xs.size(); ii++) {
        sortxs[ii] = shifted_xs[sort_indices[ii]];
        sortps[ii] = ps[sort_indices[ii]];
    }
    // if (verbose) {
    #ifdef DEBUG_L3
        cout << "sorted xs, size = " << sortxs.size() << endl;
        std::copy(sortxs.begin(), sortxs.end(), std::ostream_iterator<double>(cout, ", "));
        cout << endl << "sorted ps, size = "  << sortps.size() << endl;
        std::copy(sortps.begin(), sortps.end(), std::ostream_iterator<double>(cout, ", "));
        cout << endl;
    #endif /* DEBUG_L3 */
    // }

    stop = high_resolution_clock::now();
    // duration = duration_cast<microseconds>(stop - start);
    // if (print_profile) {
    //     cout << "sort time " << duration.count() << " microseconds" << endl;
    // }
    
    auto search_start = high_resolution_clock::now();
    std::vector<double> sortfs(xs.size() );

    std::vector<int> leaf_panel_of_points(xs.size() );
    std::vector<std::vector<int> > point_in_leaf_panels_by_inds(old_panels.size() );

#ifdef DEBUG
// counting how many xs to ps
// int ind_debug = 0;
// std::vector<int> size_p;
// std::vector<double> unique_ps;
// size_p.reserve(np);
// unique_ps.reserve(np);
// while (ind_debug < sortxs.size() ) {
//     int cntr = 0;
//     double debug_p = sortps[ind_debug];
//     unique_ps.push_back(debug_p);
//     cntr++;
//     int ind_debug2 = ind_debug + 1;
//     while (fabs(sortps[ind_debug2] - debug_p) < sort_threshold) {
//         ind_debug2++;
//         cntr++;
//     }
//     ind_debug += cntr;
//     size_p.push_back(cntr);
// }
// cout << "np = " << np << ", and there are " << size_p.size() << " distinct p values" << endl;
// for (int ii = 0; ii < size_p.size(); ++ii) {
//     cout << "at v= " << unique_ps[ii] << " there are " << size_p[ii] << " points" << endl;
// }
// verbose = false;
cout << "finding panel of first point" << endl;
// verbose=true;
#endif
    bool beyond_boundary = false;
    int leaf_ind = find_leaf_containing_xp_recursively(sortxs[0], sortps[0], beyond_boundary, 0, verbose);
#ifdef DEBUG
cout << "found first panel" << endl;
#endif
    std::vector<int> first_column_leaf_inds(np);
    // std::vector<int> first_column_ind_in_old_mesh(np);

#ifdef DEBUG 
cout << "searching first column" << endl;
// if (iter_num >= 4) { verbose = true;} 
    // if (verbose) {
        std::cout << "nx x np= " << nx << " x " << np << endl;
        cout << "xs size: " << xs.size() << endl;
    // }
// cout << "First column points" << endl;
// for (int ii = 0; ii < np; ++ii) {
//     int point_ind = ii * np;
//     cout << "sort point " << point_ind << " (x,p)=(" << sortxs[point_ind] << ", " << sortps[point_ind] << ")" << endl;
// }
#endif

    // if (bcs == periodic_bcs) {
        for (int ii =0; ii < np; ++ii) {
            beyond_boundary = false;
            int point_ind = ii * nx;
            std::set<int> history;
            history.emplace(leaf_ind);
            #ifdef DEBUG_L2
            cout << "testing point " << point_ind << ", x= " << sortxs[point_ind] << ", p= " << sortps[point_ind] << endl;
            #endif
            leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortps[point_ind], beyond_boundary, 
                                                                leaf_ind, history, verbose);
            first_column_leaf_inds[ii] = leaf_ind;
            // point_in_leaf_panels_by_inds[leaf_ind].push_back(point_ind);
            if (beyond_boundary) {
                leaf_panel_of_points[point_ind] = 0;
            } else {
                leaf_panel_of_points[point_ind] = leaf_ind;
            }
            // cout << "point " << sort_indices[point_ind] << " in panel " << leaf_ind << " (sorted ind " << point_ind << ")" << endl;
        }
//     } else { // open bcs
//         for (int ii =0; ii < np; ++ii) {
//             int jj = 0;
//             int point_ind = ii * np;
//             std::set<int> history;
//             history.emplace(leaf_ind);
//             int temp_leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortps[point_ind], leaf_ind, history, verbose);
//             leaf_panel_of_points[point_ind] = temp_leaf_ind;

//             while (temp_leaf_ind == 0 && jj < nx) {
// #ifdef DEBUG
// cout << "temp_leaf_ind was assigned 0!" << endl;
// #endif
//                 jj++;
//                 point_ind++;
//                 temp_leaf_ind = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortps[point_ind], leaf_ind, history, verbose);
//                 leaf_panel_of_points[point_ind] = temp_leaf_ind;
//             }
//             leaf_ind = temp_leaf_ind;
//             first_column_ind_in_old_mesh[ii] = jj;
//         }    
//     }
#ifdef DEBUG
verbose = false;
cout << "after first column" << endl;
// cout << "first column ind in old mesh" << endl;
// std::copy(first_column_ind_in_old_mesh.begin(), first_column_ind_in_old_mesh.end(),
//         std::ostream_iterator<int>(cout, " "));
// cout << endl;

//     cout << "leaf_panel_of_points " << endl;
//     for (int ii = 0; ii < xs.size(); ++ii) {
//         cout << "point (sorted ind) " << ii <<", unsorted ind " << sort_indices[ii] << ": (x,p)=(" << sortxs[ii] << ", " << sortps[ii] << ") is in panel " << leaf_panel_of_points[ii] << endl;
//     }
    // if (verbose) {
        cout << "panel search" << endl;
    // }
#endif
// #pragma omp parallel
// {
    // int num_threads = omp_get_num_threads();
    // if (omp_get_thread_num() == 0) {
    //     cout << "Number of threads from calc_E: " << num_threads << endl;
    // }
    // #pragma omp for
    for (int ii = 0; ii < np; ++ii) {
        // int jj0 = first_column_ind_in_old_mesh[ii];
        int jj0 = 0;
        int point_ind = ii * nx + jj0;
        int leaf_ind_c = first_column_leaf_inds[ii];//leaf_panel_of_points[point_ind];
        for (int jj = jj0+1; jj < nx; ++jj) {
            beyond_boundary = false;
            point_ind++;
            #ifdef DEBUG_L2
            // if (verbose) {
            cout << "Testing point " << point_ind << ", (x,p)= (" << sortxs[point_ind] << ", " << sortps[point_ind] << ")" << endl;
            // }
            #endif
            std::set<int> history;
            history.emplace(leaf_ind_c);
            leaf_ind_c = find_leaf_containing_point_from_neighbor(sortxs[point_ind], sortps[point_ind], beyond_boundary, leaf_ind_c, history, verbose);
            // point_in_leaf_panels_by_inds[leaf_ind_c].push_back(point_ind);
            if (beyond_boundary) {
                leaf_panel_of_points[point_ind] = 0;
            } else {
                leaf_panel_of_points[point_ind] = leaf_ind_c;
            }
            #ifdef DEBUG_L3
            cout << "sorted ind " << point_ind << ", ind " << sort_indices[point_ind] << " is in panel " << leaf_ind_c << endl;
            #endif
        }
    }    
// } // end omp parallel
#ifdef DEBUG
cout << "done with panel search" << endl;
#endif
#ifdef DEBUG_L3
    cout << "leaf_panel_of_points " << endl;
    for (int ii = 0; ii < xs.size(); ++ii) {
        cout << "point (sorted ind) " << ii <<", unsorted ind " << sort_indices[ii] << ": (x,p)=(" << sortxs[ii] << ", " << sortps[ii] << ") is in panel " << leaf_panel_of_points[ii] << endl;
    }
#endif

    for (int ii = 0; ii < leaf_panel_of_points.size(); ++ii) {
        point_in_leaf_panels_by_inds[leaf_panel_of_points[ii]].push_back(ii);
    }

#ifdef DEBUG_L3

    cout << "point in leaf panels by inds " << endl;
    for (auto ii = 0; ii < point_in_leaf_panels_by_inds.size(); ++ii) {
        std::vector<int> p_leaf = point_in_leaf_panels_by_inds[ii];
        if (p_leaf.size() > 0 ) {
            cout << "points in panel " << ii << ": ";
            for (auto jj = 0; jj < p_leaf.size(); ++jj) {
                cout << p_leaf[jj] << ", ";
            }
            cout << endl;
        }
    }
    // std::copy(leaf_panel_of_points.begin(), leaf_panel_of_points.end(), std::ostream_iterator<int>(cout, " "));
// debug
    // cout << "sort_indices[498] = " << sort_indices[498] << endl;
    // cout << "leaf_panel_of_points[498] = " << leaf_panel_of_points[498] << endl;
    // cout << "leaf_panel_of_points[sort_indicies[498]] = " << leaf_panel_of_points[sort_indices[498]] << endl;
//
#endif
    auto search_stop = high_resolution_clock::now();
    add_time(search_time, duration_cast<duration<double>>(search_stop - search_start) );
    
#ifdef DEBUG
cout << "eval interpolant" << endl;
#endif
    start = high_resolution_clock::now();
    for (int panel_ind = 0; panel_ind < old_panels.size(); panel_ind++) {
        if (point_in_leaf_panels_by_inds[panel_ind].size() > 0) {
            interpolate_from_panel_to_points(sortfs,sortxs,sortps,point_in_leaf_panels_by_inds[panel_ind], panel_ind, use_limiter, limit_val);
        }
    }
    for (int ii = 0; ii < fs.size(); ii++) {
        // debugging
        // if (sort_indices[ii] == 498) {
            // cout << "(x,p)_498 = (" << xs[498] << ", " << ps[498] << ")" << endl;
            // cout << "fs(498) = " << fs[498] << endl;
        // }
        // end debugging
        fs[sort_indices[ii]] = sortfs[ii];
    }
    // debugging
    // cout << "done unpacking sort" << endl;
    // cout << "something's up with sort_indices: " << endl;
    // std::copy(sort_indices.begin(), sort_indices.end(), std::ostream_iterator<int>(cout, " "));
    // cout << endl;
    // cout << "(x,p)_498 = (" << xs[498] << ", " << ps[498] << ")" << endl;
    // cout << "fs(498) = " << fs[498] << endl;
    // end debug element
    stop = high_resolution_clock::now();
    add_time(eval_time, duration_cast<duration<double>>(stop - start) );

    #ifdef DEBUG
    cout << "Done evaluating interpolant and done interpolating onto grid" << endl;
    #endif
}

void AMRStructure::interpolate_from_panel_to_points(
    std::vector<double>& values, std::vector<double>& xs, std::vector<double>& ps,
    std::vector<int>& point_inds, int panel_ind, bool use_limiter, double limit_val) 
{
    if (panel_ind == 0) { // if we are extrapolating beyond boundaries; assume 0
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            values[point_inds[ii]] = 0;
        }
    }
    else {
        bool verbose = false;
        #ifdef DEBUG
        cout << "in interpolate from panel to points" << endl;
        // verbose=true;
        #endif
        // troubleshooting interpolation trouble with amr
        // if (panel_ind == 2459) { verbose = true; }
        // end troubleshooting verbosity change
        Panel* panel = &(old_panels[panel_ind]);
        const int* panel_point_inds = panel->point_inds;
        double panel_xs[9], panel_ps[9], panel_fs[9];

        // if (verbose) {
        #ifdef DEBUG_L2
            cout << "From panel " << panel_ind << endl;
            cout << "Getting panel points" << endl;
        #endif
        // }
        for (int ii = 0; ii < 9; ++ii) {
            // Particle* part = &particles[vertex_inds[ii]];
            int pind = panel_point_inds[ii];
            // panel_xs[ii] = part->get_x();
            // panel_ps[ii] = part->get_p();
            // panel_fs[ii] = part->get_f();
            panel_xs[ii] = old_xs[pind];
            panel_ps[ii] = old_ps[pind];
            panel_fs[ii] = old_fs[pind];
        }

        if (do_unshear) {
            for (int ii = 0; ii < 9; ++ii) {
                #ifdef DEBUG
                cout << "panel_xs[ii] " << panel_xs[ii] << "  ";
                #endif
                panel_xs[ii] -= dt * (panel_ps[ii] - panel_ps[4]); // assumes remesh frequency = 1
#ifdef DEBUG

cout << "unshear panel x " << panel_xs[ii] << ", dt " << dt << ", p " << panel_ps[ii] << ", pmid " << panel_ps[4] << endl;
#endif
            }
            for (int ii = 0; ii < xs.size(); ++ii) {
                #ifdef DEBUG
                cout << "xs[ii]" << xs[ii] << "  ";
                #endif
                xs[ii] -= dt * (ps[ii] - panel_ps[4]);
#ifdef DEBUG

cout << "unshear xs[" << ii << "] " << xs[ii] << ", dt " << dt << ", ps[ii] " << ps[ii] << ", pmid " << panel_ps[4] << endl;
#endif
            }
            //std::transform(xs.begin(), xs.end(), ps.begin(), xs.begin(), [&](double a, double b) {return a - (b - panel_ps[4])*dt; });
        }
#ifdef DEBUG
    // if (iter_num >= 240) {
        // for (int ii = 0; ii < 9; ++ii) {
        //     if (panel_ps[ii] >= 0.013) {
        //         if (panel_xs[ii] >= 0.005 && panel_xs[ii] <= 0.015) {
        //             // cout << "(x,p,f)_" << ii << "=(" << xs[ii] << ", " << ps[ii] << ", " << fs[ii] <<")"<<endl;
        //             if (panel_fs[ii] > 15 || panel_fs[ii] < -0.5) {
        //                 verbose = true;
        //             }
        //         }
        //     }
        // }
        // for (int ii =0; ii < xs.size(); ++ii) {
        //     if (point_inds[ii] == 9172 || point_inds[ii] == 9272 || point_inds[ii] == 9275  || point_inds[ii] == 9234) {verbose = true; }
        // }
        
        // for (int ii =0; ii < xs.size(); ++ii) {
        //     if (point_inds[ii] == 8896) {
        //             cout << "(x,p,f)_" << ii << "=(" << xs[ii] << ", " << ps[ii] << ", " << fs[ii] <<") is in panel " << panel_ind<<endl;
        //         }
        // }
    // }

#endif

        // if (verbose) {
#ifdef DEBUG_L3
            std::cout << "Interpolating from " << std::endl;
            for (int ii = 0; ii < 9; ++ii) {
                std::cout << "(x,p,f)=(" << panel_xs[ii] << ", " << panel_ps[ii] << ", " << panel_fs[ii] << ")" << std::endl;
            }
        
            std::cout << "onto ";
            for (int ii = 0; ii < point_inds.size(); ii++) {
                int pind = point_inds[ii];
                std::cout << "point " << pind << ": (" << xs[pind] << ", " << ps[pind] << ")\n";
            }
#endif
        // }

        // double dx, dp, 
        double panel_dx[9], panel_dp[9];
        // dx = x - panel_xs[4];
        // dp = p - panel_ps[4];
        std::vector<double> dxs(point_inds.size()), dps(point_inds.size());
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            int pind = point_inds[ii];
            dxs[ii] = xs[pind] - panel_xs[4];
            dps[ii] = ps[pind] - panel_ps[4];
        }
        for (int ii = 0; ii < 9; ii ++) {
            panel_dx[ii] = panel_xs[ii] - panel_xs[4];
            panel_dp[ii] = panel_ps[ii] - panel_ps[4];
        }

        if (verbose) {
            std::cout << "test point distance from midpoint:" << std::endl;
            for (int ii = 0; ii < point_inds.size(); ++ii) {
                std::cout << ii <<": dx=" << dxs[ii] <<", dp=" << dps[ii] << std::endl;
            } 
            std::cout << "panel vertex distances from midpoint:" << std::endl;
            for (int ii = 0; ii < 9; ii++ ) {
                std::cout << ii << ": " << panel_dx[ii] << ", " << panel_dp[ii] << std::endl;
            }
        }

        Eigen::Matrix<double,9,9> A;
        for (int ii = 0; ii < 9; ++ii) {
            A(ii,0) = 1; A(ii,1) = panel_dx[ii];
            A(ii,2) = panel_dx[ii] * panel_dp[ii];
            A(ii,3) = panel_dp[ii];
            A(ii,4) = panel_dx[ii] * panel_dx[ii];
            A(ii,5) =  panel_dx[ii] * panel_dx[ii] * panel_dp[ii];
            A(ii,6) = panel_dx[ii] * panel_dx[ii] * panel_dp[ii] * panel_dp[ii];
            A(ii,7) = panel_dx[ii] * panel_dp[ii] * panel_dp[ii];
            A(ii,8) = panel_dp[ii] * panel_dp[ii];
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
            double dpsq = dps[ii] * dps[ii];
            Dx(ii,0) = 1; Dx(ii,1) = dxs[ii];
            Dx(ii,2) = dxs[ii] * dps[ii];
            Dx(ii,3) = dps[ii];
            Dx(ii,4) = dxsq;
            Dx(ii,5) =  dxsq * dps[ii];
            Dx(ii,6) = dxsq * dpsq;
            Dx(ii,7) = dxs[ii] * dpsq;
            Dx(ii,8) = dpsq;
        }
        Eigen::Matrix<double, Dynamic,1> interp_vals = Dx * c;
        if (verbose) {
            std::cout << "Here is the result:" << std::endl << interp_vals << std::endl;
        }
        for (int ii = 0; ii < point_inds.size(); ++ii) {
            values[point_inds[ii]] = interp_vals(ii);
        }

        if (use_limiter) {
            for (int ii = 0; ii < values.size(); ++ii) {
                if (values[ii] < 0) { values[ii] = limit_val; }
            }
        }
    }

    // return c(0) + c(1)*dx + c(2) * dx*dp + c(3) * dp +
    //         c(4) * dx*dx + c(5) * dx*dx*dp + c(6) * dx*dx*dp*dp +
    //         c(7) * dx*dp*dp + c(8) * dp*dp;

}

double AMRStructure::interpolate_from_panel(double x, double p, int panel_ind, bool verbose) {
    if (panel_ind == 0) { return 0.0; }
    else {
        Panel* panel = &(old_panels[panel_ind]);
        const int* point_inds = panel->point_inds;
        double panel_xs[9], panel_ps[9], panel_fs[9];

        for (int ii = 0; ii < 9; ++ii) {
            int pind = point_inds[ii];
            panel_xs[ii] = old_xs[pind];
            panel_ps[ii] = old_ps[pind];
            panel_fs[ii] = old_fs[pind];
        }
        if (do_unshear) {

            double gamma = sqrt(1 + p*p);
            double v = p * q / qm / gamma; 
            double panel_gammas[9], panel_vs[9];
            for (int ii = 0; ii < 9; ++ii) {
                panel_gammas[ii] = sqrt(1 + panel_ps[ii]*panel_ps[ii]);
                panel_vs[ii] = panel_ps[ii] * q / qm / panel_gammas[ii];
                panel_xs[ii] -= dt * (panel_vs[ii] - panel_vs[4]); // assumes remesh frequency = 1
#ifdef DEBUG

cout << "unshear panel x " << panel_xs[ii] << ", dt " << dt << ", p " << p << ", pmid " << panel_ps[4] << endl;
#endif
            }
            x -= dt * (v - panel_vs[4]);
#ifdef DEBUG

cout << "unshear x " << x << ", dt " << dt << ", p " << p << ", pmid " << panel_ps[4] << endl;
#endif

        }
        #ifdef DEBUG
        // if (iter_num >= 240) {
        //     if (panel_ind == 2964 || panel_ind == 2816 || panel_ind == 2892) {
        //         verbose = true;
        //     }
        // }
        // if (iter_num >= 240) {
        //     for (int ii = 0; ii < 9; ++ii) {
        //         if (panel_ps[ii] >= 0.013) {
        //             if (panel_xs[ii] >= 0.005 && panel_xs[ii] <= 0.015) {
        //                 // cout << "(x,p,f)_" << ii << "=(" << xs[ii] << ", " << ps[ii] << ", " << fs[ii] <<")"<<endl;
        //                 if (panel_fs[ii] > 15 || panel_fs[ii] < -0.5) {
        //                     verbose = true;
        //                 }
        //             }
        //         }
        //     }
            // for (int ii =0; ii < xs.size(); ++ii) {
            //     if (point_inds[ii] == 9172 || point_inds[ii] == 9272 || point_inds[ii] == 9275  || point_inds[ii] == 9234) {verbose = true; }
            // }
            
        //     for (int ii =0; ii < xs.size(); ++ii) {
        //         if (point_inds[ii] == 8896) {
        //                 cout << "(x,p,f)_" << ii << "=(" << xs[ii] << ", " << ps[ii] << ", " << fs[ii] <<") is in panel " << panel_ind<<endl;
        //             }
        //     }
        // }
        #endif /*DEBUG*/

        if (verbose) {
            std::cout << "Interpolating from " << std::endl;
            for (int ii = 0; ii < 9; ++ii) {
                std::cout << "(x,p,f)=(" << panel_xs[ii] << ", " << panel_ps[ii] << ", " << panel_fs[ii] << ")" << std::endl;
            }
        
            std::cout << "onto (" << x << ", " << p << ")\n";
        }

        double dx, dp, panel_dx[9], panel_dp[9];
        dx = x - panel_xs[4];
        dp = p - panel_ps[4];
        for (int ii = 0; ii < 9; ii ++) {
            panel_dx[ii] = panel_xs[ii] - panel_xs[4];
            panel_dp[ii] = panel_ps[ii] - panel_ps[4];
        }

        if (verbose) {
            std::cout << "test point distance from midpoint: dx=" << dx <<", dp=" << dp << std::endl;
            std::cout << "panel vertex distances from midpoint:" << std::endl;
            for (int ii = 0; ii < 9; ii++ ) {
                std::cout << ii << ": " << panel_dx[ii] << ", " << panel_dp[ii] << std::endl;
            }
        }

        Eigen::Matrix<double,9,9> A;
        for (int ii = 0; ii < 9; ++ii) {
            A(ii,0) = 1; A(ii,1) = panel_dx[ii];
            A(ii,2) = panel_dx[ii] * panel_dp[ii];
            A(ii,3) = panel_dp[ii];
            A(ii,4) = panel_dx[ii] * panel_dx[ii];
            A(ii,5) =  panel_dx[ii] * panel_dx[ii] * panel_dp[ii];
            A(ii,6) = panel_dx[ii] * panel_dx[ii] * panel_dp[ii] * panel_dp[ii];
            A(ii,7) = panel_dx[ii] * panel_dp[ii] * panel_dp[ii];
            A(ii,8) = panel_dp[ii] * panel_dp[ii];
        }
        Eigen::Map<Eigen::Matrix<double,9,1>> b(panel_fs);
        Eigen::Matrix<double,9,1> c = A.lu().solve(b);


        if (verbose) {
            std::cout << "Here is the matrix A:\n" << A << std::endl;
            std::cout << "Here is the f vector b:\n" << b << std::endl;
            std::cout << "Here is the coefficient vector c:\n" << c << std::endl;
        }

        double val = c(0) + c(1)*dx + c(2) * dx*dp + c(3) * dp +
                c(4) * dx*dx + c(5) * dx*dx*dp + c(6) * dx*dx*dp*dp +
                c(7) * dx*dp*dp + c(8) * dp*dp;
        if (verbose) { cout << "Result = " << val << endl; }
        return val;
    }
}

double AMRStructure::interpolate_from_mesh(double x, double p, bool verbose) {

#ifdef DEBUG
cout << "testing point (x,p)=(" << x <<", " << p << ")" << endl;
#endif

    // probably need to shift xs
    std::vector<double> xs(1,x);
    std::vector<double> shifted_xs(1,x);
    std::vector<double> ps(1,p);
    if (bcs==periodic_bcs) {
        shift_xs(shifted_xs, xs, ps);
    }
    double shifted_x = shifted_xs[0];

    bool beyond_boundary = false;
    int leaf_containing = find_leaf_containing_xp_recursively(shifted_x,p,beyond_boundary,0, verbose);
    if (beyond_boundary) {
        leaf_containing = 0;
    }
#ifdef DEBUG
cout << "in panel " << leaf_containing << endl;

if (iter_num >= 240) {
    if (fabs(x - 0.0118) < initial_dx/5 && fabs(p - 0.0143) < initial_dp/5) {
        cout << "(x,p)=(" << x << ", " << p << ") is in panel " << leaf_containing << endl;
        verbose = true;
    }
    if (fabs(x - 0.0054) < initial_dx/5 && fabs(p - 0.0145) < initial_dp/5) {
        cout << "(x,p)=(" << x << ", " << p << ") is in panel " << leaf_containing << endl;
        verbose = true;
    }
    if (fabs(x - 0.0093) < initial_dx/5 && fabs(p - 0.01337) < initial_dp/5) {
        cout << "(x,p)=(" << x << ", " << p << ") is in panel " << leaf_containing << endl;
        verbose = true;
    }
}
// find problem panel:
// if (iter_num >= 240) {
// }
#endif /* DEBUG */
    double val = interpolate_from_panel(shifted_x,p,leaf_containing, verbose);
    if (verbose) {
        cout << "(" << shifted_x << ", " << p << ") is in panel " << leaf_containing << ", f_interpolated(x,p) = " << val << endl;
    }
    #ifdef DEBUG
    cout << "f interpolated = " << val << endl;
    #endif
    return val;
}

void AMRStructure::interpolate_from_mesh(std::vector<double>& values, std::vector<double>& xs, std::vector<double>& ps, bool verbose) {
    std::vector<double> shifted_xs(xs.size());
    if (bcs == periodic_bcs) {
        shift_xs(shifted_xs, xs, ps);
    } else {
        shifted_xs = xs;
    }

    std::vector<int> leaves(xs.size());
    std::vector<std::vector<int> > point_in_leaf_panels_by_inds(old_panels.size() );
    for (int ii = 0; ii < xs.size(); ++ii) {
        bool beyond_boundary = false;
        int leaf_ind = find_leaf_containing_xp_recursively(shifted_xs[ii], ps[ii], beyond_boundary, 0, verbose);
        if (beyond_boundary) {
            leaf_ind = 0;
        } else {
            leaves[ii] = leaf_ind;
        }
        point_in_leaf_panels_by_inds[leaf_ind].push_back(ii);
    }

    for (int panel_ind = 0; panel_ind < old_panels.size(); panel_ind++) {
        if (point_in_leaf_panels_by_inds[panel_ind].size() > 0) {
            interpolate_from_panel_to_points(values,shifted_xs,ps,point_in_leaf_panels_by_inds[panel_ind], panel_ind, use_limiter, limit_val);
        }
    }

    // cout << endl << "interpolated fs: ";
    // std::copy(values.begin(), values.end(), std::ostream_iterator<double>(cout, " "));
}

void AMRStructure::interpolate_from_mesh_slow(std::vector<double>& values, std::vector<double>& xs, std::vector<double>& ps, bool verbose) {
    std::vector<double> shifted_xs;
    shift_xs(shifted_xs, xs, ps);

    std::vector<int> leaves;
    // std::vector<double> values;
    values = std::vector<double> ();
    for (int ii = 0; ii < xs.size(); ++ii) {
        bool beyond_boundary = false;
        leaves.push_back(find_leaf_containing_xp_recursively(shifted_xs[ii], ps[ii],beyond_boundary, 0, verbose));
        values.push_back(interpolate_from_panel(shifted_xs[ii], ps[ii], leaves[ii], verbose));
    }
    // cout << endl << "interpolated fs: ";
    // std::copy(values.begin(), values.end(), std::ostream_iterator<double>(cout, " "));
}

// double AMRStructure::InterpolateDistribution::operator() (double x, double p, bool verbose) {
//     // int leaf_containing = ::find_leaf_containing_xp_recursively(x,p,0, verbose);
//     // return ::interpolate_from_panel(x,p,leaf_containing, verbose);
//     this->find_leaf_containing_xp_recursively(x,p,0, verbose);
//     return 1.0;
// }
