#include "AMRStructure.hpp"


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
    shifted_xs = std::vector<double> (xs.size() );
    double x_bl, x_tl, x_br, x_tr;
    double v_bl, v_tl, v_br, v_tr;
    x_bl = this->xs[0]; v_bl = this->vs[0];
    x_tl = this->xs[2]; v_tl = this->vs[2];
    x_br = this->xs[6]; v_br = this->vs[6];
    x_tr = this->xs[8]; v_tr = this->vs[8];

    for (int ii = 0; ii < xs.size(); ++ii) {
        double v=vs[ii];
        double x_temp = xs[ii];

        bool ineq_00_left = (x_tl - x_bl) * (v - v_bl) <= (v_tl - v_bl) * (x_temp - x_bl);
        while (not ineq_00_left) {
            x_temp += Lx;
            ineq_00_left = (x_tl - x_bl) * (v - v_bl) <= (v_tl - v_bl) * (x_temp - x_bl);
        }
        bool ineq_00_right = (x_tr - x_br) * (v - v_br) >= (v_tr - v_br) * (x_temp - (x_br));
        while (not ineq_00_right) {
            x_temp -= Lx;
            ineq_00_right = (x_tr - x_br) * (v - v_br) >= (v_tr - v_br) * (x_temp - (x_br));
        }
        shifted_xs.push_back(x_temp);
    }
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
    // std::vector<double> values;
    values = std::vector<double> ();
    for (int ii = 0; ii < xs.size(); ++ii) {
        leaves.push_back(find_leaf_containing_xv_recursively(shifted_xs[ii], vs[ii], 0, verbose));
        values.push_back(interpolate_from_panel(shifted_xs[ii], vs[ii], leaves[ii], verbose));
    }
}
// double AMRStructure::InterpolateDistribution::operator() (double x, double v, bool verbose) {
//     // int leaf_containing = ::find_leaf_containing_xv_recursively(x,v,0, verbose);
//     // return ::interpolate_from_panel(x,v,leaf_containing, verbose);
//     this->find_leaf_containing_xv_recursively(x,v,0, verbose);
//     return 1.0;
// }