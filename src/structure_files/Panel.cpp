/**
 * @file Panel.cpp
 * @author Ryan Sandberg (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-07-29
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#include "Panel.hpp"

// Constructors
Panel::Panel() 
    : panel_ind(0), level(0), parent_ind(-1), which_child(-1),
      left_nbr_ind(0), right_nbr_ind(0), top_nbr_ind(-2), bottom_nbr_ind(-2),
      is_left_bdry(false), is_right_bdry(false),
      needs_refinement(false), is_refined_xv(false), is_refined_v(false)
{
    for(int ii = 0; ii < 9; ii++) {
        point_inds[ii] = ii;
    }
}

Panel::Panel(int panel_ind, int level, int parent_ind, int which_child,
        int (&vpoint_inds)[9], 
        int ln_ind, int tn_ind, 
        int rn_ind, int bn_ind)
        : panel_ind(panel_ind), level(level), parent_ind(parent_ind), which_child(which_child),
        left_nbr_ind(ln_ind), top_nbr_ind(tn_ind), right_nbr_ind(rn_ind), bottom_nbr_ind(bn_ind),
        is_left_bdry(false), is_right_bdry(false),
        needs_refinement(false), is_refined_xv(false), is_refined_v(false)
{
    for(int ii = 0; ii < 9; ++ii) {
        this->point_inds[ii] = point_inds[ii];
    }
    child_inds_start = -1;
}

Panel::Panel(int panel_ind, int level, int parent_ind, int which_child,
            int p0, int p1, int p2, int p3, int p4, int p5, 
            int p6, int p7, int p8, 
            int ln_ind, int tn_ind, 
            int rn_ind, int bn_ind,
            bool is_left_bdry, bool is_right_bdry)
        : panel_ind(panel_ind), level(level), parent_ind(parent_ind), which_child(which_child),
        left_nbr_ind(ln_ind), top_nbr_ind(tn_ind), right_nbr_ind(rn_ind), bottom_nbr_ind(bn_ind),
        is_left_bdry(is_left_bdry), is_right_bdry(is_right_bdry),
        needs_refinement(false), is_refined_xv(false), is_refined_v(false)
{
    point_inds[0] = p0;
    point_inds[1] = p1;
    point_inds[2] = p2;
    point_inds[3] = p3;
    point_inds[4] = p4;
    point_inds[5] = p5;
    point_inds[6] = p6;
    point_inds[7] = p7;
    point_inds[8] = p8;

    child_inds_start = -1;
}
        
Panel::Panel(int panel_ind, int level, int parent_ind, int which_child,
            int ln_ind, int tn_ind, 
            int rn_ind, int bn_ind)
        : panel_ind(panel_ind), level(level), parent_ind(parent_ind), which_child(which_child),
        left_nbr_ind(ln_ind), top_nbr_ind(tn_ind), right_nbr_ind(rn_ind), bottom_nbr_ind(bn_ind),
        is_left_bdry(false), is_right_bdry(false),
        needs_refinement(false), is_refined_xv(false), is_refined_v(false)
{
    child_inds_start = -1;
}

Panel::Panel(int panel_ind, int level, int parent_ind, int which_child)
        : panel_ind(panel_ind), level(level), parent_ind(parent_ind), which_child(which_child),
        is_left_bdry(false), is_right_bdry(false),
        needs_refinement(false), is_refined_xv(false), is_refined_v(false)
{
    for(int ii = 0; ii < 8; ii++) {
        point_inds[ii] = -1;
    }

    child_inds_start = -1;
}
// End Constructors



void Panel::set_point_inds(int p0, int p1, int p2, int p3, int p4, 
                    int p5, int p6, int p7, int p8) {
    point_inds[0] = p0;
    point_inds[1] = p1;
    point_inds[2] = p2;
    point_inds[3] = p3;
    point_inds[4] = p4;
    point_inds[5] = p5;
    point_inds[6] = p6;
    point_inds[7] = p7;
    point_inds[8] = p8;
}

void Panel::set_child_inds_start(int c0) {
    is_refined_xv = true;
    needs_refinement = false;
    child_inds_start = c0;
}

void Panel::set_child_inds_start(int c0, bool refined_v) {
    is_refined_v = true;
    needs_refinement = false;
    child_inds_start = c0;
}



void Panel::print_panel() const {
    cout << "--------\n";
    cout << "Panel " << panel_ind << " is at level " << level << endl;
    cout << "It has point indices " << point_inds[0] << " " << point_inds[1] << " ";
    cout << point_inds[2] << " " << point_inds[3]  << " " << point_inds[4] << endl;
    cout << point_inds[5] << " " << point_inds[6] << " " << point_inds[7] << " " << point_inds[8] << endl;
    cout << "It is child " << which_child << " of panel " << parent_ind << endl;
    cout << "Neighbors are: left " << left_nbr_ind << ", top " << top_nbr_ind;
    cout << ", right " << right_nbr_ind << ", bottom " << bottom_nbr_ind << endl;
    if (needs_refinement) {
        cout << "Panel is flagged for refinement\n";
    }
    else {
        cout << "Panel is not flagged for refinemenment\n";
    }
    if (is_refined_v) {
        cout << "Panel is refined in v only\n";
        cout << "Children are panels " << child_inds_start << " " << child_inds_start + 1 << endl;
    } else {
        if (is_refined_xv) {
            cout << "Panel is refined\n";
            cout << "Children are panels " << child_inds_start << " " << child_inds_start + 1 << " ";
            cout << child_inds_start + 2 << " " << child_inds_start + 3 << endl;
        } else {
            cout << "Panel is not refined\n";
        }
    }
};


std::ostream& operator<<(std::ostream& os, const Panel& panel) {   cout << "--------\n";
    os << "Panel " << panel.panel_ind << " is at level " << panel.level << endl;
    os << "It has point indices ";
    for (int ii = 0; ii < 9; ii++) {
        os << panel.point_inds[ii] << ", ";
    }
    os << endl;
    os << "It is child " << panel.which_child << " of panel " << panel.parent_ind << endl;
    os << "Neighbors are: left " << panel.left_nbr_ind << ", top " << panel.top_nbr_ind;
    os << ", right " << panel.right_nbr_ind << ", bottom " << panel.bottom_nbr_ind << endl;
    if (panel.is_left_bdry) { os << "Panel is on left boundary" << endl; }
    if (panel.is_right_bdry) { os << "Panel is on right boundary" << endl; }
    if (panel.needs_refinement) {
        os << "Panel is flagged for refinement\n";
    }
    else {
        os << "Panel is not flagged for refinemenment\n";
    }
    if (panel.is_refined_v) {
        os << "Panel is refined in v only\n";
        os << "Children are " << panel.child_inds_start << " " << panel.child_inds_start + 1 << endl;
    } else if (panel.is_refined_xv) {
        os << "Panel is refined\n";
        os << "Children are " << panel.child_inds_start << " " << panel.child_inds_start + 1 << " ";
        os << panel.child_inds_start + 2 << " " << panel.child_inds_start + 3 << endl;
    } else {
        os << "Panel is not refined\n";
    }
    return os;
}
