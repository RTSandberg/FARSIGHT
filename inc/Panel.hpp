// Panel.h
//
//
//
// Ryan Sandberg
//

#ifndef PANEL_H
#define PANEL_H

#include <iostream> // cout, endl
using std::cout;
using std::endl;

/**
 * @brief collection of information for a panel in an adaptively refined mesh
 * Attributes
 * ----------
 * panel_ind :
 * level
 * parent_ind
 * which_child
 * vertex_inds
 * left_nbr_ind
 * top_nbr_ind
 * right_nbr_ind
 * bottom_nbr_ind
 * 
 * needs_refinement
 * is_refined_xv
 * is_refined_v
 * child_inds_start
 * 
 * Notes
 * -----
 * Root panel `parent_ind` and `which_child` are -1
 * if panel is not refined, `child_ind_start` is -1
 * Neighbor indices are -1 if not determined yet
 * Neighbor indices are -2 if panel is a boundary and can't have neighbors
 */
struct Panel {
    int panel_ind;
    int point_inds[9];
    /**
     * @brief ordering of points
     * 
     */
    int level;
    int parent_ind;
    int which_child;
    int left_nbr_ind, top_nbr_ind, right_nbr_ind, bottom_nbr_ind;
    bool is_left_bdry, is_right_bdry;
    bool needs_refinement;
    bool is_refined_xv;
    bool is_refined_v;
    int child_inds_start;
    /**
     * @brief ordering of points and child indices
     * 
     * Stored in the format
     *  2 ----- 5 ----- 8
     *  |  [1]  |  [3]  |
     *  1 ----- 4 ----- 7
     *  |  [0]  |  [2]  |
     *  0 ----- 3 ----- 6
     * 
     *  v refined panel
     * Stored in the format
     *  2 ----- 5 ----- 8
     *  |      [1]      |
     *  1 ----- 4 ----- 7
     *  |      [0]      |
     *  0 ----- 3 ----- 6
     * 
     */


    // Constructor definitions
    Panel();
    Panel(int panel_ind, int level, int parent_ind, int which_child,
        int (&point_inds)[9], 
        int ln_ind, int tn_ind, 
        int rn_ind, int bn_ind);
    Panel(int panel_ind, int level, int parent_ind, int which_child,
        int p0, int p1, int p2, int p3, int p4, int p5, 
        int p6, int p7, int p8, 
        int ln_ind, int tn_ind, 
        int rn_ind, int bn_ind,
        bool is_left_bdry, bool is_right_bdry);
    Panel(int panel_ind, int level, int parent_ind, int which_child,
        int ln_ind, int tn_ind, 
        int rn_ind, int bn_ind);
    Panel(int panel_ind, int level, int parent_ind, int which_child);
    // End Constructors

    void set_point_inds(int p0, int p1, int p2, int p3, int p4, 
                        int p5, int p6, int p7, int p8);
    void set_child_inds_start(int c0);
    void set_child_inds_start(int c0, bool refined_v);

    void check_if_refinement_needed();

    void print_panel() const;
};

std::ostream& operator<<(std::ostream& os, const Panel& panel);

#endif /* PANEL_H */

