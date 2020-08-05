

#include <AMRStructure.hpp>

void AMRStructure::copy_to_old() {
    old_panels = std::vector<Panel> ();
    // deep copy
    for (std::vector<Panel>::iterator it = panels.begin(); it != panels.end(); ++it) {
        old_panels.push_back(Panel(*it));
    }
    old_xs = xs;
    old_vs = vs;
    old_vs = fs;
}

void AMRStructure::remesh() {
    copy_to_old();
    generate_mesh([&] (double x, double v) { return interpolate_from_mesh(x,v,false);}, false);
}