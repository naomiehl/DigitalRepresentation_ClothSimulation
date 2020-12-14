#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include "HalfedgeBuilder.cpp"
#include "SystemSolver.cpp"

using namespace Eigen;
using namespace std;


/* ========================================================================= *

    Obs. 1:
        The function addDampingForcesStretch() in the SystemSolver solver()
        method was the one that was causing errors for us just before the
        submission, so at the time we wrote the report this function was not
        well implemented yet, we managed to make it work a little later.
        So if you want to run the exact simulation that we described in the
        report, just comment the line 28 of SystemSolver.cpp, where this
        function is called.

    Obs. 2:
        Also after the submission, we added a new simulation whose difference
        to the previous one is only the position where we initialise the
        vertices, this time we initilise them a little more far from what
        could be considered as the "rest state" of the cloth. We did this
        with the objectif of generating a simulation where the cloth would
        have a bigger movement tendency. To switch between one or another
        version, just uncomment the version desired in the Cloth class
        constructor. Both versions work better using the function
        addDampingForcesStretch(), as expected.

 * ========================================================================= */



// ------------ main program ----------------
int main(int argc, char *argv[])
{

    igl::opengl::glfw::Viewer viewer;


    // CONSTANTS
    int n = 8;
    double triangle_size = 0.5;
    double density = 0.3;
    double gravity = 9.8;
    double time_step = 0.02;
    double stiffness = 1000.0;

    // INITIALIZE CLOTH
    Cloth* cl = new Cloth(n, triangle_size);

    // INITIALIZE HALFEDGE DATA STRUCTURE
    HalfedgeBuilder *builder = new HalfedgeBuilder();
    HalfedgeDS he = (builder->createMeshWithFaces(cl->getV().rows(), cl->getF()));
    cl->setHalfedgeDS(he);

    // INITIALIZE PARTICLES
    cl->initializeParticles(density);

    // INITIALIZE SYSTEM SOLVER
    SystemSolver solver(time_step, cl, gravity, stiffness);

    // ANIMATION
    viewer.data().set_mesh(cl->getV(), cl->getF());
    viewer.core().is_animating = true;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
    {
      solver.solve();
      viewer.data().clear();
      viewer.data().set_mesh(cl->getV(), cl->getF());
      viewer.data().set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
      return false;
    };

    viewer.launch();

}
