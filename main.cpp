#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>

#include <igl/gaussian_curvature.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include "HalfedgeBuilder.cpp"
#include "SphereGeneration.cpp"
#include "LoopSubdivision.cpp"
#include "SystemSolver.cpp"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;

// This function is called every time a keyboard button is pressed
//bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
//{
//	if (key == '1')
//	{
//		HalfedgeBuilder *builder = new HalfedgeBuilder();
//		HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
//		SphereGeneration *generator = new SphereGeneration(V, F, he);   //
//        generator->subdivide();										 // perform one round subdivision
//        generator->print(0);

//		// update the current mesh
//        V = generator->getVertexCoordinates(); // update vertex coordinates
//        F = generator->getFaces();
//        viewer.data().clear();
//		viewer.data().set_mesh(V, F);
//		return true;
//	}
//	if (key == '2')
//	{
//		HalfedgeBuilder *builder = new HalfedgeBuilder();
//		HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
//		LoopSubdivision *loop = new LoopSubdivision(V, F, he);   //
//		loop->subdivide();										 // perform one round subdivision
//		loop->print(0);

//		// update the current mesh
//        V = loop->getVertexCoordinates(); // update vertex coordinates
//        F = loop->getFaces();
//        viewer.data().clear();
//		viewer.data().set_mesh(V, F);
//		return true;
//	}

//	if (key == 'S' || key == 's') // write the mesh to file (OFF format)
//	{
//		igl::writeOFF("../data/output.off", V, F);
//		return true;
//	}

//	return false;
//}

/**
 * Create a triangle mesh corresponding to an octagon inscribed in the unit circle
 */
void createOctagon(MatrixXd &Vertices, MatrixXi &Faces)
{
    int n = 4;
    int n_vertices = n*n;
    int n_faces = (n-1)*(n-1)*2;
    double side_size = 1.;

    Vertices = MatrixXd(n_vertices, 3);

    Faces = MatrixXi(n_faces, 3);

    for (int row = 0; row < n; row++) {
        for(int pt = row*n; pt < row*n + n; pt++) {
            Vertices.row(pt) = Vector3d((pt % n) * side_size, -row * side_size, 0.);
        }
    }

    int it_faces = 0;

    for(int row = 0; row < n-1; row++) {
        for(int pt = row*n; pt < row*n + n - 1; pt++) {
            Faces.row(it_faces) = Vector3i(pt, pt+1, pt+n);
            it_faces++;
            Faces.row(it_faces) = Vector3i(pt+1, pt+n, pt+n+1);
            it_faces++;
        }
    }

}


int vertexDegree(int v, HalfedgeDS& he)
{
    int degree = 0;

        for (int e = 0; e < he.sizeOfHalfedges(); e++) {
            if (he.getTarget(e) == v && he.getFace(e) != -1) {
                degree++;
            }
        }

    return degree;
}



// ------------ main program ----------------
int main(int argc, char *argv[])
{

	if (argc < 2)
	{
		std::cout << "Creating an octagon" << std::endl;
		createOctagon(V, F);
	}
	else
	{
		std::cout << "reading input file: " << argv[1] << std::endl;
		igl::readOFF(argv[1], V, F);
	}

	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	std::cout << "Press '1' for one round sphere generation" << std::endl
			  << "Press '2' for one round Loop subdivision" << std::endl
			  << "Press 'S' save the current mesh to file" << std::endl;


    //HalfedgeBuilder* builder = new HalfedgeBuilder();
    //HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation

//    for(int i = 0; i < he.sizeOfVertices(); i++) {
//        cout << "Vertex" << i << ": " << V.row(i) << ", degree: " << vertexDegree(i, he) << endl;
//    }

//    cout << "Nb of faces: " << he.sizeOfFaces() << endl;
//    cout << "Nb of halfedges: " << he.sizeOfHalfedges() << endl;
//    he.print();




    // CONSTANTS
    int n = 4;
    double triangle_size = 1.0;
    double density = 0.3;
    double gravity = 0.8;
    double time_step = 0.02;
    double k = 1.0;

    // INITIALIZE CLOTH
    Cloth* cl = new Cloth(n, triangle_size);

    // INITIALIZE HALFEDGE DATA STRUCTURE
    HalfedgeBuilder *builder = new HalfedgeBuilder();
    HalfedgeDS he = (builder->createMeshWithFaces(cl->getV().rows(), cl->getF()));
    cl->setHalfedgeDS(he);

    // INITIALIZE PARTICLES
    cl->initializeParticles(density);

    for(int i = 0; i < he.sizeOfVertices(); i++) {
        cout << "Vertex" << i << ": " << cl->getV().row(i) << ", degree: " << vertexDegree(i, he) << endl;
    }
    cout << "Nb of faces: " << cl->getHe()->sizeOfFaces() << endl;
    cout << "Nb of halfedges: " << cl->getHe()->sizeOfHalfedges() << endl;
    cl->getHe()->print();

    // INITIALIZE SYSTEM SOLVER
    SystemSolver solver(time_step, cl, gravity);
    solver.solve();

    // ANIMATION
    //viewer.callback_key_down = &key_down;
    viewer.data().set_mesh(cl->getV(), cl->getF());
    viewer.core().is_animating = true;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & )->bool // run animation
    {
      solver.solve();
      viewer.data().clear(); // Clear should be called before drawing the mesh
      viewer.data().set_mesh(cl->getV(), cl->getF());
      viewer.data().set_colors(Eigen::RowVector3d(0.3, 0.8, 0.3));
      return false;
    };

    viewer.launch(); // run the viewer

}
