#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "Particle.cpp"

#ifndef HALFEDGE_DS_HEADER
  #define HALFEDGE_DS_HEADER
  #include "HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;


class Cloth
{
public:

    Cloth(int cloth_dim, double triangle_side){

        // CREATE MESH

        n = cloth_dim;
        n_vertices = n*n;
        n_faces = (n-1)*(n-1)*2;
        side_size = triangle_side;
        a = side_size*side_size/2;

        V = MatrixXd(n_vertices, 3);

        F = MatrixXi(n_faces, 3);

        // Building V
        for (int row = 0; row < n; row++) {
            for(int pt = row*n; pt < row*n + n; pt++) {
                V.row(pt) = Vector3d((pt % n) * side_size, -row * side_size, 0.);
            }
        }

        // Building F
        int it_faces = 0;
        for (int row = 0; row < n-1; row++) {
            for (int pt = row*n; pt < row*n + n - 1; pt++) {
                F.row(it_faces) = Vector3i(pt, pt+1, pt+n);
                it_faces++;
                F.row(it_faces) = Vector3i(pt+1, pt+n, pt+n+1);
                it_faces++;
            }
        }

        // INITIALIZE HALFEDGE DATA STRUCTURE
        //HalfedgeBuilder *builder = new HalfedgeBuilder();
//        he = (builder->createMeshWithFaces(V.rows(), F));
//        for(int i = 0; i < he.sizeOfVertices(); i++) {
//            cout << "Vertex" << i << ": " << V.row(i) << ", degree: " << vertexDegree(i) << endl;
//        }
//        cout << "Nb of faces: " << he.sizeOfFaces() << endl;
//        cout << "Nb of halfedges: " << he.sizeOfHalfedges() << endl;
//        he.print();

    }

    void setHalfedgeDS(HalfedgeDS &mesh){
        he = &mesh;
    }

    void initializeParticles(double density){
        for (int i = 0; i < n_vertices; i++) {
            int degree = vertexDegree(i);
            double mass = (degree * density * a) * 1.0/3.0;
            Vector3d pos_3d = V.row(i);
            Vector2d pos_plane(pos_3d[0], pos_3d[1]);
            particles.push_back(Particle(mass, pos_plane, pos_3d, n_vertices));
        }
        particles[0].setNdof(0);
        particles[n-1].setNdof(0);
    }

    MatrixXd& getV(){
        return V;
    }

    MatrixXi& getF(){
        return F;
    }

    HalfedgeDS* getHe(){
        return he;
    }

    int getN(){
        return n;
    }

    vector<Particle>& getParticles(){
        return particles;
    }


private:

    int n, n_vertices, n_faces;
    double side_size;
    double a;
    MatrixXd V;
    MatrixXi F;
    HalfedgeDS* he;
    vector<Particle> particles;


    int vertexDegree(int v)
    {
        int degree = 0;

            for (int e = 0; e < he->sizeOfHalfedges(); e++) {
                if (he->getTarget(e) == v && he->getFace(e) != -1) {
                    degree++;
                }
            }

        return degree;
    }
};
