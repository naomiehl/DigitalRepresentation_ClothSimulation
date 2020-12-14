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

        // *** Original V (1st animation)
        for (int row = 0; row < n; row++) {
            for(int pt = row*n; pt < row*n + n; pt++) {
                V.row(pt) = Vector3d((pt % n) * side_size, -row * side_size, 0.);
            }
        }
        // *****************************
        // ** New V (2nd animation)
//        for (int row = 0; row < n; row++) {
//            for(int pt = row*n; pt < row*n + n; pt++) {
//                V.row(pt) = Vector3d((pt % n) * side_size * 0.8, -row * side_size * 0.8, 0.);
//            }
//        }
        // *****************************

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
    }

    void setHalfedgeDS(HalfedgeDS &mesh){
        he = &mesh;
    }

    void initializeParticles(double density){
        for (int row = 0; row < n; row++) {
            for(int pt = row*n; pt < row*n + n; pt++) {
                int degree = vertexDegree(pt);
                double mass = (degree * density * a) * 1.0/3.0;
                Vector2d pos_plane((pt % n) * side_size, -row * side_size);
                particles.push_back(Particle(mass, pos_plane, n_vertices));
            }
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

    double getA(){
        return a;
    }


private:

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


    int n, n_vertices, n_faces;
    double side_size;
    double a;

    MatrixXd V;
    MatrixXi F;
    HalfedgeDS* he;
    vector<Particle> particles;

};
