#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

using namespace Eigen;
using namespace std;


class Particle
{
public:

    Particle(double mass, Vector2d pos_plane, int nb_of_particles){
        m = mass;
        liberty_degrees = 3;
        n = nb_of_particles;
        position_plane = pos_plane;
        velocity << 0.0, 0.0, 0.0;
        y_correction << 0.0, 0.0, 0.0;
        z_correction << 0.0, 0.0, 0.0;

        f << 0.0, 0.0, 0.0;

        for (int i = 0; i < n; i++) {
            MatrixXd zeros = MatrixXd::Zero(3, 3);
            del_f_del_p.push_back(zeros);
            del_f_del_v.push_back(zeros);
        }
    }

    double mass(){
        return m;
    }

    double u()
    {
        return position_plane[0];
    }

    double v()
    {
        return position_plane[1];
    }

    Vector3d& getVelocity(){
        return velocity;
    }

    Vector3d& y_correc(){
        return y_correction;
    }

    Vector3d& z_correc(){
        return z_correction;
    }

    Vector3d& getForce(){
        return f;
    }

    vector<MatrixXd>& getForceDerivativePosition(){
        return del_f_del_p;
    }

    vector<MatrixXd>& getForceDerivativeVelocity(){
        return del_f_del_v;
    }

    int ndof(){
        return liberty_degrees;
    }

    void setNdof(int ndof){
        liberty_degrees = ndof;
    }

    void setVelocity(Vector3d new_velocity) {
        velocity = new_velocity;
    }

    void addForce(Vector3d new_f) {
        f += new_f;
    }

    void addForceDerivativePosition(int i, Matrix3d new_deriv_p) {
        del_f_del_p[i] += new_deriv_p;
    }

    void addForceDerivativeVelocity(int i, Matrix3d new_deriv_v) {
        del_f_del_v[i] += new_deriv_v;
    }


private:

    double m;
    int liberty_degrees;
    int n;

    Vector2d position_plane;
    Vector3d velocity;
    Vector3d y_correction;
    Vector3d z_correction;

    Vector3d f;
    vector<MatrixXd> del_f_del_p;
    vector<MatrixXd> del_f_del_v;

};
