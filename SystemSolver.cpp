#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "Cloth.cpp"

using namespace Eigen;
using namespace std;


class SystemSolver
{
public:

    SystemSolver(double time_step, Cloth* cloth, double gravity){
        h = time_step;
        cl = cloth;
        n = cl->getV().rows();
        particles = cl->getParticles();
        g = gravity;
    }

    void solve(){
        cout << "entrei no solver" << endl;
        addWeights();

        delta_p = MatrixXd::Zero(3*n, 1);
        delta_v = MatrixXd::Zero(3*n, 1);

        calculate_M();
        calculate_p_zero();
        calculate_v_zero();
        calculate_f_zero();

        calculate_del_f_del_p();
        calculate_del_f_del_v();

        calculate_y_correction();
        calculate_z_correction();

        modified_pcg();
        updatePositions();
        updateVelocities();

    }

private:

    void addWeights(){
        for (int i = 0; i < n; i++) {
            double m = particles[i].mass();
            particles[i].addForce(Vector3d(0.0, -m*g, 0.0));
        }
    }

    void calculate_M(){
        M = MatrixXd::Zero(3*n, 3*n);
        for (int i = 0; i < n; i++) {
            M(3*i, 3*i) = particles[i].mass();
            M(3*i + 1, 3*i + 1) = particles[i].mass();
            M(3*i + 2, 3*i + 2) = particles[i].mass();
        }
    }

    void calculate_p_zero(){
        p_zero = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d pi = particles[i].getPosition();
            p_zero.block<3,1>(3*i,0) = pi;
        }
    }

    void calculate_v_zero(){
        v_zero = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d vi = particles[i].getVelocity();
            v_zero.block<3,1>(3*i,0) = vi;
        }
    }

    void calculate_f_zero(){
        f_zero = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d fi = particles[i].getForce();
            f_zero.block<3,1>(3*i,0) = fi;
        }
    }

    void calculate_del_f_del_p(){
        del_f_del_p = MatrixXd::Zero(3*n, 3*n);
        for (int j = 0; j < n; j++) {
            vector<MatrixXd> del_f_del_pj = particles[j].getForceDerivativePosition();
            for (int i = 0; i < n; i++) {
                del_f_del_p.block<3,3>(3*i, 3*j) = del_f_del_pj[i];
            }
        }

    }

    void calculate_del_f_del_v(){
        del_f_del_v = MatrixXd::Zero(3*n, 3*n);
        for (int j = 0; j < n; j++) {
            vector<MatrixXd> del_f_del_vj = particles[j].getForceDerivativeVelocity();
            for (int i = 0; i < n; i++) {
                del_f_del_v.block<3,3>(3*i, 3*j) = del_f_del_vj[i];
            }
        }
    }

    void calculate_y_correction(){
        y_correction = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d yi = particles[i].y_correc();
            y_correction.block<3,1>(3*i,0) = yi;
        }
    }

    void calculate_z_correction(){
        z_correction = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d zi = particles[i].z_correc();
            z_correction.block<3,1>(3*i,0) = zi;
        }
    }

    void modified_pcg(){
        // Initialize variables
        MatrixXd A, b, P, P_inverse, r, c, q, s;
        double delta_zero, delta_old, delta_new, alfa, eps = 0.1;

        A = M - h * del_f_del_v - h * h * del_f_del_p;
        b = h*(f_zero + h * del_f_del_p * v_zero + del_f_del_p * y_correction);
        calculate_P_and_P_inverse(P, P_inverse, A);

        // Procedure
        delta_v = z_correction;
        delta_zero = (filter(b).transpose() * P * filter(b)).determinant();
        r = filter(b - A * delta_v);
        c = filter(P_inverse * r);
        delta_new = (r.transpose() * c).determinant();
        while (delta_new > eps * eps * delta_zero) {
            q = filter(A * c);
            alfa = delta_new/((c.transpose() * q).determinant());
            delta_v += alfa * c;
            r -= alfa * q;
            s = P_inverse * r;
            delta_old = delta_new;
            delta_new = (r.transpose() * s).determinant();
            c = filter(s + (delta_new/delta_old) * c);
        }
    }

    void calculate_P_and_P_inverse(MatrixXd& P, MatrixXd& P_inverse, MatrixXd& A){
        P = MatrixXd::Zero(3*n, 3*n);
        P_inverse = MatrixXd::Zero(3*n, 3*n);
        for (int i = 0; i < 3*n; i++) {
            P(i,i) = 1./A(i,i);
            P_inverse(i,i) = A(i,i);
        }
    }

    MatrixXd filter(MatrixXd a){
        MatrixXd Si;
        for (int i = 0; i < n; i++) {
            if (particles[i].ndof() == 3) {
                Si = MatrixXd::Identity(3,3);
            } else if (particles[i].ndof() == 0) {
                Si = MatrixXd::Zero(3,3);
            }
            a.block<3,1>(3*i,0) = Si * a.block<3,1>(3*i,0);
        }
        return a;
    }

    void updateVelocities(){
        MatrixXd v = v_zero + delta_v;
        for (int i = 0; i < n; i++) {
            particles[i].setVelocity(v.block<3,1>(3*i,0));
        }
        cout << "v_zero: " << v_zero.block<3,1>(0,0) << endl;
        cout << "delta_v: " << delta_v.block<3,1>(0,0) << endl;
        cout << "v: " << v.block<3,1>(0,0) << endl;
        cout << "v(0): " << particles[0].getVelocity() << endl;
    }

    void updatePositions(){
        delta_p = h * (v_zero + delta_v) + y_correction;
        MatrixXd p = p_zero + delta_p;
        for (int i = 0; i < n; i++) {
            particles[i].setPosition3d(p.block<3,1>(3*i,0));
            cl->getV().row(i) = p.block<3,1>(3*i,0);
        }
        cout << "p(0): " << particles[0].getPosition() << endl;
        cout << endl;
    }


    double h;
    Cloth* cl;
    int n;
    double g;
    vector<Particle> particles;

    MatrixXd M;
    MatrixXd p_zero;
    MatrixXd v_zero;
    MatrixXd f_zero;

    MatrixXd del_f_del_p;
    MatrixXd del_f_del_v;

    MatrixXd y_correction;
    MatrixXd z_correction;

    MatrixXd delta_p;
    MatrixXd delta_v;
};
