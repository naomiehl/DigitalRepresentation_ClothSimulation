#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "Cloth.cpp"

using namespace Eigen;
using namespace std;


class SystemSolver
{
public:

    SystemSolver(double time_step, Cloth* cloth, double gravity, double stiffness){
        h = time_step;
        cl = cloth;
        n = cl->getV().rows();
        g = gravity;
        k_s = stiffness;
    }

    void solve(){

        vector<Particle>& particles = cl->getParticles();

        addWeights(particles);
        addStretchForces(particles);
        addShearingForces(particles);
        addDampingForcesStretch(particles);

        delta_p = MatrixXd::Zero(3*n, 1);
        delta_v = MatrixXd::Zero(3*n, 1);

        calculate_M(particles);
        calculate_p_zero();
        calculate_v_zero(particles);
        calculate_f_zero(particles);

        calculate_del_f_del_p(particles);
        calculate_del_f_del_v(particles);

        calculate_y_correction(particles);
        calculate_z_correction(particles);

        modified_pcg(particles);

        updatePositions();
        updateVelocities(particles);

    }

private:

    void addWeights(vector<Particle>& particles){
        for (int i = 0; i < n; i++) {
            double m = particles[i].mass();
            particles[i].addForce(Vector3d(0.0, -m*g, 0.0));
        }
    }

    void addStretchForces(vector<Particle>& particles){

        for (int face = 0; face < cl->getF().rows(); face++) {

            Vector3i triangle = cl->getF().row(face);
            int i = triangle[0], j = triangle[1], k = triangle[2];

            double det, facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v, facteur_k_u, facteur_k_v;
            Vector3d w_u, w_v, w_u_point, w_v_point;

            double bu = 1., bv = 1.;

            calculate_utilities_triangle(particles, i, j, k, w_u, w_v, w_u_point, w_v_point, det, facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v,
                                         facteur_k_u, facteur_k_v);

            // Calculate forces

            Vector3d fi, fj, fk;

            fi = calculate_fp_stretch(w_u, w_v, det, bu, bv, facteur_i_u, facteur_i_v);
            particles[i].addForce(fi);


            fj = calculate_fp_stretch(w_u, w_v, det, bu, bv, facteur_j_u, facteur_j_v);
            particles[j].addForce(fj);


            fk = calculate_fp_stretch(w_u, w_v, det, bu, bv, facteur_k_u, facteur_k_v);
            particles[k].addForce(fk);

            // Calculate force derivatives (position)

            MatrixXd del_fi_del_i, del_fi_del_j, del_fi_del_k;
            MatrixXd del_fj_del_i, del_fj_del_j, del_fj_del_k;
            MatrixXd del_fk_del_i, del_fk_del_j, del_fk_del_k;


            del_fi_del_i = calculate_delf_delp_stretch(facteur_i_u, facteur_i_v, facteur_i_u, facteur_i_v,
                                                       w_u, w_v, w_u_point, w_v_point, det, bu, bv);
            particles[i].addForceDerivativePosition(i, del_fi_del_i);


            del_fi_del_j = calculate_delf_delp_stretch(facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v,
                                                       w_u, w_v, w_u_point, w_v_point, det, bu, bv);
            particles[j].addForceDerivativePosition(i, del_fi_del_j);


            del_fi_del_k = calculate_delf_delp_stretch(facteur_i_u, facteur_i_v, facteur_k_u, facteur_k_v,
                                                       w_u, w_v, w_u_point, w_v_point, det, bu, bv);
            particles[k].addForceDerivativePosition(i, del_fi_del_k);


            del_fj_del_i = del_fi_del_j.transpose();
            particles[i].addForceDerivativePosition(j, del_fj_del_i);


            del_fj_del_j = calculate_delf_delp_stretch(facteur_j_u, facteur_j_v, facteur_j_u, facteur_j_v,
                                                       w_u, w_v, w_u_point, w_v_point, det, bu, bv);
            particles[j].addForceDerivativePosition(j, del_fj_del_j);


            del_fj_del_k = calculate_delf_delp_stretch(facteur_j_u, facteur_j_v, facteur_k_u, facteur_k_v,
                                                       w_u, w_v, w_u_point, w_v_point, det, bu, bv);
            particles[k].addForceDerivativePosition(j, del_fj_del_k);


            del_fk_del_i = del_fi_del_k.transpose();
            particles[i].addForceDerivativePosition(k, del_fk_del_i);


            del_fk_del_j = del_fj_del_k.transpose();
            particles[j].addForceDerivativePosition(k, del_fk_del_j);


            del_fk_del_k = calculate_delf_delp_stretch(facteur_k_u, facteur_k_v, facteur_k_u, facteur_k_v,
                                                       w_u, w_v, w_u_point, w_v_point, det, bu, bv);
            particles[k].addForceDerivativePosition(k, del_fk_del_k);
        }
    }

    Vector3d calculate_fp_stretch(Vector3d w_u, Vector3d w_v, double det, double bu, double bv,
                                  double facteur_p_u, double facteur_p_v){

        MatrixXd C = calculate_C_stretch(w_u, w_v, bu, bv);
        MatrixXd delC_delp = calculate_delC_delp_stretch(w_u, w_v, det, facteur_p_u, facteur_p_v);
        Vector3d fp = - k_s * delC_delp * C;

        return fp;
    }

    MatrixXd calculate_C_stretch(Vector3d w_u, Vector3d w_v, double bu, double bv){
        double norm_w_u = w_u.norm();
        double norm_w_v = w_v.norm();

        MatrixXd C = MatrixXd(2,1);
        C << norm_w_u - bu,
             norm_w_v - bv;
        C *= cl->getA();

        return C;
    }

    MatrixXd calculate_C_point_stretch(Vector3d w_u, Vector3d w_v,
                                       Vector3d w_u_point, Vector3d w_v_point){
        double norm_w_u = w_u.norm();
        double norm_w_v = w_v.norm();

        MatrixXd C_point = MatrixXd(2,1);
        C_point << 1./norm_w_u * w_u.dot(w_u_point),
                   1./norm_w_v * w_v.dot(w_v_point);
        C_point *= cl->getA();

        return C_point;
    }

    MatrixXd calculate_delC_delp_stretch(Vector3d w_u, Vector3d w_v, double det, double facteur_p_u, double facteur_p_v){
        double norm_w_u = w_u.norm();
        double norm_w_v = w_v.norm();

        double coeff_w_u = 1./(norm_w_u * det);
        double coeff_w_v = 1./(norm_w_v * det);

        double w_u_x = w_u[0], w_u_y = w_u[1], w_u_z = w_u[2];
        double w_v_x = w_v[0], w_v_y = w_v[1], w_v_z = w_v[2];

        MatrixXd delC_delp = MatrixXd(3,2);
        delC_delp << coeff_w_u * w_u_x * facteur_p_u, coeff_w_v * w_v_x * facteur_p_v,
                     coeff_w_u * w_u_y * facteur_p_u, coeff_w_v * w_v_y * facteur_p_v,
                     coeff_w_u * w_u_z * facteur_p_u, coeff_w_v * w_v_z * facteur_p_v;
        delC_delp *= cl->getA();

        return delC_delp;
    }

    MatrixXd calculate_delf_delp_stretch(double facteur_p1_u, double facteur_p1_v, double facteur_p2_u, double facteur_p2_v,
                                         Vector3d w_u, Vector3d w_v, Vector3d w_u_point, Vector3d w_v_point, double det, double bu, double bv){

        MatrixXd delC_delp1 = calculate_delC_delp_stretch(w_u, w_v, det, facteur_p1_u, facteur_p1_v);
        MatrixXd delC_delp2 = calculate_delC_delp_stretch(w_u, w_v, det, facteur_p2_u, facteur_p2_v);

        MatrixXd second_matrix = MatrixXd(3,3);

        second_matrix = second_matrix_delf_stretching(w_u, w_v, w_u_point, w_v_point, det, bu, bv,
                                                      facteur_p1_u, facteur_p1_v, facteur_p2_u, facteur_p2_v, 0);

        MatrixXd delf_delp = MatrixXd(3,3);

        delf_delp = delC_delp1 * delC_delp2.transpose() + second_matrix;

        delf_delp *= - k_s;

        return delf_delp;
    }

    Matrix3d second_matrix_delf_stretching(Vector3d w_u, Vector3d w_v, Vector3d w_u_point, Vector3d w_v_point,
                                           double det, double bu, double bv,
                                           double facteur_p1_u, double facteur_p1_v,
                                           double facteur_p2_u, double facteur_p2_v,
                                           int indicator){

        double w_u_x = w_u[0], w_u_y = w_u[1], w_u_z = w_u[2];
        double w_v_x = w_v[0], w_v_y = w_v[1], w_v_z = w_v[2];

        double norm_w_u = w_u.norm();
        double norm_w_v = w_v.norm();

        MatrixXd C_or_C_point;

        if (indicator == 0) {
            C_or_C_point = calculate_C_stretch(w_u, w_v, bu, bv);
        } else {
            C_or_C_point = calculate_C_point_stretch(w_u, w_v, w_u_point, w_v_point);
        }

        MatrixXd del2C_delp1_delp2x = MatrixXd(3,2);
        MatrixXd del2C_delp1_delp2y = MatrixXd(3,2);
        MatrixXd del2C_delp1_delp2z = MatrixXd(3,2);
        MatrixXd second_matrix = MatrixXd(3,3);

        del2C_delp1_delp2x(0,0) = ((- w_u_x * w_u_x)/pow(norm_w_u, 3) + 1./norm_w_u) * (facteur_p2_u * facteur_p1_u)/(det * det);
        del2C_delp1_delp2x(1,0) = -1./pow(norm_w_u, 3) * w_u_x * w_u_y * facteur_p2_u * facteur_p1_u / (det * det);
        del2C_delp1_delp2x(2,0) = -1./pow(norm_w_u, 3) * w_u_x * w_u_z * facteur_p2_u * facteur_p1_u / (det * det);
        del2C_delp1_delp2x(0,1) = ((- w_v_x * w_v_x)/pow(norm_w_v, 3) + 1./norm_w_v) * (facteur_p2_v * facteur_p1_v)/(det * det);
        del2C_delp1_delp2x(1,1) = -1./pow(norm_w_v, 3) * w_v_x * w_v_y * facteur_p2_v * facteur_p1_v / (det * det);
        del2C_delp1_delp2x(2,1) = -1./pow(norm_w_v, 3) * w_v_x * w_v_z * facteur_p2_v * facteur_p1_v / (det * det);
        del2C_delp1_delp2x *= cl->getA();

        del2C_delp1_delp2y(0,0) = -1./pow(norm_w_u, 3) * w_u_x * w_u_y * facteur_p2_u * facteur_p1_u / (det * det);
        del2C_delp1_delp2y(1,0) = ((- w_u_y * w_u_y)/pow(norm_w_u, 3) + 1./norm_w_u) * (facteur_p2_u * facteur_p1_u)/(det * det);
        del2C_delp1_delp2y(2,0) = -1./pow(norm_w_u, 3) * w_u_z * w_u_y * facteur_p2_u * facteur_p1_u / (det * det);
        del2C_delp1_delp2y(0,1) = -1./pow(norm_w_v, 3) * w_v_x * w_v_y * facteur_p2_v * facteur_p1_v / (det * det);
        del2C_delp1_delp2y(1,1) = ((- w_v_y * w_v_y)/pow(norm_w_v, 3) + 1./norm_w_v) * (facteur_p2_v * facteur_p1_v)/(det * det);
        del2C_delp1_delp2y(2,1) = -1./pow(norm_w_v, 3) * w_v_z * w_v_y * facteur_p2_v * facteur_p1_v / (det * det);
        del2C_delp1_delp2y *= cl->getA();

        del2C_delp1_delp2z(0,0) = -1./pow(norm_w_u, 3) * w_u_x * w_u_z * facteur_p2_u * facteur_p1_u / (det * det);
        del2C_delp1_delp2z(1,0) = -1./pow(norm_w_u, 3) * w_u_z * w_u_y * facteur_p2_u * facteur_p1_u / (det * det);
        del2C_delp1_delp2z(2,0) = ((- w_u_z * w_u_z)/pow(norm_w_u, 3) + 1./norm_w_u) * (facteur_p2_u * facteur_p1_u)/(det * det);
        del2C_delp1_delp2z(0,1) = -1./pow(norm_w_v, 3) * w_v_x * w_v_z * facteur_p2_v * facteur_p1_v / (det * det);
        del2C_delp1_delp2z(1,1) = -1./pow(norm_w_v, 3) * w_v_z * w_v_y * facteur_p2_v * facteur_p1_v / (det * det);
        del2C_delp1_delp2z(2,1) = ((- w_v_z * w_v_z)/pow(norm_w_v, 3) + 1./norm_w_v) * (facteur_p2_v * facteur_p1_v)/(det * det);
        del2C_delp1_delp2z *= cl->getA();

        second_matrix.col(0) = del2C_delp1_delp2x * C_or_C_point;
        second_matrix.col(1) = del2C_delp1_delp2y * C_or_C_point;
        second_matrix.col(2) = del2C_delp1_delp2z * C_or_C_point;

        return second_matrix;
    }

    void addShearingForces(vector<Particle>& particles){

        for (int face = 0; face < cl->getF().rows(); face++) {

            Vector3i triangle = cl->getF().row(face);
            int i = triangle[0], j = triangle[1], k = triangle[2];

            double det, facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v, facteur_k_u, facteur_k_v;
            Vector3d w_u, w_v, w_u_point, w_v_point;

            calculate_utilities_triangle(particles, i, j, k, w_u, w_v, w_u_point, w_v_point, det, facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v,
                                         facteur_k_u, facteur_k_v);

            // Calculate forces

            Vector3d fi, fj, fk;

            fi = calculate_fp_shearing(w_u, w_v, det, facteur_i_u, facteur_i_v);
            particles[i].addForce(fi);

            fj = calculate_fp_shearing(w_u, w_v, det, facteur_j_u, facteur_j_v);
            particles[j].addForce(fj);

            fk = calculate_fp_shearing(w_u, w_v, det, facteur_k_u, facteur_k_v);
            particles[k].addForce(fk);

            // Calculate force derivatives (position)

            MatrixXd del_fi_del_i, del_fi_del_j, del_fi_del_k;
            MatrixXd del_fj_del_i, del_fj_del_j, del_fj_del_k;
            MatrixXd del_fk_del_i, del_fk_del_j, del_fk_del_k;

            del_fi_del_i = calculate_delf_delp_shearing(facteur_i_u, facteur_i_v, facteur_i_u, facteur_i_v, w_u, w_v, det);
            particles[i].addForceDerivativePosition(i, del_fi_del_i);

            del_fi_del_j = calculate_delf_delp_shearing(facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v, w_u, w_v, det);
            particles[j].addForceDerivativePosition(i, del_fi_del_j);

            del_fi_del_k = calculate_delf_delp_shearing(facteur_i_u, facteur_i_v, facteur_k_u, facteur_k_v, w_u, w_v, det);
            particles[k].addForceDerivativePosition(i, del_fi_del_k);

            del_fj_del_i = del_fi_del_j.transpose();
            particles[i].addForceDerivativePosition(j, del_fj_del_i);

            del_fj_del_j = calculate_delf_delp_shearing(facteur_j_u, facteur_j_v, facteur_j_u, facteur_j_v, w_u, w_v, det);
            particles[j].addForceDerivativePosition(j, del_fj_del_j);

            del_fj_del_k = calculate_delf_delp_shearing(facteur_j_u, facteur_j_v, facteur_k_u, facteur_k_v, w_u, w_v, det);
            particles[k].addForceDerivativePosition(j, del_fj_del_k);

            del_fk_del_i = del_fi_del_k.transpose();
            particles[i].addForceDerivativePosition(k, del_fk_del_i);

            del_fk_del_j = del_fj_del_k.transpose();
            particles[j].addForceDerivativePosition(k, del_fk_del_j);

            del_fk_del_k = calculate_delf_delp_shearing(facteur_k_u, facteur_k_v, facteur_k_u, facteur_k_v, w_u, w_v, det);
            particles[k].addForceDerivativePosition(k, del_fk_del_k);
        }
    }

    Vector3d calculate_fp_shearing(Vector3d w_u, Vector3d w_v, double det, double facteur_p_u, double facteur_p_v){

        double C = calculate_C_shearing(w_u, w_v);
        MatrixXd delC_delp = calculate_delC_delp_shearing(w_u, w_v, det, facteur_p_u, facteur_p_v);
        Vector3d fp = - k_s * delC_delp * C;

        return fp;
    }

    double calculate_C_shearing(Vector3d w_u, Vector3d w_v){
        double C = cl->getA() * w_u.dot(w_v);
        return C;
    }

    Vector3d calculate_delC_delp_shearing(Vector3d w_u, Vector3d w_v, double det, double facteur_p_u, double facteur_p_v){
        Vector3d delC_delp;

        double w_u_x = w_u[0], w_u_y = w_u[1], w_u_z = w_u[2];
        double w_v_x = w_v[0], w_v_y = w_v[1], w_v_z = w_v[2];

        delC_delp << (facteur_p_u * w_v_x)/det + (w_u_x * facteur_p_v)/det,
                     (facteur_p_u * w_v_y)/det + (w_u_y * facteur_p_v)/det,
                     (facteur_p_u * w_v_z)/det + (w_u_z * facteur_p_v)/det;

        delC_delp *= cl->getA();

        return delC_delp;
    }

    MatrixXd calculate_delf_delp_shearing(double facteur_p1_u, double facteur_p1_v, double facteur_p2_u, double facteur_p2_v,
                                         Vector3d w_u, Vector3d w_v, double det){
        double C = calculate_C_shearing(w_u, w_v);

        Vector3d delC_delp1 = calculate_delC_delp_shearing(w_u, w_v, det, facteur_p1_u, facteur_p1_v);
        Vector3d delC_delp2 = calculate_delC_delp_shearing(w_u, w_v, det, facteur_p2_u, facteur_p2_v);

        MatrixXd del2C_delp1_delp2 = MatrixXd(3,3);

        del2C_delp1_delp2 << (facteur_p1_u * facteur_p2_v + facteur_p2_u * facteur_p1_v)/(det * det), 0, 0,
                             0, (facteur_p1_u * facteur_p2_v + facteur_p2_u * facteur_p1_v)/(det * det), 0,
                             0, 0, (facteur_p1_u * facteur_p2_v + facteur_p2_u * facteur_p1_v)/(det * det);
        del2C_delp1_delp2 *= cl->getA();

        MatrixXd delf_delp = MatrixXd(3,3);

        delf_delp = delC_delp1 * delC_delp2.transpose() + del2C_delp1_delp2 * C;

        delf_delp *= - k_s;

        return delf_delp;
    }

    void addDampingForcesStretch(vector<Particle>& particles){

        for (int face = 0; face < cl->getF().rows(); face++){

            Vector3i triangle = cl->getF().row(face);
            int i = triangle[0], j = triangle[1], k = triangle[2];

            double det, facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v, facteur_k_u, facteur_k_v;
            Vector3d w_u, w_v, w_u_point, w_v_point;

            double bu = 1., bv = 1.;

            calculate_utilities_triangle(particles, i, j, k, w_u, w_v, w_u_point, w_v_point, det, facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v,
                                         facteur_k_u, facteur_k_v);

            // Calculate forces

            Vector3d fi, fj, fk;

            fi = calculate_fp_damping_stretch(w_u, w_v, w_u_point, w_v_point, det, bu, bv, facteur_i_u, facteur_i_v);
            particles[i].addForce(fi);


            fj = calculate_fp_damping_stretch(w_u, w_v, w_u_point, w_v_point, det, bu, bv, facteur_j_u, facteur_j_v);
            particles[j].addForce(fj);


            fk = calculate_fp_damping_stretch(w_u, w_v, w_u_point, w_v_point, det, bu, bv, facteur_k_u, facteur_k_v);
            particles[k].addForce(fk);

            // Calculate force derivatives (position)

            MatrixXd del_fi_del_i, del_fi_del_j, del_fi_del_k;
            MatrixXd del_fj_del_i, del_fj_del_j, del_fj_del_k;
            MatrixXd del_fk_del_i, del_fk_del_j, del_fk_del_k;

            del_fi_del_i = - k_s * second_matrix_delf_stretching(w_u, w_v, w_u_point, w_v_point, det, bu, bv,
                                                         facteur_i_u, facteur_i_v, facteur_i_u, facteur_i_v, 1);
            particles[i].addForceDerivativePosition(i, del_fi_del_i);


            del_fi_del_j = - k_s * second_matrix_delf_stretching(w_u, w_v, w_u_point, w_v_point, det, bu, bv,
                                                         facteur_i_u, facteur_i_v, facteur_j_u, facteur_j_v, 1);
            particles[j].addForceDerivativePosition(i, del_fi_del_j);


            del_fi_del_k = - k_s * second_matrix_delf_stretching(w_u, w_v, w_u_point, w_v_point, det, bu, bv,
                                                         facteur_i_u, facteur_i_v, facteur_k_u, facteur_k_v, 1);
            particles[k].addForceDerivativePosition(i, del_fi_del_k);


            del_fj_del_i = del_fi_del_j.transpose();

            particles[i].addForceDerivativePosition(j, del_fj_del_i);


            del_fj_del_j = - k_s * second_matrix_delf_stretching(w_u, w_v, w_u_point, w_v_point, det, bu, bv,
                                                         facteur_j_u, facteur_j_v, facteur_j_u, facteur_j_v, 1);
            particles[j].addForceDerivativePosition(j, del_fj_del_j);


            del_fj_del_k = - k_s * second_matrix_delf_stretching(w_u, w_v, w_u_point, w_v_point, det, bu, bv,
                                                         facteur_j_u, facteur_j_v, facteur_k_u, facteur_k_v, 1);
            particles[k].addForceDerivativePosition(j, del_fj_del_k);


            del_fk_del_i = del_fi_del_k.transpose();

            particles[i].addForceDerivativePosition(k, del_fk_del_i);


            del_fk_del_j = del_fj_del_k.transpose();

            particles[j].addForceDerivativePosition(k, del_fk_del_j);


            del_fk_del_k = - k_s * second_matrix_delf_stretching(w_u, w_v, w_u_point, w_v_point, det, bu, bv,
                                                         facteur_k_u, facteur_k_v, facteur_k_u, facteur_k_v, 1);
            particles[k].addForceDerivativePosition(k, del_fk_del_k);

            // Calculate force derivatives (velocity)

            MatrixXd del_fi_del_vi, del_fi_del_vj, del_fi_del_vk;
            MatrixXd del_fj_del_vi, del_fj_del_vj, del_fj_del_vk;
            MatrixXd del_fk_del_vi, del_fk_del_vj, del_fk_del_vk;

            del_fi_del_vi = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_i_u, facteur_i_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_i_u, facteur_i_v).transpose();
            particles[i].addForceDerivativeVelocity(i, del_fi_del_vi);

            del_fi_del_vj = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_i_u, facteur_i_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_j_u, facteur_j_v).transpose();
            particles[j].addForceDerivativeVelocity(i, del_fi_del_vj);

            del_fi_del_vk = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_i_u, facteur_i_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_k_u, facteur_k_v).transpose();
            particles[k].addForceDerivativeVelocity(i, del_fi_del_vk);


            del_fj_del_vi = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_j_u, facteur_j_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_i_u, facteur_i_v).transpose();
            particles[i].addForceDerivativeVelocity(j, del_fj_del_vi);

            del_fj_del_vj = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_j_u, facteur_j_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_j_u, facteur_j_v).transpose();
            particles[j].addForceDerivativeVelocity(j, del_fj_del_vj);

            del_fj_del_vk = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_j_u, facteur_j_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_k_u, facteur_k_v).transpose();
            particles[k].addForceDerivativeVelocity(j, del_fj_del_vk);


            del_fk_del_vi = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_k_u, facteur_k_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_i_u, facteur_i_v).transpose();
            particles[i].addForceDerivativeVelocity(k, del_fk_del_vi);

            del_fk_del_vj = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_k_u, facteur_k_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_j_u, facteur_j_v).transpose();
            particles[j].addForceDerivativeVelocity(k, del_fk_del_vj);

            del_fk_del_vk = - k_s * calculate_delC_delp_stretch(w_u, w_v, det, facteur_k_u, facteur_k_v) *
                            calculate_delC_delp_stretch(w_u, w_v, det, facteur_k_u, facteur_k_v).transpose();
            particles[k].addForceDerivativeVelocity(k, del_fk_del_vk);
        }
    }

    Vector3d calculate_fp_damping_stretch(Vector3d w_u, Vector3d w_v, Vector3d w_u_point, Vector3d w_v_point, double det, double bu, double bv,
                                  double facteur_p_u, double facteur_p_v){

        MatrixXd C_point = calculate_C_point_stretch(w_u, w_v, w_u_point, w_v_point);
        MatrixXd delC_delp = calculate_delC_delp_stretch(w_u, w_v, det, facteur_p_u, facteur_p_v);
        Vector3d fp = - k_s * delC_delp * C_point;

        return fp;
    }

    void calculate_utilities_triangle(vector<Particle>& particles, int i, int j, int k, Vector3d& w_u, Vector3d& w_v,
                                      Vector3d& w_u_point, Vector3d& w_v_point, double& det,
                                      double& facteur_i_u, double& facteur_i_v, double& facteur_j_u,
                                      double& facteur_j_v, double& facteur_k_u, double& facteur_k_v){

        Vector3d pi = cl->getV().row(i);
        Vector3d pj = cl->getV().row(j);
        Vector3d pk = cl->getV().row(k);

        Vector3d vi = particles[i].getVelocity();
        Vector3d vj = particles[j].getVelocity();
        Vector3d vk = particles[k].getVelocity();

        double i_u = particles[i].u(), i_v = particles[i].v();
        double j_u = particles[j].u(), j_v = particles[j].v();
        double k_u = particles[k].u(), k_v = particles[k].v();

        double delta_u_1 = j_u - i_u;
        double delta_u_2 = k_u - i_u;
        double delta_v_1 = j_v - i_v;
        double delta_v_2 = k_v - i_v;

        Vector3d delta_p_1 = pj - pi;
        Vector3d delta_p_2 = pk - pi;
        MatrixXd delta_ps = MatrixXd(3,2);
        delta_ps.block<3,1>(0,0) = delta_p_1;
        delta_ps.block<3,1>(0,1) = delta_p_2;

        Vector3d delta_vel_1 = vj - vi;
        Vector3d delta_vel_2 = vk - vi;
        MatrixXd delta_vels = MatrixXd(3,2);
        delta_vels.block<3,1>(0,0) = delta_vel_1;
        delta_vels.block<3,1>(0,1) = delta_vel_2;

        det = delta_u_1 * delta_v_2 - delta_u_2 * delta_v_1;
        MatrixXd inv_right = MatrixXd(2,2);
        inv_right << delta_v_2, -delta_u_2,
                    -delta_v_1,  delta_u_1;
        inv_right *= 1./det;

        MatrixXd ws = delta_ps * inv_right;
        w_u = ws.block<3,1>(0,0);
        w_v = ws.block<3,1>(0,1);

        MatrixXd ws_point = delta_vels * inv_right;
        w_u_point = ws_point.block<3,1>(0,0);
        w_v_point = ws_point.block<3,1>(0,1);

        facteur_i_u = delta_v_1 - delta_v_2;
        facteur_i_v = delta_u_2 - delta_u_1;
        facteur_j_u = delta_v_2;
        facteur_j_v = - delta_u_2;
        facteur_k_u = - delta_v_1;
        facteur_k_v = delta_u_1;
    }

    void calculate_M(vector<Particle>& particles){
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
            p_zero.block<3,1>(3*i,0) = cl->getV().row(i);
        }
    }

    void calculate_v_zero(vector<Particle>& particles){
        v_zero = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d vi = particles[i].getVelocity();
            v_zero.block<3,1>(3*i,0) = vi;
        }
    }

    void calculate_f_zero(vector<Particle>& particles){
        f_zero = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d fi = particles[i].getForce();
            f_zero.block<3,1>(3*i,0) = fi;
        }
    }

    void calculate_del_f_del_p(vector<Particle>& particles){
        del_f_del_p = MatrixXd::Zero(3*n, 3*n);
        for (int j = 0; j < n; j++) {
            vector<MatrixXd> del_f_del_pj = particles[j].getForceDerivativePosition();
            for (int i = 0; i < n; i++) {
                del_f_del_p.block<3,3>(3*i, 3*j) = del_f_del_pj[i];
            }
        }

    }

    void calculate_del_f_del_v(vector<Particle>& particles){
        del_f_del_v = MatrixXd::Zero(3*n, 3*n);
        for (int j = 0; j < n; j++) {
            vector<MatrixXd> del_f_del_vj = particles[j].getForceDerivativeVelocity();
            for (int i = 0; i < n; i++) {
                del_f_del_v.block<3,3>(3*i, 3*j) = del_f_del_vj[i];
            }
        }
    }

    void calculate_y_correction(vector<Particle>& particles){
        y_correction = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d yi = particles[i].y_correc();
            y_correction.block<3,1>(3*i,0) = yi;
        }
    }

    void calculate_z_correction(vector<Particle>& particles){
        z_correction = MatrixXd::Zero(3*n, 1);
        for (int i = 0; i < n; i++) {
            Vector3d zi = particles[i].z_correc();
            z_correction.block<3,1>(3*i,0) = zi;
        }
    }

    void modified_pcg(vector<Particle>& particles){
        // Initialize variables
        MatrixXd A, b, P, P_inverse, r, c, q, s;
        double delta_zero, delta_old, delta_new, alfa, eps = 0.1;

        A = M - h * del_f_del_v - h * h * del_f_del_p;
        b = h*(f_zero + h * del_f_del_p * v_zero + del_f_del_p * y_correction);
        calculate_P_and_P_inverse(P, P_inverse, A);

        // Procedure
        delta_v = z_correction;
        delta_zero = (filter(b, particles).transpose() * P * filter(b, particles)).determinant();
        r = filter(b - A * delta_v, particles);
        c = filter(P_inverse * r, particles);
        delta_new = (r.transpose() * c).determinant();
        while (delta_new > eps * eps * delta_zero) {
            q = filter(A * c, particles);
            alfa = delta_new/((c.transpose() * q).determinant());
            delta_v += alfa * c;
            r -= alfa * q;
            s = P_inverse * r;
            delta_old = delta_new;
            delta_new = (r.transpose() * s).determinant();
            c = filter(s + (delta_new/delta_old) * c, particles);
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

    MatrixXd filter(MatrixXd a, vector<Particle>& particles){
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

    void updateVelocities(vector<Particle>& particles){
        MatrixXd v = v_zero + delta_v;
        for (int i = 0; i < n; i++) {
            particles[i].setVelocity(v.block<3,1>(3*i,0));
        }
    }

    void updatePositions(){
        delta_p = h * (v_zero + delta_v) + y_correction;
        MatrixXd p = p_zero + delta_p;
        for (int i = 0; i < n; i++) {
            cl->getV().row(i) = p.block<3,1>(3*i,0);
        }
    }


    int n;
    double h, g, k_s;
    Cloth* cl;

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
