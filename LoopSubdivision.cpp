/* ========================================================================= *
 *                                                                           *
 *                       Luca Castelli Aleardi                       		 *
 *           Copyright (c) 2019, Ecole Polytechnique                		 *
 *           Department of Computer Science                  				 *
 *                          All rights reserved.                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of the course material developed for		             *
 *   INF574 Digital Representation and Analysis of Shapes (2019/20)			 *
 * ========================================================================= */
#include <igl/opengl/glfw/Viewer.h>

#ifndef HALFEDGE_DS_HEADER
  #define HALFEDGE_DS_HEADER
  #include "HalfedgeDS.cpp"
#endif

using namespace Eigen;
using namespace std;

/**
 * @author Luca Castelli Aleardi (2019)
 */
class LoopSubdivision
{

public:
    /** 
	 * Initialize the data structures
	 **/
    LoopSubdivision(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
    {
        he = &mesh;
        V = &V_original;
        // F = &F_original; // NOT NEEDED if using the half-edge data structure
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = V_original.rows();         // number of vertices in the original mesh

        nVertices = n + e;
        nFaces = F_original.rows()*4;

        V1 = new MatrixXd;
        F1 = new MatrixXi;

        V1->setZero(nVertices, 3);
        F1->setZero(nFaces, 3);
    }

    /** 
	 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
	 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1' 
	 **/
    void subdivide()
    {
        std::cout << "Performing one round subdivision" << endl;
        int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
        int n = he->sizeOfVertices();      // number of vertices in the original mesh
        int F = he->sizeOfFaces();         // number of vertices in the original mesh

        // TO BE COMPLETED
        // first step: compute new midpoint vertices and assign a number, between 0..e-1, to all halfedges
        vector<int> halfEdgeMiddle(2*e, -1);

        for (int i = 0, h = 0; i < e; i++, h++) {

            while (halfEdgeMiddle[h] >= 0)    // Find a halfedge whose middle point hasn't been defined yet
                h++;

            V1->row(i) = computeEdgePoint(h);         // Calculate the middle point
            halfEdgeMiddle[h] = i;                    // Relate it with this halfedge
            halfEdgeMiddle[he->getOpposite(h)] = i;   // and its opposite halfedge
        }

        // second step: update the original points
        for (int i = 0; i < n; i++) {
            V1->row(e+i) = updateOriginalPoint(i);
        }

        // third step: set the face/vertex incidence relations
        for (int f = 0; f < F; f++) {

            int h = he->getEdgeInFace(f);
            int h_next = he->getNext(h);
            int h_prev = he->getPrev(h);
            int middle1 = halfEdgeMiddle[h], middle2 =  halfEdgeMiddle[h_next], middle3 = halfEdgeMiddle[h_prev];

            // New face 1
            F1->row(4*f) = Vector3i(middle1, he->getTarget(h)+e, middle2);

            // New face 2
            F1->row(4*f + 1) = Vector3i(middle2, he->getTarget(h_next)+e, middle3);

            // New face 3
            F1->row(4*f + 2) = Vector3i(middle3, he->getTarget(h_prev)+e, middle1);

            // New face 4
            F1->row(4*f + 3) = Vector3i(middle1, middle2, middle3);
        }

    }

    /** 
	 * Return the number of half-edges
	 **/
    MatrixXd getVertexCoordinates()
    {
        return *V1;
    }

    /** 
	 * Return the number of faces
	 **/
    MatrixXi getFaces()
    {
        return *F1;
    }

    /** 
	 * Print the combinatorial information of the subdivided mesh <b>
	 * verbosity=0: print only the number of vertices and faces <b>
	 * verbosity=1: print all incidence relations
	 **/
    void print(int verbosity)
    {
        cout << "\tn=" << nVertices << ", f=" << nFaces << endl;

        if (verbosity > 0) // print all vertex coordinates and face/vertex incidence relations
        {
            for (int i = 0; i < nVertices; i++)
            {
                cout << "v" << i << ": " << V1->row(i) << endl;
            }

            std::cout << "new faces: " << nFaces << endl;
            for (int i = 0; i < nFaces; i++)
            {
                cout << "f" << i << ": " << F1->row(i) << endl;
            }
        }
    }

private:
    /**
	 * Compute the midpoint of the given half-edge 'h=(u,v)'
	 */
    MatrixXd computeEdgePoint(int h)
    {
        int opposite = he->getOpposite(h);
        int next = he->getNext(h);
        int next_opposite = he->getNext(opposite);

        VectorXd v2 = V->row(he->getTarget(h));
        VectorXd v0 = V->row(he->getTarget(opposite));
        VectorXd v1 = V->row(he->getTarget(next_opposite));
        VectorXd v4 = V->row(he->getTarget(next));

        MatrixXd projectedMiddle = ((v0 + v2) * 3./8 + (v1 + v4) * 1./8).transpose();
        return projectedMiddle;
    }

    /**
	 * Given a vertex 'v' of the original mesh, compute and return its new coordinates
	 */
    MatrixXd updateOriginalPoint(int v)
    {
        int degree = vertexDegree(v);
        double alpha = (degree == 3) ? 3./16 : 3./(8 * degree);

        int firstIncidentEdge = he->getEdge(v);
        int oppositeEdge;
        int nextIncidentCCW = firstIncidentEdge;

        MatrixXd new_v;
        MatrixXd sum = MatrixXd::Zero(1,3);

        do {
            oppositeEdge = he->getOpposite(nextIncidentCCW);

            sum += V->row(he->getTarget(oppositeEdge));

            nextIncidentCCW = he->getPrev(oppositeEdge);

        } while (nextIncidentCCW != firstIncidentEdge);

        new_v = (1 - alpha * degree) * V->row(v) + alpha * sum;

        return new_v;
    }

    int vertexDegree(int v)
    {
        int firstIncidentEdge = he->getEdge(v);
        int degree = 0;
        int oppositeEdge;
        int nextIncidentCCW = firstIncidentEdge;

        do {

            if (he->getFace(nextIncidentCCW) != -1) {
                degree++;
            }

            oppositeEdge = he->getOpposite(nextIncidentCCW);
            nextIncidentCCW = he->getPrev(oppositeEdge);

        } while (nextIncidentCCW != firstIncidentEdge);

        return degree;
    }

    /** Half-edge representation of the original input mesh */
    HalfedgeDS *he;
    MatrixXd *V; // vertex coordinates of the original input mesh

	/** faces/vertex incidence relations in the original mesh */
    // MatrixXi *F; // REMARK: not needed if using the half-edge data structure

    int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
    MatrixXd *V1;          // vertex coordinates of the new subdivided mesh
    MatrixXi *F1;          // faces of the new subdivided mesh
};
