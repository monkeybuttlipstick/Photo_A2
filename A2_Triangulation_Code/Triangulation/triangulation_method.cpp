/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;
Vector3D linear_triangulation(const Vector2D& x1, const Vector2D& x2, const Matrix34& P1, const Matrix34& P2);

/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /*
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'
    

    // TODO: delete all above example code in your final submission
    */
    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...
    // implement 8 points algorithm, first we conduct normalization
    
    if (points_0.size() != points_1.size())
    {
		std::cout << "The number of points in two images are not equal!" << std::endl;
		return false;
	}   // examine whether the number of points in two images are equal
    Vector2D mean1(0, 0);
    Vector2D mean2(0, 0);
    for (int i = 0; i < points_0.size(); i++)
    {
		mean1 += points_0[i];
		mean2 += points_1[i];
	}
    mean1 /= points_0.size();
    mean2 /= points_1.size();
    double mean_distance1 = 0;
    double mean_distance2 = 0;
    for (int i = 0; i < points_0.size(); i++)
    {
        mean_distance1 += (points_0[i] - mean1).norm();
        mean_distance2 += (points_1[i] - mean2).norm();
    }
    mean_distance1 /= points_0.size();
    mean_distance2 /= points_1.size();
    double scale1 = sqrt(2) / mean_distance1;
    double scale2 = sqrt(2) / mean_distance2;
    // construct the transformation matrix
    Matrix33 T1(scale1, 0, -scale1 * mean1[0],
        				0, scale1, -scale1 * mean1[1],
        				0, 0, 1); //transformation matrix for the first image
    Matrix33 T2(scale2, 0, -scale2 * mean2[0],
                0, scale2, -scale2 * mean2[1],
        		0, 0, 1); //transformation matrix for the second image
    // normalize the points
    std::vector<Vector3D> points1_normalized(points_0.size());
    std::vector<Vector3D> points2_normalized(points_1.size());
    for (int i = 0; i < points_0.size(); i++)
    {
        // normalize the points and store them in the vector 3d
        Vector3D points0_tmp(points_0[i][0], points_0[i][1], 1);
        Vector3D points1_tmp(points_1[i][0], points_1[i][1], 1);
        //std::cout << points1_tmp << std::endl;
        //std::cout << points0_tmp << std::endl;
        //std::cout << "**************************************************" << std::endl;
		points1_normalized[i] = T1 * points0_tmp;
        points2_normalized[i] = T2 * points1_tmp;
        //std::cout << points1_normalized[i] << std::endl;
        //std::cout << points2_normalized[i] << std::endl;
        //std::cout << "..............................................." << std::endl;
	}
    // construct matrix W, Wf = 0,solve the equation using SVD, F is the last column of V, we recover fundamental matrix from F,

    Matrix W(points_0.size(), 9, 0.0);
    
    //std::cout <<"------------------------------"<<std::endl;
    int num_rows = points_0.size();
    for (int i = 0; i < num_rows; i++)

    {
        W(i, 0) = points1_normalized[i][0] * points2_normalized[i][0];
        W(i, 1) = points1_normalized[i][1] * points2_normalized[i][0];
        W(i, 2) = points2_normalized[i][0];
        W(i, 3) = points1_normalized[i][0] * points2_normalized[i][1];
        W(i, 4) = points1_normalized[i][1] * points2_normalized[i][1];
        W(i, 5) = points2_normalized[i][1];
        W(i, 6) = points1_normalized[i][0];
        W(i, 7) = points1_normalized[i][1];
        W(i, 8) = 1;
    }
    std::cout<<"W line 195 is correct"<<std::endl;
    // solve the equation using SVD   // w = n*9, U = n*n, S = n*9, V = 9*9

    Matrix U(points_0.size(), points_0.size(), 0.0);
    Matrix S(points_0.size(), 9, 0.0);
    Matrix V(9, 9, 0.0);

    svd_decompose(W,U,S,V);
   
    std::cout<<"S line 204 is correct"<<std::endl;
    // get the last column of V
    Vector F = V.get_column(8);
    // construct the fundamental matrix, recovered from normalized points
    Matrix33 F_normalized(F[0], F[1], F[2],
        				  F[3], F[4], F[5],
        				  F[6], F[7], F[8]);
    // enforce the rank 2 constraint
    Matrix33 U_f(0.0);
    Matrix33 S_f(0.0);
    Matrix33 V_f(0.0);

    svd_decompose(F_normalized, U_f, S_f, V_f);
    S_f(2, 2) = 0;
    Matrix33 F_f = U_f * S_f * V_f.transpose();
    // recover the fundamental matrix, denormalized
    std::cout << F_f << std::endl;
    std::cout << T2 << std::endl;
    std::cout << T1 << std::endl;
    std::cout << "********************above is F_f, t2,t1**********************" << std::endl;
    Matrix33 F_denormalized = T2.transpose() * F_f * T1;
    std::cout<<F_denormalized<<std::endl;
    std::cout<<"********************above is f noralized **********************"<<std::endl;
    // obtained the fundamental matrix, now we recover the relative pose
    //first we compute the essential matrix
    // compute the intrinsic matrix, using given focal length and principle point, and skewness assume to be 0
    Matrix33 K1(0.0);
    Matrix33 K2(0.0);
    // fx, fy, cx, cy is given
    K1(0, 0) = fx; K1(1, 1) = fy; K1(0, 2) = cx; K1(1, 2) = cy; K1(2, 2) = 1;
    K2(0, 0) = fx;K2(1, 1) = fy;K2(0, 2) = cx; K2(1, 2) = cy; K2(2, 2) = 1;

    std::cout << K1 << std::endl;
    std::cout << K2 << std::endl;
    std::cout<<"********************above is K1, K2 **********************"<<std::endl;


    Matrix33 E = K2.transpose() * F_denormalized * K1;
    std::cout << E << std::endl;
    std::cout << "********************above is E **********************" << std::endl;
    // confirm F and E is crrect

    // Compute the SVD of E
    Matrix33 U_e(0.0);
    Matrix33 D_e(0.0);
    Matrix33 V_e(0.0);
    svd_decompose(E, U_e, D_e, V_e);
    // compute the rotation matrix
    Matrix33 W_e(0.0);
    W_e(0, 1) = -1;
    W_e(1, 0) = 1;
    W_e(2, 2) = 1;
    Matrix33 R1 = U_e * W_e * V_e.transpose();
    double det_R1 = determinant(R1);
    // two possible rotation matrix
    Matrix33 R2 = U_e * W_e.transpose() * V_e.transpose();
    double det_R2 = determinant(R2);
    // compute the translation matrix

    R1 =  det_R1* R1;
    R2 = det_R2 * R2;

    // two possible translation matrix
    Vector3D t1 = U_e.get_column(2);
    Vector3D t2 = -U_e.get_column(2);
    // compute the 3D points and judge whether the points are in front of the camera to determine the correct rotation and translation
    // first we need to construct the projection matrix
    std::cout<<"project matrix"<<std::endl;
    Matrix34 P1(0.0);
    Matrix34 P2(0.0);
    Matrix34 P3(0.0);
    Matrix34 P4(0.0);

    P1.set_column(0, R1.get_column(0));
    P1.set_column(1, R1.get_column(1));
    P1.set_column(2, R1.get_column(2));
    P1.set_column(3, t1);
    P2.set_column(0, R1.get_column(0));
    P2.set_column(1, R1.get_column(1));
    P2.set_column(2, R1.get_column(2));
    P2.set_column(3, t2);
    P3.set_column(0, R2.get_column(0));
    P3.set_column(1, R2.get_column(1));
    P3.set_column(2, R2.get_column(2));
    P3.set_column(3, t1);
    P4.set_column(0, R2.get_column(0));
    P4.set_column(1, R2.get_column(1));
    P4.set_column(2, R2.get_column(2));
    P4.set_column(3, t2);
    

    // given K1, K2, P1, P2, P3, P4, we can compute the 3D points
    // first we need to construct the matrix A

    P1 = K2 * P1;
    P2 = K2 * P2;
    P3 = K2 * P3;
    P4 = K2 * P4;
    Matrix34 P(0.0); // matrix P is the projection matrix of camera 1
    P(0,0) = 1; P(1, 1) = 1; P(2, 2) = 1;
    P = K1 * P;

    // compute the 3D points
    int positive_count1 = 0;
    int positive_count2 = 0;
    int positive_count3 = 0;
    int positive_count4 = 0;




    for (int i = 0; i < points_0.size(); i++)   // linear method to triangulate the 3D points AP = 0, where A is 4*4 matrix, P is 4*1 matrix
 {
        Vector2D x1(points_1[i].x(), points_1[i].y());
        Vector2D x0(points_0[i].x(), points_0[i].y());


        // compute the 3D points

       
        // linear method to triangulate the 3D points AP = 0, where A is 4*3 matrix, P is 3*1 matrix, we can use SVD to solve this problem

        Vector3D X1 = linear_triangulation(x0,x1,P,P1);
        Vector3D X2 = linear_triangulation(x0,x1,P,P2);
        Vector3D X3 = linear_triangulation(x0,x1,P,P3);
        Vector3D X4 = linear_triangulation(x0,x1,P,P4);

        // judge whether the points are in front of the camera
        if (X1[2] > 0) {
			positive_count1++;
		}
        if (X2[2] > 0) {
            positive_count2++;
        }
        if (X3[2] > 0) {
			positive_count3++;
		}
        if (X4[2] > 0) {
            positive_count4++;
        }
    }
    Matrix final_P(3, 4, 0.0);
    std::cout<< "positive_couts 1" << positive_count1 << std::endl;
    std::cout<< "positive_couts 2" << positive_count2 << std::endl;
    std::cout<< "positive_couts 3" << positive_count3 << std::endl;
    std::cout<< "positive_couts 4" << positive_count4 << std::endl;
    std::cout<<"projection matrix after"<<std::endl;
    std::cout << P1 << std::endl;
    std::cout << "********************above is P1 **********************" << std::endl;
    std::cout << P2 << std::endl;
    std::cout << "********************above is P2 **********************" << std::endl;
    std::cout << P3 << std::endl;
    std::cout << "********************above is P3 **********************" << std::endl;
    std::cout << P4 << std::endl;
    std::cout << "********************above is P4 **********************" << std::endl;
    std::cout << P << std::endl;
    std::cout << "********************above is P **********************" << std::endl;
    // choose the correct rotation and translation matrix
    if (!(positive_count1 < positive_count2) && !(positive_count1 <positive_count3) && !(positive_count1 < positive_count4)) {
		R = R1;
		t = t1;
        final_P = P1;
		std::cout << "final p1 : " << final_P << std::endl;
    }
    else if (!(positive_count2 < positive_count1) && !(positive_count2 < positive_count3) && !(positive_count2 < positive_count4)) {
		R = R1;
		t = t2;
		final_P = P2;
		std::cout << "final p2 : " << final_P << std::endl;
    }
    else if (!(positive_count3 < positive_count1) && !(positive_count3 < positive_count2) && !(positive_count3 < positive_count4)) {
		R = R2;
		t = t1;
        final_P = P3; // wait why it is p2 but not p3
        std::cout << "final p3 : " << final_P << std::endl;
    }
    else {
		R = R2;
		t = t2;
        final_P = P4;
        std::cout << "final p4 : " << final_P << std::endl;
    }
    std::cout<<"final R "<<R<<std::endl;
    std::cout<<"final t "<<t<<std::endl;
    // test wheter the r and t is correct
    
    // compute the 3D points, recover the 3D points from the 2D points
    for (int i = 0; i < points_0.size(); i++) {
        double pt0_x = points_0[i].x();
        double pt0_y = points_0[i].y();
        double pt1_x = points_1[i].x();
        double pt1_y = points_1[i].y();
        Vector2D x1(pt1_x, pt1_y);
        Vector2D x0(pt0_x, pt0_y);
        // compute the 3D points


        // linear method to triangulate the 3D points AP = 0, where A is 4*3 matrix, P is 3*1 matrix, we can use SVD to solve this problem
        Vector3D X1 = linear_triangulation(x0, x1, P, final_P);
        //std::cout<< X1 << std::endl;

        points_3d.push_back(X1);
        std::cout<<"points_3d "<<points_3d[i]<<std::endl;
       // show the error by examining th 2d points'error 
        Vector4D X1_4D(X1[0],X1[1],X1[2],1);
        Vector3D x1_hat = final_P * X1_4D;
        x1_hat = x1_hat / x1_hat[2];
        Vector3D x2_hat = P * X1_4D;
        x2_hat = x2_hat / x2_hat[2];
        // mean error on corresponding points

        double error1 = (x1_hat[0] - x1[0]);
        double error2 = (x1_hat[1] - x1[1]);
        double error3 = (x2_hat[0] - x0[0]);
        double error4 = (x2_hat[1] - x0[1]);
        std::cout<<"error1 "<<error1<<std::endl;
        std::cout<<"error2 "<<error2<<std::endl;
        std::cout<<"error3 "<<error3<<std::endl;
        std::cout<<"error4 "<<error4<<std::endl;

       
    }






    // TODO: check if the input is valid (always good because you never known how others will call your function).

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}
// linear method
Vector3D linear_triangulation(const Vector2D& x1, const Vector2D& x2, const Matrix34& P1, const Matrix34& P2) {
	// construct the matrix A

	Matrix A(4, 4, 0.0);
  
    Vector row1 = P1.get_row(2);
    Vector row2 = P1.get_row(2);
    Vector row3 = P2.get_row(2);
    Vector row4 = P2.get_row(2);
    
    row1 = row1 * x1[0];
    row2 = row2 * x1[1];
    row3 = row3 * x2[0];
    row4 = row4 * x2[1];
    


	A.set_row(0, row1-P1.get_row(0));
	A.set_row(1, row2-P1.get_row(1));
	A.set_row(2, row3-P2.get_row(0));
	A.set_row(3, row4-P2.get_row(1));

    
	// compute the SVD of A
    //solve_least_squares(const Matrix &A, const std::vector<double> &b, std::vector<double> &x) the method has been introduced
    // remember we calculate svd of a and then use the last column of V to get the 3D points

    //solve_least_squares(A, b, x);
    //SVD svd(A);
    //svd solve the least square problem
    // u,S,V

    Matrix44 U;
    Matrix S(4, 4, 0.0);
    Matrix44 V;
    svd_decompose(A, U, S, V);

    std::vector<double> b = {0,0,0,0};
    std::vector<double> x;
    Vector4D X_4D = V.get_column(3);
    Vector3D X;
    X[0] = X_4D[0]/ X_4D[3];
    X[1] = X_4D[1]/ X_4D[3];
    X[2] = X_4D[2]/ X_4D[3];

    
   
   


	return X;
}
