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


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */

//error evaluation
double calculate_error(const std::vector<Vector3D> &points_3d,const std::vector<Vector2D> &x0,const std::vector<Vector2D> &x1,const Matrix& final_P,const Matrix& P){
    double error = 0;
    for(int i = 0;i<points_3d.size();++i){
        Vector4D X1_4D(points_3d[i][0],points_3d[i][1],points_3d[i][2],1);
        Vector3D x1_hat3 = final_P * X1_4D;
        x1_hat3 = x1_hat3 / x1_hat3[2];
        Vector2D x1_hat(x1_hat3[0],x1_hat3[1]);
        Vector3D x0_hat3 = P * X1_4D;
        x0_hat3 = x0_hat3 / x0_hat3[2];
        Vector2D x0_hat(x0_hat3[0],x0_hat3[1]);
        error+= (x1_hat-x1[i]).norm()+(x0_hat-x0[i]).norm();
    }
    return error;
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
    Matrix44 U;
    Matrix S(4, 4, 0.0);
    Matrix44 V;
    svd_decompose(A, U, S, V);
    Vector4D X_4D = V.get_column(3);
    Vector3D X;
    X[0] = X_4D[0]/ X_4D[3];
    X[1] = X_4D[1]/ X_4D[3];
    X[2] = X_4D[2]/ X_4D[3];

    return X;
}

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
    if (points_0.size() != points_1.size()&&points_0.size()<8||points_1.size()<8)
    {
        std::cout << "The number of points in two images are not equal!" << std::endl;
        return false;
    }

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
        points1_normalized[i] = T1 * points0_tmp;
        points2_normalized[i] = T2 * points1_tmp;
    }
    // construct matrix W, Wf = 0,solve the equation using SVD, F is the last column of V, we recover fundamental matrix from F,
    int num_points = points_0.size();
    Matrix W(num_points, 9, 0.0);
    for (int i = 0; i < num_points; i++)
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
    // solve the equation using SVD   // w = n*9, U = n*n, S = n*9, V = 9*9
    Matrix U(num_points, num_points, 0.0);
    Matrix S(num_points, 9, 0.0);
    Matrix V(9, 9, 0.0);

    svd_decompose(W,U,S,V);

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
    Matrix33 F_denormalized = T2.transpose() * F_f * T1;
    // obtained the fundamental matrix, now we recover the relative pose
    //first we compute the essential matrix
    // compute the intrinsic matrix, using given focal length and principle point, and skewness assume to be 0
    Matrix33 K1(0.0);
    Matrix33 K2(0.0);
    // fx, fy, cx, cy is given
    K1(0, 0) = fx; K1(1, 1) = fy; K1(0, 2) = cx; K1(1, 2) = cy; K1(2, 2) = 1;
    K2(0, 0) = fx;K2(1, 1) = fy;K2(0, 2) = cx; K2(1, 2) = cy; K2(2, 2) = 1;


    Matrix33 E = K2.transpose() * F_denormalized * K1;
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
    R1 =  det_R1*R1;
    R2 = det_R2*R2 ;
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


    // compute the 3D points

    // linear method to triangulate the 3D points AP = 0, where A is 4*4 matrix, P is 4*1 matrix
    for (int i = 0; i < points_0.size(); i++) {
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

    // choose the correct rotation and translation matrix
    if (!(positive_count1 < positive_count2) && !(positive_count1 <positive_count3) && !(positive_count1 < positive_count4)) {
        R = R1;
        t = t1;
        final_P = P1;
    }
    else if (!(positive_count2 < positive_count1) && !(positive_count2 < positive_count3) && !(positive_count2 < positive_count4)) {
        R = R1;
        t = t2;
        final_P = P2;
    }
    else if (!(positive_count3 < positive_count1) && !(positive_count3 < positive_count2) && !(positive_count3 < positive_count4)) {
        R = R2;
        t = t1;
        final_P = P3; // wait why it is p2 but not p3
    }
    else {
        R = R2;
        t = t2;
        final_P = P4;
    }
    // compute the 3D points, recover the 3D points from the 2D points
    std::vector<Vector2D> x0_pts;
    std::vector<Vector2D> x1_pts;
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
        points_3d.push_back(X1);
        x0_pts.push_back(x0);
        x1_pts.push_back(x1);
    }
    double error_evaluation = calculate_error(points_3d,x0_pts,x1_pts,final_P,P);
    std::cout<< "error_result "<<error_evaluation<<std::endl;
    return points_3d.size() > 0;
}
