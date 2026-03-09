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

#include "calibration.h"
#include "matrix_algo.h" // h-files are headers with underlying .cpp source code


using namespace easy3d; // checks if unknown functions are included in easy3d library, used mostly in .cpp

/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by fx, fy, cx, cy, skew, R, and t).
 */

// =============================== HELPER FUNCTIONS ====================================


Matrix construct_P_matrix(const std::vector<Vector3D>& points_3d, // "const...&" implies read-only and don't copy, for speed
                          const std::vector<Vector2D>& points_2d) {

    int n_points = points_3d.size();
    Matrix P(2 * n_points, 12, 0.0);    // empty matrix

    for (int i = 0; i < n_points; i++) {    // start, end, iteration i++ means i = i+1
        const Vector3D& wc = points_3d[i];  // world coordinate
        const Vector2D& pc = points_2d[i];  // pixel coordinate

        double u = pc.x();                  // for readability
        double v = pc.y();
        double X = wc.x();
        double Y = wc.y();
        double Z = wc.z();

        // .x(), .y(), and .z() are methods for Vector2D and Vector3D

        P.set_row(2 * i, {              // populate every even row
            X, Y, Z, 1.0, 0.0, 0.0, 0.0, 0.0, -u * X, -u * Y, -u * Z, -u
        });

        P.set_row(2 * i + 1, {          // populate every odd row
            0.0, 0.0, 0.0, 0.0, X, Y, Z, 1.0, -v * X, -v * Y, -v * Z, -v
        });
    }
    return P;
}

Matrix34 solve_M_from_P(const Matrix& P) {
    const int r = P.rows(); // number of coordinate pairs!
    const int c = 12;

    // Dimensions of matrices in SVD:
    Matrix U(r, r, 0.0);
    Matrix S(r, c, 0.0);
    Matrix V(c, c, 0.0);

    svd_decompose(P, U, S, V);

    // the last column (the solution vector m)
    // columns - 1 is index 11
    Vector m = V.get_column(11);

    // From LiangLiang's helper notes! Much faster.
    return Matrix34(m[0], m[1], m[2], m[3],
                    m[4], m[5], m[6], m[7],
                    m[8], m[9], m[10], m[11]);
}

bool intrinsics_from_M(const Matrix34& M,
                                    double& fx, double& fy, double& cx, double& cy, double& s,
                                    double& rho)
{

    // Vector3D constructor
    Vector3D a1(M(0,0), M(0,1), M(0,2));
    Vector3D a2(M(1,0), M(1,1), M(1,2));
    Vector3D a3(M(2,0), M(2,1), M(2,2));

    // Tool: .length() and .length2() (squared length)
    double len3_sq = a3.length2();  // .length2() gives squared magnitude
    if (len3_sq < 1e-12) return false;

    rho = 1.0 / std::sqrt(len3_sq);
    double rho2 = rho * rho;        // rho squared

    cx = rho2 * dot(a1, a3);
    cy = rho2 * dot(a2, a3);

    Vector3D a1xa3 = cross(a1, a3);
    Vector3D a2xa3 = cross(a2, a3);

    double n1 = a1xa3.length();     // magnitude of normal to a1 and a3
    double n2 = a2xa3.length();     // magnitude of normal to a2 and a3
    if (n1 < 1e-12 || n2 < 1e-12) return false;     // degeneracy check

    double cos_theta = -dot(a1xa3, a2xa3) / (n1 * n2);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta)); // Clamp for safety

    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);      // trigonometric identity
    if (sin_theta < 1e-12) return false;                            // degeneracy check

    fx = rho2 * n1 * sin_theta;
    fy = rho2 * n2;
    s  = -fx * (cos_theta / sin_theta);

    return true;
}

bool extrinsics_from_M(const Matrix34& M,
                                    double fx, double fy, double cx, double cy, double s, double rho,
                                    Matrix33& R, Vector3D& t)
{
    // b is the 4th column of M (the translation-related part)
    Vector3D b(M(0,3), M(1,3), M(2,3));

    // We need a2 and a3 again to compute r1
    Vector3D a2(M(1,0), M(1,1), M(1,2));
    Vector3D a3(M(2,0), M(2,1), M(2,2));

    // Calculate r1 (row 1 of R)
    // Tool: cross() and .normalize()
    Vector3D r1 = cross(a2, a3);
    r1.normalize();

    // Calculate r3 and r2
    Vector3D r3 = a3;
    r3.normalize(); // r3 = rho * a3

    // r2 is the cross product of r3 and r1 (ensures orthogonality)
    Vector3D r2 = cross(r3, r1);

    // Build R matrix
    // Tool: Matrix33 constructor using the components of our vectors
    R = Matrix33(r1.x(), r1.y(), r1.z(),
                 r2.x(), r2.y(), r2.z(),
                 r3.x(), r3.y(), r3.z());

    // double detR = R(0,0)*(R(1,1)*R(2,2) - R(1,2)*R(2,1))
    //         - R(0,1)*(R(1,0)*R(2,2) - R(1,2)*R(2,0))
    //         + R(0,2)*(R(1,0)*R(2,1) - R(1,1)*R(2,0));

     if (b.z() < 0) {
        rho = - rho;
    //    r1 = r1;   // this row is unaffected by the sign of rho
        r2 = -r2;
        r3 = -r3;
        R = Matrix33(r1.x(), r1.y(), r1.z(),
                     r2.x(), r2.y(), r2.z(),
                     r3.x(), r3.y(), r3.z());
    }

    // Solve for t: t = rho * K_inv * b
    // Using back-substitution (K is upper triangular)
    double t3 = b.z();
    double t2 = (b.y() - cy * t3) / fy;
    double t1 = (b.x() - s * t2 - cx * t3) / fx;

    // Final Translation Vector
    t = Vector3D(t1, t2, t3);
    t = rho * t;

    return true;
}


// ============================= MAIN CALIBRATION FUNCTION ============================



bool Calibration::calibration(                  /// Calibration CLASS defined in header, logic included cpp file
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points. const promises read-only!
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,  /// output: focal length (i.e., K[0][0]).
        double& fy,  /// output: focal length (i.e., K[1][1]).
        double& cx,  /// output: x component of the principal point (i.e., K[0][2]).
        double& cy,  /// output: y component of the principal point (i.e., K[1][2]).
        double& s,   /// output: skew factor (i.e., K[0][1]), which is s = -alpha * cot(theta).
        Matrix33& R, /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t) /// output：a 3D vector encoding camera translation.
{
    // 1. check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    if (points_3d.size() != points_2d.size()) { // .size() is a method of std::vector
        std::cerr << "Error: Number of 3D points and 2D points must match!" << std::endl; // flushes and starts new line
        return false; // Exit the function early with failure, flushing guarantees that error reaches console if program crashes
    }

    if (points_3d.size() < 6) {
        std::cerr << "Error: At least 6 point pairs are required for calibration." << std::endl;
        return false;
    }

    // 2. construct the P matrix (so P * m = 0).

    Matrix P_matrix = construct_P_matrix(points_3d, points_2d); // Auxiliary matrix made up of point data

    // 3. solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.

    Matrix34 M_matrix = solve_M_from_P(P_matrix); // the camera parameter matrix!

    // Check if the matrix projects 3D point to pixel coordinate! Can be due to bad measurements or degenerate sample
    double threshold = 10.0; // this is an arbitrary choice. Improvement: relate to dataset

    double SSE = 0;
    for (size_t i = 0; i < points_3d.size(); ++i) {
        Vector3D p_homo = M_matrix * points_3d[i].homogeneous();
        Vector2D p_proj = p_homo.cartesian();
        double error = (p_proj - points_2d[i]).length(); // magnitude of the vector!
        SSE += error * error;
        if (error > threshold) {
            std::cerr << "Calibration failed. Results are invalid." << std::endl;
            return false;
        }
    }

    double MSE = SSE / points_3d.size();

    // 4. extract intrinsic parameters from M.

    double rho = 0.0; // Declare a temporary variable for rho (needed for extrinsics)
    if (!intrinsics_from_M(M_matrix, fx, fy, cx, cy, s, rho)) {
        std::cerr << "Error: Failed to decompose intrinsic parameters." << std::endl;
        return false;
    }

    // 5. extract extrinsic parameters from M.

    if (!extrinsics_from_M(M_matrix, fx, fy, cx, cy, s, rho, R, t)) {
        std::cerr << "Error: Failed to decompose extrinsic parameters." << std::endl;
        return false;
    }

    // 6. make sure the recovered parameters are passed to the corresponding variables (fx, fy, cx, cy, s, R, and t)

    std::cout << "\n=======================================\n";
    std::cout << "      CALIBRATION SUCCESSFUL!          \n";
    std::cout << "=======================================\n";

    std::cout << "\n--- Intrinsic Parameters ---\n";
    std::cout << "Focal length (fx, fy): " << fx << ", " << fy << "\n";
    std::cout << "Principal pt (cx, cy): " << cx << ", " << cy << "\n";
    std::cout << "Skew factor  (s)     : " << s << "\n";

    std::cout << "\n--- Extrinsic Parameters ---\n";
    std::cout << "Rotation Matrix (R):\n" << R << "\n";

    // For the translation vector, we can format it nicely with brackets
    std::cout << "Translation Vector (t):\n";
    std::cout << "[" << t.x() << ", " << t.y() << ", " << t.z() << "]\n";

    std::cout <<"The Mean Squared Error for this dataset is: " << MSE << std::endl;

    std::cout << "=======================================\n\n";

    return true;
}


















