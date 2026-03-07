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
    fy = rho2 * n2 * sin_theta;
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

    double detR = determinant(R);

    if (detR < 0) {
        rho = -rho;
        r1 = -r1;   // negate all three rows
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
    std::cout << "\nTODO: implement the 'calibration()' function in the file 'Calibration/calibration_method.cpp'\n\n";
    // these are just print statements in C++ language
    std::cout << "[Liangliang]:\n"
                 "\tIn this assignment, two essential data structures, 'Matrix' and 'Vector', are provided for the\n"
                 "\tmanipulation and storage of matrices and vectors. These data structures are defined in:\n"
                 "\t    - Calibration/matrix.h: handles matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/vector.h: manages vectors of arbitrary sizes and related functions.\n"
                 "\tCamera calibration requires computing the SVD and inverse of matrices. These functions, along\n"
                 "\twith several other relevant ones, are provided in:\n"
                 "\t    - Calibration/matrix_algo.h: contains functions for determinant, inverse, SVD, linear least-squares...\n"
                 "\tIn the 'Calibration::calibration(...)' function, code snippets are provided for your reference.\n"
                 "\tFor more details about these data structures and a complete list of related functions, please\n"
                 "\trefer to the header files mentioned above.\n\n"
                 "\tFor your final submission, adhere to the following guidelines:\n"
                 "\t    - submit ONLY the 'Calibration/calibration_method.cpp' file.\n"
                 "\t    - remove ALL unrelated test code, debugging code, and comments.\n"
                 "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" << std::flush;
    // std:flush tells computer not to wait for buffer to fill but simply to print as soon as ready.

    /// Below are a few examples showing some useful data structures and functions.
    ///
    /// THE FOLLOWING IS EFFECTIVELY "C++ NUMPY"..............................................

    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. Works like python list. In case you may want to learn more about 'std::vector'
    // check here: https://en.cppreference.com/w/cpp/container/vector
    std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
    array.push_back(5); // append 5 to the array (so the size will increase by 1).
    array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).

    /// To access the value of an element.
    double a = array[2];

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D c(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = c.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// the length of a vector
    double len = p.length();
    /// the squared length of a vector
    double sqr_len = p.length2();

    /// the dot product of two vectors
    double dot_prod = dot(p, q);

    /// the cross product of two vectors
    Vector cross_prod = cross(c, q);

    /// normalize this vector
    cross_prod.normalize();

    // Define an m-by-n double valued matrix.
    // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
    const int m = 6, n = 5;
    Matrix A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
    std::cout << "M: \n" << A << std::endl;

    /// define a 3 by 4 matrix (and all elements initialized to 0.0)
    Matrix M(3, 4, 0.0);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 B;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 P(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

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

    Matrix U(m, m, 0.0);   // initialized with 0s
    Matrix S(m, n, 0.0);   // initialized with 0s
    Matrix V(n, n, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    svd_decompose(A, U, S, V);

    // Now let's check if the SVD result is correct

    // Check 1: U is orthogonal, so U * U^T must be identity
    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
    std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;

    // Compute the inverse of a matrix
    Matrix invT;
    inverse(T, invT);
    // Let's check if the inverse is correct
    std::cout << "T * invT: \n" << T * invT << std::endl;

    // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
    //       final submission.

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    std::cout << "\n[Liangliang]:\n"
                 "\tThis function takes two arrays as input parameters:\n"
                 "\t\t- points_3d: An array of 3D points representing the scene\n"
                 "\t\t- points_2d: An array of 2D image points corresponding to the 3D points\n"
                 "\tThe function should return either 'true' upon successful calibration or 'false' otherwise.\n"
                 "\tUpon success, the following parameters must be stored in the specified variables:\n"
                 "\t\t- fx and fy: focal lengths along the x and y axes, respectively\n"
                 "\t\t- cx and cy: coordinates of the principal point\n"
                 "\t\t- s: the skew factor, i.e., s = -alpha * cot(theta)\n"
                 "\t\t- R: the 3x3 rotation matrix encoding camera orientation\n"
                 "\t\t- t: a 3D vector encoding camera location.\n"
                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)

    if (points_3d.size() != points_2d.size()) { // .size() is a method of std::vector
        std::cerr << "Error: Number of 3D points and 2D points must match!" << std::endl; // flushes and starts new line
        return false; // Exit the function early with failure, flushing guarantees that error reaches console if program crashes
    }

    if (points_3d.size() < 6) {
        std::cerr << "Error: At least 6 point pairs are required for calibration." << std::endl;
        return false;
    }

    // TODO: construct the P matrix (so P * m = 0).

    Matrix P_matrix = construct_P_matrix(points_3d, points_2d); // Auxiliary matrix made up of point data

    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    Matrix34 M_matrix = solve_M_from_P(P_matrix); // the camera parameter matrix!

    // Check if the matrix projects 3D point to pixel coordinate! Can be due to bad measurements or degenerate sample
    double threshold = 10.0;

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

    double MSE = SSE / static_cast<double>(points_3d.size());

    // TODO: extract intrinsic parameters from M.

    // 1. Declare a temporary variable for rho (needed for extrinsics)

    double rho = 0.0;
    if (!intrinsics_from_M(M_matrix, fx, fy, cx, cy, s, rho)) {
        std::cerr << "Error: Failed to decompose intrinsic parameters." << std::endl;
        return false;
    }

    // TODO: extract extrinsic parameters from M.


    if (!extrinsics_from_M(M_matrix, fx, fy, cx, cy, s, rho, R, t)) {
        std::cerr << "Error: Failed to decompose extrinsic parameters." << std::endl;
        return false;
    }

    // TODO: make sure the recovered parameters are passed to the corresponding variables (fx, fy, cx, cy, s, R, and t)

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


















