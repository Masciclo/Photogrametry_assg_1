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
#include "matrix_algo.h"


using namespace easy3d;


    // Return double to handle precise coordinates
    std::vector<double> read_file(std::string path) {
        std::ifstream myfile(path);
        std::vector<double> all_data;
        double val; // Changed to double

        if (myfile.is_open()) {
            while (myfile >> val) {
                all_data.push_back(val);
            }
            myfile.close();
        } else {
            std::cerr << "Error: Could not open file at: " << path << std::endl;
        }
        return all_data;
    }

    // Define coord structure to output two vectors in flat_to matrix function
    struct coord {
        std::vector<std::vector<double>> pw; // Changed to double for precision
        std::vector<std::vector<double>> pc;
    };

    coord flat_to_matrix(std::vector<double> flat_data) {
        // 1. Catch the data as doubles

        coord result; // This specify the result format of the data.

        // 2. Loop through the flat vector 5 elements at a time
        for (size_t i = 0; i + 4 < flat_data.size(); i += 5) {
            // Create 3D point {X, Y, Z}
            std::vector<double> point3D = { flat_data[i], flat_data[i+1], flat_data[i+2], 1 };
            result.pw.push_back(point3D);

            // Create 2D point {u, v}
            std::vector<double> point2D = { flat_data[i+3], flat_data[i+4] };
            result.pc.push_back(point2D);
        }

        return result;
    }

    Matrix vector_to_matrix(std::vector<std::vector<double>> pw, std::vector<std::vector<double>> pc) {
        int n_points = pw.size();

        // The DLT Design Matrix 'L' usually has 2*n rows and 11 columns
        Matrix L(2 * n_points, 12, 0.0);

        for (int i = 0; i < n_points; i++) {
            // 1. Extract our point values for readability
            double X = pw[i][0];
            double Y = pw[i][1];
            double Z = pw[i][2];
            double u = pc[i][0];
            double v = pc[i][1];

            // 2. Build the first row (for the u-coordinate)
            // Structure: [ X, Y, Z, 1, 0, 0, 0, 0, -u*X, -u*Y, -u*Z ]
            L.set_row(2 * i, {X, Y, Z, 1.0, 0.0, 0.0, 0.0, 0.0, -u*X, -u*Y, -u*Z, -u});

            // 3. Build the second row (for the v-coordinate)
            // Structure: [ 0, 0, 0, 0, X, Y, Z, 1, -v*X, -v*Y, -v*Z ]
            L.set_row(2 * i + 1, {0.0, 0.0, 0.0, 0.0, X, Y, Z, 1.0, -v*X, -v*Y, -v*Z, -u});
        }

        return L;
    }


    Matrix solve_M_from_P(const Matrix& P) {
            // P should be (2n x 12)
            if (P.cols() != 12) {
                std::cerr << "Expected P to have 12 columns, got " << P.cols() << "\n";
                return Matrix(); // empty matrix as "error"
            }

            const int r = P.rows();
            const int c = P.cols(); // 12

            Matrix U(r, r, 0.0);
            Matrix S(r, c, 0.0);
            Matrix V(c, c, 0.0);

            svd_decompose(P, U, S, V);

            // Solution is last column of V (smallest singular value)
            Vector m = V.get_column(V.cols() - 1); // length 12

            // Reshape m -> 3x4 projection matrix M
            Matrix M(3, 4, 0.0);
            M(0,0) = m[0];   M(0,1) = m[1];   M(0,2) = m[2];   M(0,3) = m[3];
            M(1,0) = m[4];   M(1,1) = m[5];   M(1,2) = m[6];   M(1,3) = m[7];
            M(2,0) = m[8];   M(2,1) = m[9];   M(2,2) = m[10];  M(2,3) = m[11];

            return M;
        }

    bool check_projection_simple(const coord& data, const Matrix& M) {
            double tol = 0.1;

            for (int i = 0; i < (int)data.pw.size(); i++) {
                double X = data.pw[i][0], Y = data.pw[i][1], Z = data.pw[i][2];
                double u = data.pc[i][0], v = data.pc[i][1];

                double x = M(0,0)*X + M(0,1)*Y + M(0,2)*Z + M(0,3);
                double y = M(1,0)*X + M(1,1)*Y + M(1,2)*Z + M(1,3);
                double w = M(2,0)*X + M(2,1)*Y + M(2,2)*Z + M(2,3);

                if (std::abs(w) < 1e-12) return false;   // can't project -> fail

                double up = x / w;
                double vp = y / w;

                if (std::abs(up - u) > tol || std::abs(vp - v) > tol) {
                    return false; // as soon as one point is bad -> fail
                }
            }

            return true; // all points were within tolerance
        }

/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by fx, fy, cx, cy, skew, R, and t).
 */
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
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


    // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
    //       final submission.

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...



    std::string myPath = R"(D:\TU Delft\GEO1016 Photogrametry\Photogrametry_assg_1\A1_Calibration_Code\resources\data\test_data_1(6_points)-normal.txt)";
    std::vector<double> flat_data = read_file(myPath);

    coord coord_result  = flat_to_matrix(flat_data);

    std::vector<std::vector<double>> point3D = coord_result.pw;
    std::cout << "Number of points " << point3D.size() << std::endl;

    std::vector<std::vector<double>> point2D = coord_result.pc;

    Matrix m = vector_to_matrix(coord_result.pw, coord_result.pc);
    std::cout << m << std::endl;

    Matrix M = solve_M_from_P(m);

    bool results = check_projection_simple(coord_result, M);



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

    // TODO: construct the P matrix (so P * m = 0).

    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    // TODO: extract intrinsic parameters from M.

    // TODO: extract extrinsic parameters from M.

    // TODO: make sure the recovered parameters are passed to the corresponding variables (fx, fy, cx, cy, s, R, and t)

    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return results;
}

















