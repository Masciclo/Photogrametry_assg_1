#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Simplest function: Returns one flat vector of all integers found in the file
std::vector<int> read_file() {
    std::string path = R"(D:\TU Delft\GEO1016 Photogrametry\Assignment_n_1\Photogrametry_assg_1\A1_Calibration_Code\resources\data\test_data_1(6_points)-normal.txt)";
    std::ifstream myfile(path);

    std::vector<int> all_data;
    int val;

    if (myfile.is_open()) {
        // This operator >> automatically skips spaces and newlines
        while (myfile >> val) {
            all_data.push_back(val);
        }
        myfile.close();
    } else {
        std::cerr << "Error: Could not open file at: " << path << std::endl;
    }

    return all_data;
}
#include <iostream>
#include <vector>

struct coord {
    std::vector<std::vector<int>> pw; // 2D: each element is {x, y, z}
    std::vector<std::vector<int>> pc; // 2D: each element is {u, v}
};

int main() {
    // This is the flat vector you just read

        

    coord result;

    // Loop through the flat vector 5 elements at a time
    for (size_t i = 0; i + 4 < flat_data.size(); i += 5) {

        // 1. Create the 3D point (pw)
        std::vector<int> point3D = { flat_data[i], flat_data[i+1], flat_data[i+2] };
        result.pw.push_back(point3D);

        // 2. Create the 2D point (pc)
        std::vector<int> point2D = { flat_data[i+3], flat_data[i+4] };
        result.pc.push_back(point2D);
    }

    // --- Verification Print ---
    std::cout << "Processed " << result.pw.size() << " full points.\n";
    std::cout << "First PW row: " << result.pw[0][0] << " " << result.pw[0][1] << " " << result.pw[0][2] << "\n";
    std::cout << "First PC row: " << result.pc[0][0] << " " << result.pc[0][1] << "\n";

    return 0;
}
