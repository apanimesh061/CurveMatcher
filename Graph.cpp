#include "Graph.h"
#include <boost/filesystem.hpp>

using namespace boost::filesystem;
using namespace std;

template<typename T>
ostream &operator<<(ostream &os, vector<T> vec) {
    os << "[";
    if (vec.size() != 0) {
        copy(vec.begin(), vec.end() - 1, ostream_iterator<T>(os, ", "));
        os << vec.back();
    }
    os << "]";
    return os;
}

/**
 * Using Boost library to read all entries in a directory
 * @param dir_path
 * @return
 */
vector<string> get_files(string dir_path) {
    path p(dir_path);
    try {
        if (exists(p)) {
            if (is_directory(p)) {
                vector<string> result;
                directory_iterator it{p};
                while (it != directory_iterator{}) {
                    string entry_path = it->path().string();
                    result.push_back(entry_path);
                    *it++;
                }
                return result;
            }
        } else {
            cerr << p << " does not exist" << endl;
        }
    } catch (const filesystem_error &ex) {
        cerr << ex.what() << endl;
    }
    return {};
}

/**
 *
 * @param str
 * @param c
 * @return
 */
vector<string> split(const char *str, char c = ',') {
    vector<string> result;
    do {
        const char *begin = str;
        while (*str != c && *str)
            str++;
        result.push_back(string(begin, str));
    } while (0 != *str++);
    return result;
}

/**
 * Read a CSV file and convert it to a Graph object.
 * @param file_path
 * @return
 */
Graph *convert_file_to_graph_input(string file_path) {
    Graph *graph = new Graph();
    std::ifstream file(file_path);
    string current_line;
    if (file.good()) {
        vector<long double> x_axis;
        string titles;
        getline(file, titles);
        vector<string> splits = split(titles.c_str());
        graph->setX_axis_title(splits[0]);
        vector<string> y_axis_titles(splits.begin() + 1, splits.end());
        graph->setY_axes_titles(y_axis_titles);

        size_t no_of_y_graphs = splits.size() - 1;
        vector<vector<long double>> y_axes(no_of_y_graphs);
        while (getline(file, current_line)) {
            if (current_line.length() >= 1) {
                splits = split(current_line.c_str());
                x_axis.push_back(stod(splits[0]));
                vector<string> ys(splits.begin() + 1, splits.end());
                for (int i = 0; i < no_of_y_graphs; ++i) {
                    y_axes[i].push_back(stod(ys[i]));
                }
            }
        }
        graph->setX_axis(x_axis);
        graph->setY_axes(y_axes);
    } else {
        return NULL;
    }
    return graph;
}

/**
 * Main Function
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "No path to the data directory given!" << endl;
        cerr << "Usage: " << argv[0] << " <path to directory>" << endl;
        return 1;
    }
    vector<string> files = get_files(argv[1]);
    if (files.size() < 1)
        return 1;

    cout << "ID.\tFILENAME" << endl;
    for (int i = 0; i < files.size(); ++i)
        cout << i << ".\t" << files[i] << endl;
    cout << endl;

    int reference_file_index;
    int test_file_index;

    cout << "Enter the file ID for Reference: ";
    cin >> reference_file_index;

    cout << "Enter the file ID for Test: ";
    cin >> test_file_index;

    cout << endl;

    Graph *reference = convert_file_to_graph_input(files[reference_file_index]);
    Graph *test = convert_file_to_graph_input(files[test_file_index]);
    if (reference == NULL || test == NULL)
        return 1;

    cout << "Reference File: " << files[reference_file_index] << endl;
    cout << "Reference X_Axis Title: " << reference->getX_axis_title() << endl;
    cout << "Reference Y_Axis Title(s): " << reference->getY_axes_titles() << endl;
    cout << endl;

    cout << "Test File: " << files[test_file_index] << endl;
    cout << "Test X_Axis Title: " << test->getX_axis_title() << endl;
    cout << "Test Y_Axis Title(s): " << test->getY_axes_titles() << endl;
    cout << endl;

    int processing_index;
    cout << "Enter the title number for processing for Reference and Test data: ";
    cin >> processing_index;
    cout << endl;

    string x_axis_title = reference->getX_axis_title();
    string y_axis_title = reference->getY_axes_titles()[processing_index];

    cout << "Processing for [" << x_axis_title << " by " << y_axis_title << "]...." << endl;

    cout << endl;

    if (reference->is_valid_for_comparison(test)) {
        reference->setProcessing_index(processing_index);
        reference->process();

        test->setProcessing_index(processing_index);
        test->process();

        cout << "Peaks/Troughs for Reference data (in " << y_axis_title << "): ";
        cout << reference->getPeaks() << endl;

        cout << "Peaks/Troughs for Test data (in " << y_axis_title << "): ";
        cout << test->getPeaks() << endl;

        long double error = reference->relative_error(test->getY_axes()[processing_index]);
        long double corr = reference->correlation(test->getY_axes()[processing_index]);

        cout << "Normalized Squared Error: " << error * 100.0 << " %" << endl;
        cout << "Pearson Correlation: " << corr << endl;

    } else {
        cerr << "Check if both the inputs have same number of data points?" << endl;
    }

    return 0;
}
