#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include "filter.h"

using namespace std;

/**
 * This class represents an input for processing
 */
class Graph {
private:
    vector<long double> x_axis;
    vector<vector<long double>> y_axes;
    string x_axis_title;
    vector<string> y_axes_titles;
    vector<long double> peaks;
    int processing_index = -1;

public:
    Graph();

    virtual ~Graph();

    const string &getX_axis_title() const;

    void setX_axis_title(const string &x_axis_title);

    const vector<string> &getY_axes_titles() const;

    void setY_axes_titles(const vector<string> &y_axes_titles);

    const vector<long double> &getX_axis() const;

    void setX_axis(const vector<long double> &x_axis);

    const vector<vector<long double>> &getY_axes() const;

    void setY_axes(const vector<vector<long double>> &y_axes);

    const vector<long double> &getPeaks() const;

    void setProcessing_index(int processing_index);

    void process();

    bool is_valid_for_comparison(Graph *input);

    long double relative_error(vector<long double> other);

    long double correlation(vector<long double> other);
};

Graph::Graph() {}

const vector<long double> &Graph::getX_axis() const {
    return x_axis;
}

void Graph::setX_axis(const vector<long double> &x_axis) {
    Graph::x_axis = x_axis;
}

const vector<vector<long double>> &Graph::getY_axes() const {
    return y_axes;
}

void Graph::setY_axes(const vector<vector<long double>> &y_axes) {
    Graph::y_axes = y_axes;
}

Graph::~Graph() {
    Graph::x_axis.clear();
    Graph::y_axes.clear();
}

/**
 * Returns true iff x-axis if same for both graphs
 * @param test  Graph
 * @return bool
 */
bool Graph::is_valid_for_comparison(Graph *test) {
    return (Graph::x_axis == test->getX_axis());
}

const string &Graph::getX_axis_title() const {
    return x_axis_title;
}

void Graph::setX_axis_title(const string &x_axis_title) {
    Graph::x_axis_title = x_axis_title;
}

const vector<string> &Graph::getY_axes_titles() const {
    return y_axes_titles;
}

void Graph::setY_axes_titles(const vector<string> &y_axes_titles) {
    Graph::y_axes_titles = y_axes_titles;
}

/**
 * Return all peaks/troughs of the graph
 * @return
 */
const vector<long double> &Graph::getPeaks() const {
    return peaks;
}

/**
 * This is core of the whole processing.
 *
 * This functions sets all the required features of the input.
 */
void Graph::process() {
    vector<long double> y = getY_axes()[Graph::processing_index];
    vector<long double> y_norm = normalize(y);

    vector<long double> y_smoothed = apply_gaussian_filter(y_norm);

    /*
     * This part ensures that the structuring element for top hat filter
     * if odd in size.
     */
    int possible_window_size = (int) (0.1 * y_norm.size());
    if (possible_window_size % 2 == 0)
        possible_window_size++;
    possible_window_size = max(possible_window_size, 3);

    // Apply all top hat filter operations
    vector<long double> y_bth = apply_black_tophat_filter(y_smoothed, possible_window_size);
    vector<long double> y_wth = apply_white_tophat_filter(y_smoothed, possible_window_size);

    // Apply threshold to the filtered versions of input.
    // This will act as High Pass Filter for the input.
    vector<long double> y_bthreshed = apply_threshold(y_bth, average(y_bth));
    vector<long double> y_wthreshed = apply_threshold(y_wth, average(y_wth));

    vector<int> peak_indices;
    for (int peak : get_peak_indices(y_bthreshed))
        peak_indices.push_back(local_search(y_norm, peak));
    for (int peak : get_peak_indices(y_wthreshed))
        peak_indices.push_back(local_search(y_norm, peak));

    for (auto index : peak_indices)
        Graph::peaks.push_back(y[index]);
}

/**
 * @param processing_index
 */
void Graph::setProcessing_index(int processing_index) {
    Graph::processing_index = processing_index;
}

/**
 * This function calculates the relative sqaured error of current graph with
 * another graph.
 * @param other         vector<long double>
 * @return long double (0.0 to 1.0)
 */
long double Graph::relative_error(vector<long double> other) {
    vector<long double> ref_y = normalize(Graph::y_axes[Graph::processing_index]);
    vector<long double> other_norm = normalize(other);

    long double error = 0.0;
    long double factor = 0.0;
    for (int i = 0; i < ref_y.size(); ++i) {
        long double diff = ref_y[i] - other_norm[i];
        error += pow(fabs(diff), 2);
        factor += pow(fabs(ref_y[i]), 2);
    }
    return error / factor;
}

/**
 * Calculates the Pearson Correlation Coefficient of current graph with another graph.
 *
 * @param other         vector<long double>
 * @return long double (-1.0 to 1.0)
 */
long double Graph::correlation(vector<long double> other) {
    vector<long double> ref_y = normalize(Graph::y_axes[Graph::processing_index]);
    vector<long double> other_norm = normalize(other);

    int n = (int) ref_y.size();

    long double ref_y_avg = average(ref_y);
    long double other_avg = average(other_norm);

    long double numerator = 0.0;
    long double var_ref_y = 0.0;
    long double var_other = 0.0;
    for (int i = 0; i < n; ++i) {
        numerator += ((ref_y[i] - ref_y_avg) * (other_norm[i] - other_avg));
        var_ref_y += pow((ref_y[i] - ref_y_avg), 2);
        var_other += pow((other_norm[i] - other_avg), 2);
    }

    return numerator / (sqrt(var_ref_y) * sqrt(var_other));
}
