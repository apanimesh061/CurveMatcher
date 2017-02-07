#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <deque>

#define MAX_ITER 10

using namespace std;

/**
 * Perform min-max normalisation on the input vector.
 * This is done so that the data points are scaled between 0 and 1 in
 * the first quadrant.
 * @param input     vector<long double>
 * @return output   vector<long double>
 */
vector<long double> normalize(vector<long double> input) {
    vector<long double> output;
    long double max = *max_element(input.begin(), input.end());
    long double min = *min_element(input.begin(), input.end());
    for (int i = 0; i < input.size(); ++i) {
        long double value = (input[i] - min) / (max - min);
        output.push_back(value);
    }
    return output;
}

/**
 * Returns average of the input vector
 * @param input     vector<long double>
 * @return t        long double
 */
long double average(vector<long double> input) {
    int n = (int) input.size();
    long double sum = 0.0;
    for (int i = 0; i < input.size(); ++i)
        sum += input[i];
    long double t = sum / n;
    return t;
}

/**
 * This is sort of a high-pass filter that filters values
 * lower than a threshold t
 * @param input     vector<long double>
 * @param t         long double
 * @return output   vector<long double>
 */
vector<long double> apply_threshold(vector<long double> input, long double t) {
    vector<long double> output;
    for (int i = 0; i < input.size(); ++i)
        output.push_back((input[i] >= t) ? input[i] : 0.0);
    return output;
}

/**
 * Returns a smoothed version of the input vector.
 * The smoothing is done using a 1D Gaussion Kernel [1, 4, 6, 4, 1]
 *
 * The input is convolved with the kernel successively MAX_ITER times.
 *
 * @param input     vector<long double>
 * @return output   vector<long double>
 */
vector<long double> apply_gaussian_filter(vector<long double> input) {
    int current_iter = 0;
    while (current_iter <= MAX_ITER) {
        vector<long double> output;
        output.reserve(input.size());
        int window = 5;

        output.push_back((6 * input[0] + 4 * input[1] + input[2]) / (long double) window);
        output.push_back((4 * input[0] + 6 * input[1] + 4 * input[2] + input[3]) / (long double) window);

        for (int i = 2; i < input.size() - 2; ++i)
            output.push_back(
                    (input[i - 2] + 4 * input[i - 1] + 6 * input[i] + 4 * input[i + 1] + input[i + 2]) /
                    (long double) window);

        output.push_back((4 * input[input.size() - 1] + 6 * input[input.size() - 2] + 4 * input[input.size() - 3] +
                          input[input.size() - 4]) / (long double) window);
        output.push_back(
                (6 * input[input.size() - 1] + 4 * input[input.size() - 2] + input[input.size() - 3]) /
                (long double) window);
        input = output;
        output.clear();
        current_iter++;
    }
    return input;
}

/**
 * Gets the maximum/minimum sliding window output of an input vector.
 * The input is first padded on the end and start according to size of
 * the input window size. The window size will always be odd.
 *
 *
 * @param input                                         vector<long double>
 * @param w                                             int
 * @param is_max                                        bool
 * @return eroded or dilated version of input vector    vector<long double>
 */
vector<long double> sliding_window(vector<long double> input, int w, bool is_max) {
    vector<long double> padding((size_t) ((w - 1) / 2));
    long double padding_value = numeric_limits<long double>::infinity();
    fill(padding.begin(), padding.end(), (is_max) ? -padding_value : padding_value);

    vector<long double> padded_input;
    padded_input.insert(padded_input.end(), padding.begin(), padding.end());
    padded_input.insert(padded_input.end(), input.begin(), input.end());
    padded_input.insert(padded_input.end(), padding.begin(), padding.end());

    size_t n = padded_input.size();
    vector<long double> output(n - w + 1);
    deque<long double> dq;
    for (int i = 0; i < w; ++i) {
        if (is_max) {
            while (!dq.empty() && padded_input[i] >= padded_input[dq.back()])
                dq.pop_back();
        } else {
            while (!dq.empty() && padded_input[i] <= padded_input[dq.back()])
                dq.pop_back();
        }
        dq.push_back(i);
    }
    for (int i = w; i < n; ++i) {
        output[i - w] = padded_input[dq.front()];
        if (is_max) {
            while (!dq.empty() && padded_input[i] >= padded_input[dq.back()])
                dq.pop_back();
        } else {
            while (!dq.empty() && padded_input[i] <= padded_input[dq.back()])
                dq.pop_back();
        }
        while (!dq.empty() && dq.front() <= (i - w))
            dq.pop_front();
        dq.push_back(i);
    }
    output[n - w] = padded_input[dq.front()];
    return output;
}

/**
 * This function implements the Erosion morphological operation
 * @param input
 * @param w
 * @return
 */
vector<long double> erosion(vector<long double> input, int w) {
    return sliding_window(input, w, false);
}

/**
 * This function implements the Dilation morphological operation
 * @param input
 * @param w
 * @return
 */
vector<long double> dilation(vector<long double> input, int w) {
    return sliding_window(input, w, true);
}

/**
 * This is done using a filter called Top Hat Filter that ulitizes the
 * max or min sliding window method.
 *
 * structuring_element_size is window size which is used as input window
 * size in calculating the min/max sliding window output
 *
 * @param input
 * @param structuring_element_size
 * @return
 */
vector<long double> apply_white_tophat_filter(vector<long double> input, int structuring_element_size) {
    vector<long double> result;
    vector<long double> temp = dilation(erosion(input, structuring_element_size), structuring_element_size);
    for (int i = 0; i < input.size(); ++i)
        result.push_back(input[i] - temp[i]);
    return result;
}

/**
 * This is done using a filter called Top Hat Filter that ulitizes the
 * max or min sliding window method.
 *
 * structuring_element_size is window size which is used as input window
 * size in calculating the min/max sliding window output
 *
 * @param input
 * @param structuring_element_size
 * @return
 */
vector<long double> apply_black_tophat_filter(vector<long double> input, int structuring_element_size) {
    vector<long double> result;
    vector<long double> temp = erosion(dilation(input, structuring_element_size), structuring_element_size);
    for (int i = 0; i < input.size(); ++i)
        result.push_back(temp[i] - input[i]);
    return result;
}

/**
 * Performs a local search on a peak or trough point in order to find any near by maxima/minima.
 *
 * When a Gaussian Filter is applied on vector, it results in smoothing which is technically data
 * loss. So, when we try to find features on the smoothed version we might we having some error
 * margin.
 *
 * Here, start is the index of a feature on the processed version of the input.
 * We first check if the index is actaully a peak or trough on the original normalized version.
 * If not, we check all points, within a window size of 10 on left as well as right side of start,
 * whether or not they are peaks or troughs.
 *
 * We return the nearest index of all the 20 indices checked.
 *
 * @param input                 vector<long double>
 * @param start                 int
 * @return final feature index  int
 */
int local_search(vector<long double> input, int start) {
    if ((input[start - 1] > input[start] && input[start + 1] > input[start]) ||
        (input[start - 1] < input[start] && input[start + 1] < input[start]))
        return start;

    int local_search_window = 10;
    int current = start;
    int nearest = start;
    int diff = numeric_limits<int>::max();

    int left_end = start - local_search_window;
    while (current > left_end && current > 0 && current < input.size() - 1) {
        if ((input[current - 1] > input[current] && input[current + 1] > input[current]) ||
            (input[current - 1] < input[current] && input[current + 1] < input[current])) {
            if ((abs(current - start) < diff)) {
                diff = (int)abs(current - start);
                nearest = current;
            }
        }
        current--;
    }

    current = start;
    int right_end = start + local_search_window;
    while (current < right_end && current > 0 && current < input.size() - 1) {
        if ((input[current - 1] > input[current] && input[current + 1] > input[current]) ||
            (input[current - 1] < input[current] && input[current + 1] < input[current])) {
            if ((abs(current - start) < diff)) {
                diff = (int)abs(current - start);
                nearest = current;
            }
        }
        current++;
    }

    return nearest;
}

/**
 * Calculates the First Order Derivative of a 1D Input
 * For a 1D input the derivative of successive differences of
 * elements in input.
 *
 * @param input                             vector<long double>
 * @return first order derivative of input  vector<long double>
 */
vector<long double> first_order_derivative(vector<long double> input) {
    vector<long double> result;
    result.push_back(0.0);
    for (int i = 1; i < input.size(); ++i) {
        result.push_back(input[i] - input[i - 1]);
    }
    return result;
}

/**
 * Checks for zero crossings on the First Order Derivative output from
 * first_order_derivative(vector<long double> input)
 *
 * Zero crossings tell us that there was a minima or maxima on that index.
 *
 * @param input                 vector<long double>
 * @return vector of indices    vector<int>
 */
vector<int> get_peak_indices(vector<long double> input) {
    vector<long double> fod = first_order_derivative(input);
    vector<int> peaks;
    for (int i = 1; i < fod.size(); ++i) {
        if (fod[i] * fod[i - 1] < 0)
            peaks.push_back(i - 1);
    }
    return peaks;
}
