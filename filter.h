#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <deque>

#define MAX_ITER 10

using namespace std;

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

long double average(vector<long double> input) {
    int n = (int) input.size();
    long double sum = 0.0;
    for (int i = 0; i < input.size(); ++i)
        sum += input[i];
    long double t = sum / n;
    return t;
}

vector<long double> apply_threshold(vector<long double> input, long double t) {
    vector<long double> output;
    for (int i = 0; i < input.size(); ++i)
        output.push_back((input[i] >= t) ? input[i] : 0.0);
    return output;
}

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

vector<long double> erosion(vector<long double> input, int w) {
    return sliding_window(input, w, false);
}

vector<long double> dilation(vector<long double> input, int w) {
    return sliding_window(input, w, true);
}

vector<long double> apply_white_tophat_filter(vector<long double> input, int structuring_element_size) {
    vector<long double> result;
    vector<long double> temp = dilation(erosion(input, structuring_element_size), structuring_element_size);
    for (int i = 0; i < input.size(); ++i)
        result.push_back(input[i] - temp[i]);
    return result;
}

vector<long double> apply_black_tophat_filter(vector<long double> input, int structuring_element_size) {
    vector<long double> result;
    vector<long double> temp = erosion(dilation(input, structuring_element_size), structuring_element_size);
    for (int i = 0; i < input.size(); ++i)
        result.push_back(temp[i] - input[i]);
    return result;
}

vector<long double> first_order_derivative(vector<long double> input) {
    vector<long double> result;
    result.push_back(0.0);
    for (int i = 1; i < input.size(); ++i) {
        result.push_back(input[i] - input[i - 1]);
    }
    return result;
}

int local_search(vector<long double> input, int start) {
    if ((input[start - 1] > input[start] && input[start + 1] > input[start]) ||
        (input[start - 1] < input[start] && input[start + 1] < input[start]))
        return start;

    int local_search_window = 10;
    int current = start;
    int nearest = start;
    int diff = numeric_limits<int>::max();

    int left_end = start - local_search_window;
    while (current > left_end) {
        if ((input[current - 1] > input[current] && input[current + 1] > input[current]) ||
            (input[current - 1] < input[current] && input[current + 1] < input[current])) {
            if ((abs(current - start) < diff)) {
                diff = (int) abs(current - start);
                nearest = current;
            }
        }
        current--;
    }

    current = start;
    int right_end = start + local_search_window;
    while (current < right_end) {
        if ((input[current - 1] > input[current] && input[current + 1] > input[current]) ||
            (input[current - 1] < input[current] && input[current + 1] < input[current])) {
            if ((abs(current - start) < diff)) {
                diff = (int) abs(current - start);
                nearest = current;
            }
        }
        current++;
    }

    return nearest;
}

vector<int> get_peak_indices(vector<long double> input) {
    vector<long double> fod = first_order_derivative(input);
    vector<int> peaks;
    for (int i = 1; i < fod.size(); ++i) {
        if (fod[i] * fod[i - 1] < 0)
            peaks.push_back(i - 1);
    }
    return peaks;
}
