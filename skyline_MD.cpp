#include <iostream>
#include <vector>
#include <limits>
#include <stack>

using namespace std;

int d, t;

struct Point {
    vector<double> coords; // d + t dimensions combined taken here
    int index;
};

struct FeatureBox {
    vector<pair<double, double>> bounds; 
};

struct QueryBox {
    vector<pair<double, double>> bounds;
};

vector<Point> points;

bool inRangeBox(const Point& p, const QueryBox& box) {
    for (int i = 0; i < d; i++) {
        if (p.coords[i] < box.bounds[i].first || p.coords[i] > box.bounds[i].second)
            return false;
    }
    return true;
}

bool inFeatureBox(const Point& p, const FeatureBox& box) {
    for (int i = 0; i < t; i++) {
        if (p.coords[d + i] < box.bounds[i].first || p.coords[d + i] > box.bounds[i].second)
            return false;
    }
    return true;
}

Point* findMinPoint(const QueryBox& qbox, const FeatureBox& fbox) {
    Point* best = nullptr;
    for (auto& p : points) {
        if (inRangeBox(p, qbox) && inFeatureBox(p, fbox)) {
            if (!best || p.coords[d + t - 1] < best->coords[d + t - 1]) {
                best = &p;
            }
        }
    }
    return best;
}

vector<FeatureBox> splitFeatureBox(const FeatureBox& fbox, const Point& p) {
    vector<FeatureBox> subBoxes;
    for (int i = 0; i < t - 1; i++) {
        if (p.coords[d + i] <= fbox.bounds[i].first) continue; 

        FeatureBox box = fbox;
        box.bounds[i].second = p.coords[d + i]; 
        box.bounds[t - 1].first = p.coords[d + t - 1] + 1e-6; 

        if (box.bounds[i].first < box.bounds[i].second && box.bounds[t - 1].first <= box.bounds[t - 1].second)
            subBoxes.push_back(box);
    }
    return subBoxes;
}

bool dominates(const Point& p, const Point& q) {
    bool strictlyBetter = false;
    for (int i = 0; i < t; i++) {
        if (p.coords[d + i] > q.coords[d + i]) return false;
        if (p.coords[d + i] < q.coords[d + i]) strictlyBetter = true;
    }
    return strictlyBetter;
}

void paperAlgorithm(const QueryBox& qbox) {
    FeatureBox fullBox;
    fullBox.bounds.resize(t);
    for (int i = 0; i < t; i++) {
        fullBox.bounds[i] = { -numeric_limits<double>::infinity(), numeric_limits<double>::infinity() };
    }

    stack<FeatureBox> Z;
    Z.push(fullBox);

    cout << "\nSkyline points (Paper Algorithm 1):\n";
    while (!Z.empty()) {
        FeatureBox R = Z.top(); Z.pop();
        Point* p = findMinPoint(qbox, R);

        if (p) {
            cout << "Point " << p->index << ": Range [ ";
            for (int i = 0; i < d; i++) cout << p->coords[i] << " ";
            cout << "] Features [ ";
            for (int i = 0; i < t; i++) cout << p->coords[d + i] << " ";
            cout << "]\n";

            auto subBoxes = splitFeatureBox(R, *p);
            for (auto it = subBoxes.rbegin(); it != subBoxes.rend(); ++it) {
                Z.push(*it);
            }
        }
    }
}

void bruteForceSkyline(const QueryBox& qbox) {
    vector<Point> candidates;
    for (auto& p : points) {
        if (inRangeBox(p, qbox)) candidates.push_back(p);
    }

    vector<Point> skyline;
    for (const auto& p : candidates) {
        bool dominated = false;
        for (const auto& q : candidates) {
            if (p.index != q.index && dominates(q, p)) {
                dominated = true;
                break;
            }
        }
        if (!dominated) skyline.push_back(p);
    }

    cout << "\nSkyline points (Brute-force method):\n";
    for (const auto& p : skyline) {
        cout << "Point " << p.index << ": Range [ ";
        for (int i = 0; i < d; i++) cout << p.coords[i] << " ";
        cout << "] Features [ ";
        for (int i = 0; i < t; i++) cout << p.coords[d + i] << " ";
        cout << "]\n";
    }
}

int main() {
    int n;
    cout << "Enter number of points: ";
    cin >> n;

    cout << "Enter range space dimension d and feature space dimension t: ";
    cin >> d >> t;

    points.resize(n);
    cout << "Enter " << d << " range + " << t << " feature coordinates per point:\n";
    for (int i = 0; i < n; i++) {
        points[i].coords.resize(d + t);
        for (int j = 0; j < d + t; j++) cin >> points[i].coords[j];
        points[i].index = i;
    }

    QueryBox qbox;
    qbox.bounds.resize(d);
    cout << "Enter query box (range space) [" << d << " intervals]:\n";
    for (int i = 0; i < d; i++) {
        cin >> qbox.bounds[i].first >> qbox.bounds[i].second;
    }

    paperAlgorithm(qbox);
    bruteForceSkyline(qbox);

    return 0;
}


/*
5
3 3
1 2 3  5 5 1
2 2 2  4 6 2
3 3 3  3 7 3
4 4 4  2 8 4
5 5 5  1 9 5
1 4
2 4
2 4
*/