#include <iostream>
#include <vector>
#include <limits>
#include <stack>
#include <algorithm>
#include <unordered_set>
using namespace std;

int d, t;

struct Point {
    vector<double> coords; // d + t dimensions combined
    int index;
};

struct QueryBox {
    vector<pair<double, double>> bounds;
};

struct FeatureBox {
    vector<pair<double, double>> bounds;
};

struct RangeTreeNode {
    int level;
    int split_val;
    RangeTreeNode *left;
    RangeTreeNode *right;
    RangeTreeNode *next_tree;
    Point *min_point;
    Point *curr_point;

    RangeTreeNode(int l, int v, RangeTreeNode *le, RangeTreeNode *ri, RangeTreeNode *ne, Point *cp, Point *mp)
        : level(l), split_val(v), left(le), right(ri), next_tree(ne), curr_point(cp), min_point(mp) {}
};

vector<Point> points;

bool inRangeBox(const Point &p, const QueryBox &box) {
    for (int i = 0; i < d; ++i) {
        if (p.coords[i] < box.bounds[i].first || p.coords[i] > box.bounds[i].second)
            return false;
    }
    return true;
}

bool inFeatureBox(const Point &p, const FeatureBox &box) {
    for (int i = 0; i < t; ++i) {
        if (p.coords[d + i] < box.bounds[i].first || p.coords[d + i] > box.bounds[i].second)
            return false;
    }
    return true;
}

bool min_all(Point *p1, Point *p2) {
    for (int i = t - 1; i >= 0; i--) {
        if (p1->coords[d + i] < p2->coords[d + i]) return true;
        if (p1->coords[d + i] > p2->coords[d + i]) return false;
    }
    return false;
}

RangeTreeNode* build_tree(vector<Point *> &pts, int level, int total_levels) {
    if (pts.empty()) return nullptr;

    int dim = (level < d) ? level : (level - d);

    sort(pts.begin(), pts.end(), [&](Point *a, Point *b) {
        return a->coords[dim + (level < d ? 0 : d)] < b->coords[dim + (level < d ? 0 : d)];
    });

    int mid = pts.size() / 2;
    Point *curr = pts[mid];
    int split_val = curr->coords[dim + (level < d ? 0 : d)];

    vector<Point *> left_pts(pts.begin(), pts.begin() + mid);
    vector<Point *> right_pts(pts.begin() + mid + 1, pts.end());

    RangeTreeNode *left = build_tree(left_pts, level, total_levels);
    RangeTreeNode *right = build_tree(right_pts, level, total_levels);
    RangeTreeNode *next = nullptr;
    if (level < total_levels - 1)
        next = build_tree(pts, level + 1, total_levels);

    Point *min_p = curr;
    if (left && min_all(left->min_point, min_p)) min_p = left->min_point;
    if (right && min_all(right->min_point, min_p)) min_p = right->min_point;
    if (next && min_all(next->min_point, min_p)) min_p = next->min_point;

    return new RangeTreeNode(level, split_val, left, right, next, curr, min_p);
}

void range_query(RangeTreeNode *node, QueryBox &qbox, FeatureBox &fbox, int level, int total_levels, vector<RangeTreeNode *> &canonicals) {
    if (!node) return;

    int dim = (level < d) ? level : (level - d);
    double l = (level < d) ? qbox.bounds[dim].first : fbox.bounds[dim].first;
    double r = (level < d) ? qbox.bounds[dim].second : fbox.bounds[dim].second;
    double val = node->curr_point->coords[dim + (level < d ? 0 : d)];

    if (val < l) {
        range_query(node->right, qbox, fbox, level, total_levels, canonicals);
    } else if (val > r) {
        range_query(node->left, qbox, fbox, level, total_levels, canonicals);
    } else {
        if (level == total_levels - 1) {
            canonicals.push_back(node);
        } else {
            range_query(node->next_tree, qbox, fbox, level + 1, total_levels, canonicals);
            range_query(node->left, qbox, fbox, level, total_levels, canonicals);
            range_query(node->right, qbox, fbox, level, total_levels, canonicals);
        }
    }
}

vector<FeatureBox> splitFeatureBox(const FeatureBox &fbox, const Point &p) {
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

void paperAlgorithmWithRangeTree(RangeTreeNode *root, const QueryBox &qbox) {
    FeatureBox fullBox;
    fullBox.bounds.resize(t);
    for (int i = 0; i < t; i++) {
        fullBox.bounds[i] = { -numeric_limits<double>::infinity(), numeric_limits<double>::infinity() };
    }

    stack<FeatureBox> Z;
    Z.push(fullBox);

    unordered_set<int>reported;
    cout << "\nSkyline points (Range Tree Based):\n";
    while (!Z.empty()) {
        FeatureBox R = Z.top(); Z.pop();
        vector<RangeTreeNode *> canonicals;
        range_query(root, const_cast<QueryBox &>(qbox), R, 0, d + t, canonicals);

        Point *best = nullptr;
        for (auto *node : canonicals) {
            if (!inRangeBox(*node->curr_point, qbox) || !inFeatureBox(*node->curr_point, R)) continue;
            if (!best || min_all(node->curr_point, best)) {
                best = node->curr_point;
            }
        }

        if (best) {
            cout << "Point " << best->index << ": Range [ ";
            for (int i = 0; i < d; i++) cout << best->coords[i] << " ";
            cout << "] Features [ ";
            for (int i = 0; i < t; i++) cout << best->coords[d + i] << " ";
            cout << "]\n";

            auto subBoxes = splitFeatureBox(R, *best);
            for (auto it = subBoxes.rbegin(); it != subBoxes.rend(); ++it) {
                Z.push(*it);
            }
            reported.insert(best->index);
        }
    }
}

bool dominates(const Point& p, const Point& q) {
    bool strictlyBetter = false;
    for (int i = t-1; i >= 0; i--) {
        if (p.coords[d + i] > q.coords[d + i]) return false;
        if (p.coords[d + i] < q.coords[d + i]) strictlyBetter = true;
    }
    return strictlyBetter;
}

void bruteForceSkyline(const QueryBox& qbox) {
    vector<Point> candidates;
    for (auto& p : points) {
        if (inRangeBox(p, qbox)) {candidates.push_back(p);cout<<p.index<<" ";}
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

    vector<Point *> point_ptrs(n);
    for (int i = 0; i < n; ++i) point_ptrs[i] = &points[i];

    RangeTreeNode *root = build_tree(point_ptrs, 0, d + t);

    QueryBox qbox;
    qbox.bounds.resize(d);
    cout << "Enter query box (range space) [" << d << " intervals]:\n";
    for (int i = 0; i < d; i++) {
        cin >> qbox.bounds[i].first >> qbox.bounds[i].second;
    }

    paperAlgorithmWithRangeTree(root, qbox);
    bruteForceSkyline(qbox);
    return 0;
}

/*
10 
4 7
12 35 15 10   5 12 35 8 15 10 10
18 45 14 9    6 18 45 9 14 9 11
15 38 17 11   7 15 38 8 17 11 9
13 50 15 13   4 13 50 7 15 13 10
19 59 13 8    8 19 59 8 13 8 12
10 40 16 14   3 10 40 7 16 14 10
11 37 15 12   5 11 37 8 15 12 10
17 60 15 9    2 17 60 6 15 9 13
22 42 14 11   6 22 42 9 14 11 9
9  33 17 10   1 9  33 5 17 10 10
10 20
30 60
5 15
5 20
*/