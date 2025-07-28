#include <iostream>
#include <vector>
#include <limits>
#include <stack>
#include <algorithm>
#include <unordered_set>
#include <fstream>
#include <queue>
#include <bits/stdc++.h>
using namespace std;

int d, t;
double eps = 1e-9;
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
    double split_val;
    int subtree_size;
    RangeTreeNode *left;
    RangeTreeNode *right;
    RangeTreeNode *next_tree;
    Point *min_point;
    Point *curr_point;

    RangeTreeNode(int l, double v, int subcnt, RangeTreeNode *le, RangeTreeNode *ri, RangeTreeNode *ne, Point *cp, Point *mp)
        : level(l), split_val(v), subtree_size(subcnt), left(le), right(ri), next_tree(ne), curr_point(cp), min_point(mp) {}
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

RangeTreeNode* build_tree(vector<Point *> &pts, int level, int last_dim, int total_levels) {
    if (pts.empty()) return nullptr;

    int dim = level >= last_dim ? level+1:level;
    if(level==total_levels-1)
        dim = last_dim;
    sort(pts.begin(), pts.end(), [&](Point *a, Point *b) {
        return a->coords[dim] < b->coords[dim];
    });
    int mid = pts.size() / 2;
    Point *curr = pts[mid];
    double split_val = curr->coords[dim];

    vector<Point *> left_pts(pts.begin(), pts.begin() + mid);
    vector<Point *> right_pts(pts.begin() + mid + 1, pts.end());

    RangeTreeNode *left = build_tree(left_pts, level, last_dim, total_levels);
    RangeTreeNode *right = build_tree(right_pts, level, last_dim, total_levels);
    RangeTreeNode *next = nullptr;
    if (level < total_levels - 1)
        next = build_tree(pts, level + 1, last_dim, total_levels);

    Point *min_p = curr;
    if (left && min_all(left->min_point, min_p)) min_p = left->min_point;
    if (right && min_all(right->min_point, min_p)) min_p = right->min_point;
    if (next && min_all(next->min_point, min_p)) min_p = next->min_point;

    return new RangeTreeNode(level, split_val, (int)(pts.size()), left, right, next, curr, min_p);
}

void range_query(RangeTreeNode *node, QueryBox &qbox, FeatureBox &fbox, int level, int total_levels, vector<RangeTreeNode *> &canonicals) {
    if (!node) return;

    int dim = (level < d) ? level : (level - d);
    double l = (level < d) ? qbox.bounds[dim].first : fbox.bounds[dim].first;
    double r = (level < d) ? qbox.bounds[dim].second : fbox.bounds[dim].second;
    double val = node->curr_point->coords[level];

    if (val < l) {
        range_query(node->right, qbox, fbox, level, total_levels, canonicals);
    } else if (val > r) {
        range_query(node->left, qbox, fbox, level, total_levels, canonicals);
    } else {
        if (level == total_levels - 1) {
            canonicals.push_back(node);
            range_query(node->left, qbox, fbox, level, total_levels, canonicals);
            range_query(node->right, qbox, fbox, level, total_levels, canonicals);
        } else {
            range_query(node->next_tree, qbox, fbox, level + 1, total_levels, canonicals);
            range_query(node->left, qbox, fbox, level, total_levels, canonicals);
            range_query(node->right, qbox, fbox, level, total_levels, canonicals);
        }
    }
}
bool inRangeQ(Point* p, QueryBox & Q, int D, int T) {
    for (int i = 0; i < D; ++i) {
        if (p->coords[i] < Q.bounds[i].first || p->coords[i] > Q.bounds[i].second)
            return false;
    }
    return true;
}

bool inRangeF(Point* p, FeatureBox & F, int D, int T) {
    for (int i = D; i < D+T; ++i) {
        if (p->coords[i] < F.bounds[i].first || p->coords[i] > F.bounds[i].second)
            return false;
    }
    return true;
}

RangeTreeNode* findSplitNode(RangeTreeNode* root, double low, double high) {
    while (root && (root->left || root->right) &&
           (high < root->split_val || low > root->split_val)) {
        if (high < root->split_val)
            root = root->left;
        else
            root = root->right;
    }
    return root;
}

void dDRangeQuery(RangeTreeNode* root,  QueryBox& Q, FeatureBox& F, int level, int total_levels, int d, int t, vector<Point*>& out, int &size) {
    if (!root || level >= total_levels) return;

    int dim = (level < d) ? level : (level - d);
    double low = (level < d) ? Q.bounds[dim].first : F.bounds[dim].first;
    double high = (level < d) ? Q.bounds[dim].second : F.bounds[dim].second;

    RangeTreeNode* split = findSplitNode(root, low, high);
    if (!split) return;
    out.push_back(split->curr_point);

    if(inRangeQ(split->curr_point, Q, d, t)){
        size++;
    }
    // Traverse left path;
    RangeTreeNode* v = split->left;
    while (v) {
        if (low <= v->split_val) {
            if(inRangeQ(v->curr_point, Q, d, t)){
                out.push_back(v->curr_point);
                size++;
            }
            if(level == total_levels-1){
                if(v->right) size+= v->right->subtree_size;
            }
            if(v->right && v->right->next_tree){
                dDRangeQuery(v->right->next_tree, Q, F, level + 1, total_levels, d, t, out, size);
            }
            v = v->left;
        } else {
            v = v->right;
        }
    }
    
    v = split->right;
    while (v) {
        if (high >= v->split_val) {
            if(inRangeQ(v->curr_point, Q, d, t)){
                out.push_back(v->curr_point);
                size++;
            }
            if(level == total_levels-1){
                if(v->left) size+= v->left->subtree_size;
            }
            if(v->left && v->left->next_tree){
                dDRangeQuery(v->left->next_tree, Q, F, level + 1, total_levels, d, t, out, size);
            }
            v = v->right;
        } else {
            v = v->left;
        }
    }
}


void side_expand_2(RangeTreeNode* root,  QueryBox& Q, FeatureBox& F, int level, int total_levels, int d, int t, vector<Point*>& out, int &size) {
    if (!root || level >= total_levels) return;

    int dim = (level < d) ? level : (level - d);
    double low = (level < d) ? Q.bounds[dim].first : F.bounds[dim].first;
    double high = (level < d) ? Q.bounds[dim].second : F.bounds[dim].second;

    RangeTreeNode* split = findSplitNode(root, low, high);
    if (!split) return;
    out.push_back(split->curr_point);

    if(inRangeQ(split->curr_point, Q, d, t)){
        size++;
    }
    // Traverse left path;
    RangeTreeNode* v = split->left;
    while (v) {
        if (low <= v->split_val) {
            if(inRangeQ(v->curr_point, Q, d, t)){
                out.push_back(v->curr_point);
                size++;
            }
            if(level == total_levels-1){
                if(v->right) size+= v->right->subtree_size;
                if(v->left && !inRangeQ(v->left->curr_point, Q, d, t)){
                    out.push_back(v->left->curr_point);
                    size++;
                }
            }
            if(v->right && v->right->next_tree){
                side_expand_2(v->right->next_tree, Q, F, level + 1, total_levels, d, t, out, size);
            }
            v = v->left;
        } else {
            v = v->right;
        }
    }
    
    v = split->right;
    while (v) {
        if (high >= v->split_val) {
            if(inRangeQ(v->curr_point, Q, d, t)){
                out.push_back(v->curr_point);
                size++;
            }
            if(level == total_levels-1){
                if(v->left) size+= v->left->subtree_size;
                if(v->right && !inRangeQ(v->right->curr_point, Q, d, t)){
                    out.push_back(v->right->curr_point);
                    size++;
                }
            }
            if(v->left && v->left->next_tree){
                side_expand_2(v->left->next_tree, Q, F, level + 1, total_levels, d, t, out, size);
            }
            v = v->right;
        } else {
            v = v->left;
        }
    }
}



vector<FeatureBox> splitFeatureBox(const FeatureBox &fbox, const Point &p) {
    vector<FeatureBox> subBoxes;
    for(int i=0; i<t-1; i++){ 
        FeatureBox box=fbox;    
        // box.bounds.resize(t);
        // for (int i = 0; i < t; i++) {
        //     box.bounds[i] = { -numeric_limits<double>::infinity(), numeric_limits<double>::infinity() };
        // }
        for(int j=0; j<i; j++){    
            box.bounds[j].first = max(p.coords[d + j],box.bounds[j].first);
        }
        box.bounds[i].second = min(p.coords[d + i] - 1e-9,box.bounds[i].second);
        box.bounds[t-1].first=max(p.coords[d + t - 1], box.bounds[t-1].first); 
        bool valid = true;
        for (int j = 0; j < t; j++) {
            if (box.bounds[j].first > box.bounds[j].second) {
                valid = false;
                break;
            }
        }
        
        if (valid) {
            subBoxes.push_back(box);
        }
    }
    return subBoxes;
}

vector<Point*> paperAlgorithmWithRangeTree(RangeTreeNode *root, const QueryBox &qbox) {
    FeatureBox fullBox;
    fullBox.bounds.resize(t);
    for (int i = 0; i < t; i++) {
        fullBox.bounds[i] = { -numeric_limits<double>::infinity(), numeric_limits<double>::infinity() };
    }

    stack<FeatureBox> Z;
    Z.push(fullBox);

    unordered_set<int>reported;
    vector<Point*> skyline;
    cout << "\nSkyline points (Range Tree Based):\n";
    while (!Z.empty()) {
        FeatureBox R = Z.top(); Z.pop();
        vector<Point *> canonicals;
        
        int size = 0;
        dDRangeQuery(root, const_cast<QueryBox &>(qbox), R, 0, d+t, d, t, canonicals,size);
        Point *best = nullptr;
        for (auto pt : canonicals) {
            if (!inRangeBox(*pt, qbox) || !inFeatureBox(*pt, R)) continue;
            if (!best || min_all(pt, best)) {
                best = pt;
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
            if(reported.count(best->index)==0){
                skyline.push_back(best);
                reported.insert(best->index);
            }
        }
    }
    return skyline;
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
        if (inRangeBox(p, qbox)) {candidates.push_back(p);}
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


vector<QueryBox> generateCornerSkylineRanges(
    const QueryBox& qbox,
    const vector<double>& best_lhs,
    const vector<double>& best_rhs) {

    int d = qbox.bounds.size();
    vector<QueryBox> corner_ranges;

    int total_masks = 1 << d;

    for (int mask = 0; mask < total_masks; ++mask) {
        QueryBox newbox = qbox;
        for (int i = 0; i < d; ++i) {
            if ((mask >> i) & 1) {
                if (best_rhs[i] < 1e17)
                    newbox.bounds[i].second = best_rhs[i];
            } else {
                if (best_lhs[i] > -1e17)
                    newbox.bounds[i].first = best_lhs[i];
            }
        }

        corner_ranges.push_back(newbox);
    }

    return corner_ranges;
}

void bfs(vector<RangeTreeNode *> &trees, QueryBox &qbox, RangeTreeNode* root){
    FeatureBox fullBox;
    fullBox.bounds.resize(t);
    for (int i = 0; i < t; i++) {
        fullBox.bounds[i] = { -numeric_limits<double>::infinity(), numeric_limits<double>::infinity() };
    }

    priority_queue<pair<double, int>>pq;
    vector<QueryBox>boxes;
    boxes.push_back(qbox);
    pq.push({1e9,0});
    int base_size=0;
    vector<Point *> base_canonicals;
    dDRangeQuery(root, const_cast<QueryBox &>(qbox), fullBox, 0, d, d, t, base_canonicals,base_size);

    while(!pq.empty()){
        double sim = pq.top().first; int ind = pq.top().second; QueryBox box = boxes[ind];
        pq.pop();
        vector<double> lhs_expands, rhs_expands;
        for (int dim = 0; dim < d; ++dim) {
            vector<Point*> out;
            int size=0;
            side_expand_2(trees[dim], box, fullBox, 0, d, d, t,  out, size);
            double lhs_bound=LONG_MIN,rhs_bound = LONG_MAX;
            double lhs_shrink_pt=box.bounds[dim].second, rhs_shrink_pt=box.bounds[dim].first;
            for(auto p: out){
                if (p->coords[dim] < box.bounds[dim].first)
                    lhs_bound = max(lhs_bound, p->coords[dim]);
                
                if (p->coords[dim] > box.bounds[dim].second)
                    rhs_bound = min(rhs_bound, p->coords[dim]);
                
                if(p->coords[dim]> box.bounds[dim].first && p->coords[dim]< box.bounds[dim].second){
                    lhs_shrink_pt = min(lhs_shrink_pt, p->coords[dim]);
                    rhs_shrink_pt = max(rhs_shrink_pt, p->coords[dim]);
                }
                lhs_bound = min(lhs_bound, p->coords[dim]);
                rhs_bound = max(rhs_bound, p->coords[dim]);
            }
            QueryBox lhs_expanded_box = box, rhs_expanded_box = box, lhs_shrink_box = box, rhs_shrink_box = box;
            vector<Point*> canonicals;
            int intersection_size, union_size;
            lhs_expands.push_back(lhs_bound);
            rhs_expands.push_back(rhs_bound);

            lhs_expanded_box.bounds[dim].first = lhs_bound;
            union_size=0;
            intersection_size = base_size;
            dDRangeQuery(trees[dim], const_cast<QueryBox &>(lhs_expanded_box), fullBox, 0, d, d, t, canonicals,union_size);
            double lhs_sim = (double)(intersection_size)/(double)(union_size);
            boxes.push_back(lhs_expanded_box);
            pq.push({lhs_sim,++ind});

            rhs_expanded_box.bounds[dim].second = rhs_bound;
            union_size=0;
            intersection_size=base_size;
            dDRangeQuery(trees[dim], const_cast<QueryBox &>(rhs_expanded_box), fullBox, 0, d, d, t, canonicals,union_size);
            double rhs_sim = (double)(intersection_size)/(double)(union_size);
            boxes.push_back(rhs_expanded_box);
            pq.push({rhs_sim, ++ind});


            lhs_shrink_box.bounds[dim].first = lhs_shrink_pt;
            union_size=base_size;
            intersection_size=0;
            dDRangeQuery(trees[dim], const_cast<QueryBox &>(lhs_shrink_box), fullBox, 0, d, d, t, canonicals,intersection_size);
            lhs_sim = (double)(intersection_size)/(double)(union_size);
            boxes.push_back(lhs_shrink_box);
            pq.push({lhs_sim, ++ind});

            rhs_shrink_box.bounds[dim].second = rhs_shrink_pt;
            union_size=base_size;
            intersection_size=0;
            dDRangeQuery(trees[dim], const_cast<QueryBox &>(rhs_shrink_box), fullBox, 0, d, d, t, canonicals,intersection_size);
            rhs_sim = (double)(intersection_size)/(double)(union_size);
            boxes.push_back(rhs_shrink_box);
            pq.push({rhs_sim, ++ind});

        }

        vector<QueryBox> corner_boxes = generateCornerSkylineRanges(box, lhs_expands, rhs_expands);

        for (const auto& corner_box : corner_boxes) {
            for(int i=0;i<d;i++){
                cout<<"{"<<corner_box.bounds[i].first<<", "<< corner_box.bounds[i].second<<"}, ";
            }
            vector<Point*> skyline = paperAlgorithmWithRangeTree(root, corner_box);
            bruteForceSkyline(corner_box);
            // bruteForceSkyline(corner_box);
            for (Point* pt : skyline) {
                QueryBox neighbor = box;
                for (int i = 0; i < d; ++i) {
                    neighbor.bounds[i].first = min(neighbor.bounds[i].first, pt->coords[i]);
                    neighbor.bounds[i].second = max(neighbor.bounds[i].second, pt->coords[i]);
                }
            }
        }
    }
    
}

int main() {
    int n;
    // ifstream cin("skyline_data_random.txt");
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

    RangeTreeNode *root = build_tree(point_ptrs, 0, d+t-1, d + t);
    vector<RangeTreeNode *> trees(d);
    for(int i=0;i<d;i++){
        trees[i]=build_tree(point_ptrs, 0, i, d);
    }
    QueryBox qbox;
    qbox.bounds.resize(d);
    cout << "Enter query box (range space) [" << d << " intervals]:\n";
    for (int i = 0; i < d; i++) {
        cin >> qbox.bounds[i].first >> qbox.bounds[i].second;
    }

    // paperAlgorithmWithRangeTree(root, qbox);
    // bruteForceSkyline(qbox);
    FeatureBox fullBox;
    fullBox.bounds.resize(t);
    for (int i = 0; i < t; i++) {
        fullBox.bounds[i] = { -numeric_limits<double>::infinity(), numeric_limits<double>::infinity() };
    }
    vector<Point*> lhs_expands, rhs_expands;
    for (int dim = 0; dim < d; ++dim) {
        // side_expand(trees[dim], qbox, 0, d, dim, true, true, lhs_expands);
        // side_expand(trees[dim], qbox, 0, d, dim, false, true, rhs_expands);
        vector<Point*> out;
        int size=0;
        side_expand_2(trees[dim], qbox, fullBox, 0, d, d, t,  out, size);
        for(auto p: out){
            if(p->coords[dim]<qbox.bounds[dim].first)lhs_expands.push_back(p);
            if(p->coords[dim]>qbox.bounds[dim].second)rhs_expands.push_back(p);
        }
    }

    vector<QueryBox> corner_boxes = generateCornerSkylineRanges(qbox, lhs_expands, rhs_expands);

    for (const auto& corner_box : corner_boxes) {
        for(int i=0;i<d;i++){
            cout<<"{"<<corner_box.bounds[i].first<<", "<< corner_box.bounds[i].second<<"}, ";
        }
        vector<Point*> skyline = paperAlgorithmWithRangeTree(root, corner_box);
        bruteForceSkyline(corner_box);
        // bruteForceSkyline(corner_box);
        for (Point* pt : skyline) {
            QueryBox neighbor = qbox;
            for (int i = 0; i < d; ++i) {
                neighbor.bounds[i].first = min(neighbor.bounds[i].first, pt->coords[i]);
                neighbor.bounds[i].second = max(neighbor.bounds[i].second, pt->coords[i]);
            }
        }
    }
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
15 15
45 45
15 15
10 10
*/

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

/*
4
2 3
10 10   5 10 15
15 15   4 12 16
12 12   6 9 14
18 18   5 10 13
10 20
10 20
*/