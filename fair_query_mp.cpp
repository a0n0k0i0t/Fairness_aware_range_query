#include<bits/stdc++.h>
using namespace std;
struct Tuple {
    vector<int>row;
};

vector<int>d_attr;
struct Point{
    Tuple T;
};

struct RangeTreeNode {
    int level;
    int val;
    RangeTreeNode *left;
    RangeTreeNode *right;
    RangeTreeNode *next_tree;
    Point min_point;
    Point* curr_pt;
    RangeTreeNode(int lev, int v, RangeTreeNode *l, RangeTreeNode *r, RangeTreeNode *next, Point p, Point *curr)
    :level{lev}, val{v}, left{l}, right{r}, next_tree{next}, min_point{p}, curr_pt{curr}{}
};
bool min_all(Point& p1, Point& p2, int t){
    for(int i=t-1;i>=0;i--){
        if(p1.T.row[i]<p2.T.row[i]) return 1;
        if(p1.T.row[i]>p2.T.row[i]) return 0;
    }
    return 0;
}
RangeTreeNode* build_tree(vector<Point>& pts, int l, int r, int level, int total_levels, int t, int d, vector<vector<int>>&sorted_indices){
    // if(pts.size()==0) return NULL;
    if(l>r){
        return NULL;
    }
    // sort(pts.begin()+l,pts.begin()+r+1,[&](Point a, Point b){
    //     if(level<d){
    //         return a.T.row[d_attr[level]] < b.T.row[d_attr[level]];
    //     }
    //     else{
    //         return a.T.row[level-d] < b.T.row[level-d];
    //     }
    // });

    int mid = (l+r)/2;
    // int split_val = (level<d) ? pts[mid].T.row[d_attr[level]] : pts[mid].T.row[level-d];
    int dim = (level<d) ? d_attr[level]:level - d;
    int split_val = (level<d) ? pts[sorted_indices[dim][mid]].T.row[dim] : pts[sorted_indices[dim][mid]].T.row[dim];
    // cout<<split_val<<" "<<level<<"\n";
    Point * curr = &pts[sorted_indices[dim][mid]];
    RangeTreeNode *left = build_tree(pts, l, mid-1, level, total_levels, t, d, sorted_indices);
    RangeTreeNode *right= build_tree(pts, mid+1, r, level, total_levels, t, d, sorted_indices);
    RangeTreeNode *next = NULL;
    if(level !=total_levels-1)
    next = build_tree(pts, l, r, level+1, total_levels, t, d, sorted_indices);

    vector<Point> candidates;
    candidates.push_back(pts[mid]);
    if (left) candidates.push_back(left->min_point);
    if (right) candidates.push_back(right->min_point);
    if (next) candidates.push_back(next->min_point);

    Point min_p = candidates[0];
    for (auto& p : candidates) {
        if (min_all(p, min_p, t))
            min_p = p;
    }

    return new RangeTreeNode(level, split_val, left, right, next, min_p,curr);
}

bool in_box(Point& p, vector<pair<int, int>>& box, int d, int t) {
    for (int i = 0; i < d + t; ++i) {
        if (p.T.row[i] < box[i].first || p.T.row[i] > box[i].second)
            return false;
    }
    return true;
}


void range_query(RangeTreeNode* node, vector<pair<int, int>>& query, vector<pair<int, int>>& box,  int level, int total_levels, set<Point*>& result,
    int t, int d) {
    if (!node) return;
    int dim = (level<d) ? level:level - d;
    int l = (level<d) ? query[dim].first:box[dim].first;
    int r = (level<d) ? query[dim].second:box[dim].second;
    // cout<<node->val<<" "<<l<<" "<<r<<" "<<level<<"\n";
    if(node->val<l){
        range_query(node->right, query, box, level, total_levels, result, t, d);
    } 
    else if(node->val>r){
        range_query(node->left, query, box, level, total_levels, result, t, d);
    } 
    else{
        if(level==total_levels-1){
            if (node->val>=l && node->val<=r) {
                result.insert(node->curr_pt);
            }
            range_query(node->left, query, box, level, total_levels, result, t, d);
            range_query(node->right, query, box, level, total_levels, result, t, d);
        } 
        else{
            range_query(node->next_tree, query, box, level + 1, total_levels, result, t, d);
            range_query(node->left, query, box, level, total_levels, result, t, d);
            range_query(node->right, query, box, level, total_levels, result, t, d);
        }
    }
}


void dfs(RangeTreeNode *node, RangeTreeNode *root,int level, int total_levels){
    if(node == NULL){
        return;
    }
    while(root!=NULL){
        dfs(node->left,root,level,total_levels);
        cout<<node->val<<" ";
        dfs(node->right,root,level,total_levels);
        if(node==root){
            root=root->next_tree;
            node=root;
            cout<<"\n";
        }
        else{
            return;
        }
    }
}
int main() {
    int d = 2;
    int t = 3;

    d_attr = {0, 2};

    vector<Point> points = {
        {{ {1, 100, 50} }},
        {{ {2, 150, 80} }},
        {{ {3, 120, 70} }},
        {{ {1, 80,  60} }}
    };
    int n=points.size();
    vector<vector<int>> sorted_indices(t,vector<int>(n));
    for (int i = 0; i < t; ++i) {
        iota(sorted_indices[i].begin(), sorted_indices[i].end(), 0);
        sort(sorted_indices[i].begin(), sorted_indices[i].end(), [&](int a, int b) {
            return points[a].T.row[i] < points[b].T.row[i];
        });
    }

    RangeTreeNode* root = build_tree(points, 0, points.size() - 1, 0, d + t, t, d, sorted_indices);
    set<Point*>result;
    vector<pair<int, int>> query = {{1, 2}, {50, 80}};
    vector<pair<int, int>> box   = {{-1e9, 1e9}, {-1e9, 1e9}, {-1e9, 1e9}};
    range_query(root, query, box, 0, d+t, result, t, d);
    
    for(auto ptr : result){
        Point P =*ptr;
        for(int i=0;i<t;i++){
            cout<<P.T.row[i]<<" ";
        }
        cout<<"\n";
    }
    // cout<<root->val;
    // dfs(root,root,0,d+t);
}