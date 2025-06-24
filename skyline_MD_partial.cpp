#include <iostream>
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

const int d = 2; // Range dimensions
const int t = 2; // Feature dimensions

typedef bg::model::point<double, d, bg::cs::cartesian> RangePoint;
typedef bg::model::box<RangePoint> RangeBox;
typedef pair<RangePoint, size_t> Value; // point + index

struct DataPoint {
    RangePoint range_coords;
    vector<double> features;
    int index;
};

bool dominates(const DataPoint& p, const DataPoint& q) {
    bool strictlyBetter = false;
    for (int i = 0; i < t; i++) {
        if (p.features[i] > q.features[i]) return false; // worse in this feature
        if (p.features[i] < q.features[i]) strictlyBetter = true;
    }
    return strictlyBetter;
}

int main() {
    int n;
    cout << "Enter number of points: ";
    cin >> n;

    vector<DataPoint> points(n);

    cout << "Enter range coords (dim=" << d << ") and features (dim=" << t << ") for each point:\n";
    for (int i = 0; i < n; i++) {
        vector<double> rc(d);
        for (int j = 0; j < d; j++) cin >> rc[j];
        RangePoint rp(rc[0], rc[1]);
        vector<double> features(t);
        for (int j = 0; j < t; j++) cin >> features[j];
        points[i] = {rp, features, i};
    }

    // Build R-tree
    bgi::rtree<Value, bgi::quadratic<16>> rtree;
    for (const auto& p : points) {
        rtree.insert(make_pair(p.range_coords, p.index));
    }

    // Read query box
    cout << "Enter query box bounds (lower and upper) for each dimension:\n";
    vector<double> lower(d), upper(d);
    for (int i = 0; i < d; i++) cin >> lower[i] >> upper[i];
    RangeBox query_box(RangePoint(lower[0], lower[1]), RangePoint(upper[0], upper[1]));

    // Query R-tree for points inside the box
    vector<Value> result_s;
    rtree.query(bgi::intersects(query_box), back_inserter(result_s)); // <-- Important fix here!

    // Collect candidates
    vector<DataPoint> candidates;
    for (auto& val : result_s) {
        candidates.push_back(points[val.second]);
    }

    // Compute skyline among candidates
    vector<DataPoint> skyline;
    for (const auto& p : candidates) {
        bool dominated_flag = false;
        for (const auto& q : candidates) {
            if (p.index != q.index && dominates(q, p)) {
                dominated_flag = true;
                break;
            }
        }
        if (!dominated_flag) skyline.push_back(p);
    }

    // Output skyline
    cout << "Skyline points in query box:\n";
    for (const auto& p : skyline) {
        cout << "Point " << p.index << ": Range [ ";
        cout << p.range_coords.get<0>() << " " << p.range_coords.get<1>() << " ] Features [ ";
        for (auto f : p.features) cout << f << " ";
        cout << "]\n";
    }

    return 0;
}
