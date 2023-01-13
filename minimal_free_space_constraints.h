#pragma once
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <limits>
#include <vector>
// You can provide your own vec3 type, but it must expose
// the following interface.
struct vec3 {
    float x, y, z;
};
vec3 operator-(vec3 a, vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
vec3 operator+(vec3 a, vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
vec3 operator*(float k, vec3 a) { return {k * a.x, k * a.y, k * a.z}; }
vec3 operator*(vec3 a, float k) { return {k * a.x, k * a.y, k * a.z}; }
vec3 operator/(vec3 a, float k) { return {a.x / k, a.y / k, a.z / k}; }
float m_dot(vec3 a, vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
float m_length(vec3 a) { return sqrtf(m_dot(a, a)); }
vec3 m_normalize(vec3 a) { return a / m_length(a); }

// You can define your own Point and Ball types, but they must
// adhere to the following definitions.
struct Point {
    vec3 p;     // position
    vec3 n;     // normal
    bool jump;  // true if the point is not on a visible surface
    float mad;  // distance to medial axis
};

struct Ball {
    vec3 c;     // ball center
    vec3 p;     // associated point on point cloud
    float r;    // ball radius
    bool jump;  // true if the associated point on point cloud is not on a
                // visible surface
};

// Instead of this, you can include loguru for better logging
#define INFO NULL
#define LOG_F(type, ...) \
    printf("..");        \
    printf(__VA_ARGS__); \
    printf("\n");
#define LOG_SCOPE_F(type, ...) \
    printf("\n");              \
    printf(__VA_ARGS__);       \
    printf("\n");
#define LOG_SCOPE_FUNCTION(type, ...)

// These #defines configure the constraint generator to be fast
#define MIN_COVER_PRIORITIZE_ZEROS
#define MIN_COVER_SIMPLIFY_ZEROS
#define INCREMENTAL_MIN_COVER

void _compute_medial_axis_distance(std::vector<Point> &points, float r_max,
                                   float min_angle) {
    LOG_SCOPE_F(INFO, "Computing medial axis distances");
    LOG_F(INFO, "%llu points", points.size());
    LOG_F(INFO, "r_max = %g", r_max);
    float cos_min_angle = cosf(min_angle);
    for (size_t i = 0; i < points.size(); i++) {
        float r = r_max;
        for (size_t j = 0; j < points.size(); j++) {
            if (i == j) continue;
            vec3 d = points[j].p - points[i].p;
            float d_dot_n = m_dot(d, points[i].n);
            if (d_dot_n <= 0.0f) continue;
            float r_ij = 0.5f * m_dot(d, d) / d_dot_n;
            if (r_ij > r) continue;
            float cos_a = m_dot(
                points[i].n,
                m_normalize(points[i].p + r_ij * points[i].n - points[j].p));
            if (cos_a > cos_min_angle) continue;
            r = r_ij;
        }
        points[i].mad = r;
    }
}

void _smooth_medial_axis_distance(std::vector<Point> &points, float r_max,
                                  float min_angle,
                                  int smoothing_neighborhood_size,
                                  float smoothing_factor) {
    LOG_SCOPE_F(INFO, "Smooth medial axis distances");
    LOG_F(INFO, "Computing local expected medial distances.");
    struct keyval {
        float d2;
        size_t j;
    };
    static std::vector<keyval> neighbors;
    static std::vector<float> local_expected_medial_distances;
    local_expected_medial_distances.resize(points.size());
    neighbors.reserve(points.size());
    float cos_min_angle = cosf(min_angle);
    for (size_t i = 0; i < points.size(); i++) {
        neighbors.resize(0);
        for (size_t j = 0; j < points.size(); j++) {
            if (i == j) continue;
            vec3 pij = points[j].p - points[i].p;
            neighbors.push_back({.d2 = m_dot(pij, pij), .j = j});
        }

        vec3 pi = points[i].p;
        vec3 ni = points[i].n;

        // find expected radius based on k nearest neighbors
        float expected = 0.0f;
        for (int ik = 0; ik < smoothing_neighborhood_size; ik++) {
            // find the ik'th nearest neighbor
            size_t best_n = 0;
            float best_d2 = std::numeric_limits<float>::max();
            for (size_t n = 0; n < neighbors.size(); n++) {
                if (neighbors[n].d2 < best_d2) {
                    best_d2 = neighbors[n].d2;
                    best_n = n;
                }
            }
            size_t j = neighbors[best_n].j;
            neighbors[best_n] = neighbors[neighbors.size() - 1];
            neighbors.pop_back();

            vec3 cj = points[j].p + points[j].mad * points[j].n;
            float r = smoothing_factor * m_dot(cj - pi, ni);
            if (r > expected) expected = r;
        }
        local_expected_medial_distances[i] = expected;
    }
    LOG_F(INFO, "Recomputing medial distances using local expectation.");
    for (size_t i = 0; i < points.size(); i++) {
        float r = r_max;
        float r_min = local_expected_medial_distances[i];
        for (size_t j = 0; j < points.size(); j++) {
            if (i == j) continue;
            vec3 d = points[j].p - points[i].p;
            float d_dot_n = m_dot(d, points[i].n);
            if (d_dot_n <= 0.0f) continue;
            float r_ij = 0.5f * m_dot(d, d) / d_dot_n;
            if (r_ij > r || r_ij < r_min) continue;
            float cos_a = m_dot(
                points[i].n,
                m_normalize(points[i].p + r_ij * points[i].n - points[j].p));
            if (cos_a > cos_min_angle) continue;
            r = r_ij;
        }
        points[i].mad = r;
    }
}

void compute_medial_axis_distance(std::vector<Point> &points, float r_max,
                                  float min_angle,
                                  int smoothing_neighborhood_size,
                                  float smoothing_factor) {
    _compute_medial_axis_distance(points, r_max, min_angle);
    if (smoothing_neighborhood_size > 1)
        _smooth_medial_axis_distance(points, r_max, min_angle,
                                     smoothing_neighborhood_size,
                                     smoothing_factor);
}

struct Cover {
    Ball *all_balls;
    int all_balls_count;
    Ball *used_balls;
    int used_balls_count;
    bool *is_used;
    bool *is_covered;
    int *is_covered_by;
    std::vector<int> *neighbors;
    std::vector<int> *potential_neighbors;
    float current_neighborhood_delta = 0.0f;

    Cover() {
        all_balls = NULL;
        all_balls_count = 0;
        used_balls = NULL;
        used_balls_count = 0;
        is_used = NULL;
        is_covered = NULL;
        is_covered_by = NULL;
        neighbors = NULL;
        potential_neighbors = NULL;
        current_neighborhood_delta = 0.0f;
    }

    void initialize(Ball *_all_balls, int _all_balls_count);
    void initialize_neighbors(float delta);
    void update_min_cover(float delta, int max_used_balls_count);
    void use_all_balls();
};

void Cover::initialize(Ball *_all_balls, int _all_balls_count) {
    LOG_SCOPE_FUNCTION(INFO);
    if (used_balls) free(used_balls);
    if (is_used) free(is_used);
    if (is_covered) free(is_covered);
    if (is_covered_by) free(is_covered_by);
    // if (neighbors) free(neighbors);
    // if (potential_neighbors) free(potential_neighbors);
    all_balls = _all_balls;
    assert(all_balls);
    all_balls_count = _all_balls_count;
    used_balls = new Ball[all_balls_count];
    assert(used_balls);
    used_balls_count = all_balls_count;
    is_used = new bool[all_balls_count];
    assert(is_used);
    is_covered = new bool[all_balls_count];
    assert(is_covered);
    is_covered_by = new int[all_balls_count];
    assert(is_covered_by);
    neighbors = new std::vector<int>[all_balls_count];
    assert(neighbors);
    potential_neighbors = new std::vector<int>[all_balls_count];
    assert(potential_neighbors);
    for (int bi = 0; bi < all_balls_count; bi++) {
        is_covered[bi] = false;
        is_covered_by[bi] = bi;
        is_used[bi] = true;
        used_balls[bi] = all_balls[bi];
    }
    current_neighborhood_delta = 0.0f;
}

#ifdef INCREMENTAL_MIN_COVER
// Here, for each ball bi, we maintain a list of coverable balls: neighbors[bi].
// Each list is trimmed at each call to update_min_cover to remain correct for
// the new given delta. We assume that the given delta is smaller or equal to
// the previous delta and that the balls are always shrinking. This ensures
// that a ball cannot acquire new coverables; i.e. each neighbors list only
// needs to be trimmed.
// Second, instead of recalculating the number of coverables contributable by
// each ball in each iteration of the greedy loop, we compute this once in
// update_min_cover, and update it efficiently in each iteration based on
// tracking which balls could cover the newly chosen ball.
//
// Actually, the above assumption does not hold. Consider balls approaching a
// concave corner from different directions. Initially, they are not neighbors
// but eventually they are. I think it's sufficient to test the balls of points
// that are within ball_delta of p, for potential new neighbors.
//
// Note: We probably should use a seperate delta for potential_neighbors of
// surface points. This is noticeable if delta is increased, then decreased,
// as we do not re-trim potential_neighbors.
void Cover::initialize_neighbors(float delta) {
    LOG_SCOPE_F(INFO, "Initializing neighbors for %d balls", all_balls_count);
    assert(all_balls);
    assert(all_balls_count > 0);
    assert(neighbors);
    current_neighborhood_delta = delta;
    int sum_num_neighbors = 0;
    int max_num_neighbors = 0;
    for (int bi = 0; bi < all_balls_count; bi++) {
        neighbors[bi].resize(0);
        for (int bj = 0; bj < all_balls_count; bj++) {
            if (bi == bj) continue;
            float ri = all_balls[bi].r;
            float rj = all_balls[bj].r;
            float d = m_length(all_balls[bj].c - all_balls[bi].c);
            if (d + rj <= ri + delta) neighbors[bi].push_back(bj);
        }
        int num_neighbors = (int)neighbors[bi].size();
        sum_num_neighbors += num_neighbors;
        if (num_neighbors > max_num_neighbors)
            max_num_neighbors = num_neighbors;
    }
    LOG_F(INFO, "Average neighbors: %.2f",
          (float)sum_num_neighbors / all_balls_count);
    LOG_F(INFO, "Maximum neighbors: %d", max_num_neighbors);

    float delta2 = delta * delta;
    int sum_num_potential_neighbors = 0;
    int max_num_potential_neighbors = 0;
    for (int bi = 0; bi < all_balls_count; bi++) {
        potential_neighbors[bi].resize(0);
        for (int bj = 0; bj < all_balls_count; bj++) {
            if (bi == bj) continue;
            vec3 pipj = all_balls[bj].p - all_balls[bi].p;
            float d2 = m_dot(pipj, pipj);
            if (d2 > delta2) continue;
            bool already_neighbors = false;
            for (auto bk : neighbors[bi])
                if (bk == bj) already_neighbors = true;
            if (already_neighbors) continue;
            potential_neighbors[bi].push_back(bj);
        }
        int num_potential_neighbors = (int)potential_neighbors[bi].size();
        sum_num_potential_neighbors += num_potential_neighbors;
        if (num_potential_neighbors > max_num_potential_neighbors)
            max_num_potential_neighbors = num_potential_neighbors;
    }
    LOG_F(INFO, "Average potential neighbors: %.2f",
          (float)sum_num_potential_neighbors / all_balls_count);
    LOG_F(INFO, "Maximum potential neighbors: %d", max_num_potential_neighbors);
}

void Cover::update_min_cover(float delta, int max_used_balls_count = 0) {
    LOG_SCOPE_F(INFO, "Updating min-cover");
    assert(all_balls_count > 0);
    assert(is_used);
    assert(is_covered);
    assert(neighbors);

    LOG_F(INFO, "delta = %g", delta);
    LOG_F(INFO, "max_used_balls_count = %d", max_used_balls_count);

    if (delta > current_neighborhood_delta) {
        LOG_F(INFO,
              "delta > current_neighborhood_delta, recomputing neighbors from "
              "scratch.");
        initialize_neighbors(delta);
    }
    current_neighborhood_delta = delta;

    for (int i = 0; i < all_balls_count; i++) {
        is_covered[i] = false;
        is_used[i] = false;
    }

    if (max_used_balls_count == 0) max_used_balls_count = all_balls_count;

    static std::vector<int> num_coverable;
    static std::vector<std::vector<int>> is_coverable_by;
    num_coverable.reserve(all_balls_count);
    is_coverable_by.resize(all_balls_count);

    for (int bi = 0; bi < all_balls_count; bi++) {
        is_coverable_by[bi].resize(0);
        is_coverable_by[bi].push_back(bi);
    }

    LOG_F(INFO, "Trimming neighbor lists...");
    int sum_removed = 0;
    for (int bi = 0; bi < all_balls_count; bi++) {
        int num_neighbors_bi = (int)neighbors[bi].size();
        int *neighbors_bi = &neighbors[bi][0];
        for (int n = 0; n < num_neighbors_bi; n++) {
            int bj = neighbors_bi[n];
            float d = m_length(all_balls[bj].c - all_balls[bi].c);
            if (d + all_balls[bj].r > all_balls[bi].r + delta) {
                neighbors_bi[n] = neighbors_bi[num_neighbors_bi - 1];
                num_neighbors_bi--;
                n--;
                sum_removed++;
            } else {
                is_coverable_by[bj].push_back(bi);
            }
        }
        neighbors[bi].resize(num_neighbors_bi);

        int num_potential_neighbors_bi = (int)potential_neighbors[bi].size();
        int *potential_neighbors_bi = &potential_neighbors[bi][0];
        for (int n = 0; n < num_potential_neighbors_bi; n++) {
            int bj = potential_neighbors_bi[n];
            float d = m_length(all_balls[bj].c - all_balls[bi].c);
            if (d + all_balls[bj].r <= all_balls[bi].r + delta) {
                neighbors[bi].push_back(bj);
                potential_neighbors_bi[n] =
                    potential_neighbors_bi[num_potential_neighbors_bi - 1];
                num_potential_neighbors_bi--;
                n--;
            }
        }
        potential_neighbors[bi].resize(num_potential_neighbors_bi);

        num_coverable[bi] =
            (int)neighbors[bi].size() + 1;  // ball can also cover itself
    }
    LOG_F(INFO, "Removed %d neighbors.", sum_removed);

#ifdef MIN_COVER_PRIORITIZE_ZEROS
    LOG_F(INFO,
          "ENABLED: Prioritize zeros. Terminating surface points will always "
          "be included.");
#else
    LOG_F(INFO,
          "DISABLED: Prioritize zeros. Terminating surface points will not "
          "necessarily be included.");
#endif

#ifdef MIN_COVER_SIMPLIFY_ZEROS
    LOG_F(INFO,
          "ENABLED: Simplify zeros. Terminating surface points will be subject "
          "to min-cover simplification.");
#else
    LOG_F(INFO,
          "DISABLED: Simplify zeros Terminating surface points will not be "
          "subject to min-cover simplification.");
#endif

    used_balls_count = 0;
    while (used_balls_count < max_used_balls_count) {
        int best_num_coverable = -1;
        int best_bi = 0;
        for (int bi = 0; bi < all_balls_count; bi++) {
            // If we skip is_covered check, we use all zero-terminated
            // constraints otherwise, we allow them to be simplified

            if (is_used[bi]) continue;
#ifdef MIN_COVER_PRIORITIZE_ZEROS
#ifdef MIN_COVER_SIMPLIFY_ZEROS
            if (is_covered[bi]) continue;
#endif
            if (all_balls[bi].r == 0.0f && !all_balls[bi].jump) continue;
            if (all_balls[bi].r == 0.0f) {
                best_num_coverable = num_coverable[bi];
                best_bi = bi;
                break;
            }
#endif
#ifndef MIN_COVER_SIMPLIFY_ZEROS
            if (is_covered[bi]) continue;
#endif
            if (num_coverable[bi] > best_num_coverable) {
                best_num_coverable = num_coverable[bi];
                best_bi = bi;
            }
        }

        if (best_num_coverable <= 0) break;

        int bi = best_bi;
        is_used[bi] = true;
        is_covered[bi] = true;
        used_balls[used_balls_count++] = all_balls[bi];

        for (int bj : neighbors[bi]) {
#ifdef MIN_COVER_PRIORITIZE_ZEROS
#ifndef MIN_COVER_SIMPLIFY_ZEROS
            if (all_balls[bj].r == 0.0f) continue;
#endif
#endif
            if (is_used[bj]) continue;
            if (is_covered[bj]) continue;
            float d = m_length(all_balls[bj].c - all_balls[bi].c);
            if (d + all_balls[bj].r <= all_balls[bi].r + delta) {
                is_covered[bj] = true;
                is_covered_by[bj] = bi;
                for (auto bi : is_coverable_by[bj]) num_coverable[bi]--;
            }
        }
    }
    LOG_F(INFO, "Used %d balls.", used_balls_count);
}
#else
void Cover::initialize_neighbors(float delta) {
    printf("Updating neighbors for cover (%d balls) ... ", all_balls_count);
    assert(all_balls);
    assert(all_balls_count > 0);
    assert(neighbors);
    // todo: can probably estimate more tight neighborhood
    // Neighborhood of B_i is the set of balls that B_i could at most cover.
    int sum_num_neighbors = 0;
    for (int bi = 0; bi < all_balls_count; bi++) {
        neighbors[bi].resize(0);
        for (int bj = 0; bj < all_balls_count; bj++) {
            if (bi == bj) continue;
            float ri = all_balls[bi].r;
            float rj = all_balls[bj].r;
            float d = m_length(all_balls[bj].c - all_balls[bi].c);
            if (d <= ri + delta) neighbors[bi].push_back(bj);
        }
        sum_num_neighbors += (int)neighbors[bi].size();
    }
    printf("done. Average neighbor count %.0f.\n",
           (float)sum_num_neighbors / all_balls_count);
}

void Cover::update_min_cover(float delta, int max_used_balls_count = 0) {
    assert(all_balls_count > 0);
    assert(is_used);
    assert(is_covered);
    assert(neighbors);

    for (int i = 0; i < all_balls_count; i++) {
        is_covered[i] = false;
        is_used[i] = false;
    }

    if (max_used_balls_count == 0) max_used_balls_count = all_balls_count;

#ifdef FAST_MIN_COVER
#pragma message("Computing fast min cover")
    printf("... Updating neighbors ... ");
    int sum_removed = 0;
    for (int bi = 0; bi < all_balls_count; bi++) {
        int num_neighbors_bi = (int)neighbors[bi].size();
        int *neighbors_bi = &neighbors[bi][0];
        for (int n = 0; n < num_neighbors_bi; n++) {
            int bj = neighbors_bi[n];
            float d = m_length(all_balls[bj].c - all_balls[bi].c);
            if (d + all_balls[bj].r > all_balls[bi].r + delta) {
                neighbors_bi[n] = neighbors_bi[num_neighbors_bi - 1];
                num_neighbors_bi--;
                n = n - 1;
                sum_removed++;
            }
        }
        neighbors[bi].resize(num_neighbors_bi);
    }
    printf("done. Removed %d neighbors.\n", sum_removed);
#endif

    used_balls_count = 0;
    while (used_balls_count < max_used_balls_count) {
        int best_num_covered = 0;
        int best_bi = 0;
        for (int bi = 0; bi < all_balls_count; bi++) {
            if (is_used[bi]) continue;
            if (is_covered[bi]) continue;

            int num_covered = 1;

#ifdef FAST_MIN_COVER
#pragma message("Computing fast min cover")
            int num_neighbors_bi = (int)neighbors[bi].size();
            int *neighbors_bi = &neighbors[bi][0];
            for (int n = 0; n < num_neighbors_bi; n++) {
                int bj = neighbors_bi[n];
                if (is_used[bj]) continue;
                if (is_covered[bj]) continue;
                float d = m_length(all_balls[bj].c - all_balls[bi].c);
                if (d + all_balls[bj].r <= all_balls[bi].r + delta)
                    num_covered++;
            }
#else
#pragma message("WARNING: Computing min cover using all balls")
            for (int bj = 0; bj < all_balls_count; bj++) {
                if (is_used[bj]) continue;
                if (is_covered[bj]) continue;
                float d = m_length(all_balls[bj].c - all_balls[bi].c);
                if (d + all_balls[bj].r <= all_balls[bi].r + delta)
                    num_covered++;
            }
#endif

            if (num_covered > best_num_covered) {
                best_num_covered = num_covered;
                best_bi = bi;
            }
        }

        if (best_num_covered == 0) break;

        int bi = best_bi;
        is_used[bi] = true;
        is_covered[bi] = true;
        used_balls[used_balls_count++] = all_balls[bi];

#ifdef FAST_MIN_COVER
#pragma message("Computing fast min cover")
        int num_neighbors_bi = (int)neighbors[bi].size();
        int *neighbors_bi = &neighbors[bi][0];
        for (int n = 0; n < num_neighbors_bi; n++) {
            int bj = neighbors_bi[n];
            if (is_used[bj]) continue;
            if (is_covered[bj]) continue;
            float d = m_length(all_balls[bj].c - all_balls[bi].c);
            if (d + all_balls[bj].r <= all_balls[bi].r + delta) {
                is_covered[bj] = true;
                is_covered_by[bj] = bi;
            }
        }
#else
#pragma message("WARNING: Computing min cover using all balls")
        for (int bj = 0; bj < all_balls_count; bj++) {
            if (is_used[bj]) continue;
            if (is_covered[bj]) continue;
            float d = m_length(all_balls[bj].c - all_balls[bi].c);
            if (d + all_balls[bj].r <= all_balls[bi].r + delta) {
                is_covered[bj] = true;
                is_covered_by[bj] = bi;
            }
        }
#endif
    }
}
#endif

void Cover::use_all_balls() {
    for (int i = 0; i < all_balls_count; i++) {
        is_used[i] = true;
        is_covered_by[i] = i;
        is_covered[i] = true;
        used_balls[i] = all_balls[i];
    }
    used_balls_count = all_balls_count;
}
