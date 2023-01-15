#include <math.h>
#include <stdio.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "minimal_free_space_constraints.h"

using Points = std::vector<Point>;
using Balls = std::vector<Ball>;
using PointsPerPoke = std::map<int, Points>;

PointsPerPoke read_freespace_pc(std::ifstream &input) {
    PointsPerPoke ret;
    while (input) {
        int num_pokes = -1, num_points = -1;
        input >> num_pokes >> num_points;
        if (num_pokes == -1) {
            break;
        }

        std::cout << "freespace: poke " << num_pokes << " with " << num_points
                  << " points\n";
        Points pts;

        for (int i = 0; i < num_points; i++) {
            float x, y, z;
            input >> x >> y >> z;

            Point pt = {};
            // freespace boundary is not necessarily part of physical surfaces
            pt.jump = true;
            pt.p = {x, y, z};
            pts.push_back(pt);
        }
        for (int i = 0; i < num_points; i++) {
            float x, y, z;
            input >> x >> y >> z;
            pts[i].n = {x, y, z};
        }
        ret[num_pokes] = pts;
    }

    return ret;
}

void read_contact_pc(std::ifstream &input, PointsPerPoke& points_per_poke) {
    while (input) {
        int num_pokes = -1, num_points = -1;
        input >> num_pokes >> num_points;
        if (num_pokes == -1) {
            break;
        }

        std::cout << "contacts: poke " << num_pokes << " with " << num_points
                  << " points\n";

        auto& pts = points_per_poke[num_pokes];

        for (int i = 0; i < num_points; i++) {
            float x, y, z;
            int in_contact;
            input >> x >> y >> z >> in_contact;
            if (in_contact != 1) {
                continue;
            }

            Point pt = {};
            // contact points are known to be on physical surfaces
            pt.jump = false;
            pt.p = {x, y, z};
            // these experiments are always pushing along +x, so estimate normal of -x
            pt.n = {-1, 0, 0};
            pts.push_back(pt);
        }
    }
}

int main(int argc, char **argv) {
    std::string freespace_file(argv[1]);
    std::string contact_file(argv[2]);
    std::string output_file(argv[3]);

    std::ifstream freespace_stream{freespace_file};
    auto points_per_poke = read_freespace_pc(freespace_stream);
    std::ifstream contact_stream{contact_file};
    read_contact_pc(contact_stream, points_per_poke);

    float r_max = 1.0f;
    float min_angle = 0.1f;
    int smoothing_neighborhood_size = 4;
    float smoothing_factor = 0.1f;
    float ball_delta = 0.02f;
    int max_balls_to_use = 0;
    // DESCRIPTION OF PARAMETERS
    // min_angle:
    //     Angle criterion described in Ma et al. §7.2.
    // r_max:
    //     Clipping distance. Should be larger than your scene radius.
    // smoothing_neighborhood, smoothing_factor:
    //     Smoothing filter applied to medial axis distances, based on
    //     k nearest neighbors.
    // ball_delta:
    //     Ball delta radius used to compute simplifying min-cover,
    //     as described in Ma et al. §7.2.
    // max_balls_to_use:
    //     Can be set to non-zero to additionally constrain
    //     min-cover simplification to use a given number of
    //     balls. Set to 0 to only use provided ball_delta.

    std::ofstream output_stream(output_file);

    for (auto& kv: points_per_poke) {
        int num_pokes = kv.first;
        Points& points = kv.second;

        std::cout << "compute ball: poke " << num_pokes << " with " << points.size() << " points" << std::endl;
        compute_medial_axis_distance(points, r_max, min_angle,
                                    smoothing_neighborhood_size, smoothing_factor);

        Balls balls;
        for (auto pt : points) {
            Ball ball;
            ball.c = pt.p + pt.mad * pt.n;
            ball.p = pt.p;
            ball.r = pt.mad;
            balls.push_back(ball);
        }

        // cover.is_used[i] will now be true if ball i is part
        // of the simplified cover of free space.

        // cover.update_min_cover can be called repeatedly after
        // this point with a smaller ball_delta for fast updates.
        // If called with a larger ball_delta, it will be slower.
        Cover cover;
        cover.initialize(&balls[0], (int)balls.size());
        cover.update_min_cover(ball_delta, max_balls_to_use);

        size_t balls_used = 0;
        for (size_t i = 0; i < balls.size(); i++) {
            const auto& b = balls[i];
            if (cover.is_used[i]) {
                balls_used += 1;
            }
        }

        std::cout << "poke " << num_pokes << " with " << balls_used << " balls used" << std::endl;
        output_stream << num_pokes << " " << balls_used << std::endl;

        for (size_t i = 0; i < balls.size(); i++) {
            const auto& b = balls[i];
            if (cover.is_used[i]) {
                output_stream << b.c.x << ' ' << b.c.y << ' ' << b.c.z << ' ' <<
                        b.p.x << ' ' << b.p.y << ' ' << b.p.z << ' ' << b.r << std::endl;
            }
        }
    }
}
