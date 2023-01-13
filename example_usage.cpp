#include <stdio.h>
#include <math.h>
#include <vector>


#include "minimal_free_space_constraints.h"

int main()
{
    // Generate a point cloud of a circle, with normals pointing
    // inward to the center. Noise is added to the radius.
    std::vector<Point> points;
    for (int i = 0; i < 250; i++)
    {
        float t = 2.0f*3.14f*(i+0.5f)/250.0f;
        float r = 1.0f + 0.2f*(rand() % 10000)/10000.0f;
        float nx = cosf(t);
        float ny = sinf(t);
        float x = r*nx;
        float y = r*ny;

        Point pt = {};
        pt.jump = false;
        pt.p = {x, y, 0.0f};
        pt.n = {-nx, -ny, 0.0f};
        points.push_back(pt);
    }

    float r_max = 10.0f;
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

    compute_medial_axis_distance(points,
        r_max,
        min_angle,
        smoothing_neighborhood_size,
        smoothing_factor);

    std::vector<Ball> balls;
    for (auto pt : points)
    {
        Ball ball;
        ball.c = pt.p + pt.mad*pt.n;
        ball.p = pt.p;
        ball.r = pt.mad;
        balls.push_back(ball);
    }

    Cover cover;
    cover.initialize(&balls[0], (int)balls.size());
    cover.update_min_cover(ball_delta, max_balls_to_use);

    // cover.is_used[i] will now be true if ball i is part
    // of the simplified cover of free space.

    // cover.update_min_cover can be called repeatedly after
    // this point with a smaller ball_delta for fast updates.
    // If called with a larger ball_delta, it will be slower.

    // Save the balls for visualization in visualize.py
    {
        FILE *file = fopen("balls.txt", "w");
        for (size_t i = 0; i < balls.size(); i++)
        {
            Ball b = balls[i];
            if (cover.is_used[i])
                fprintf(file, "%g %g %g %g %g %g %g\n", b.c.x, b.c.y, b.c.z, b.p.x, b.p.y, b.p.z, b.r);
        }
        fclose(file);
    }
}
