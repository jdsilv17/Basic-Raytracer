// BasicRayTracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "EngineMath.h"
#include "Defines.h"
#include <iostream>



struct Sphere {
    TVECTOR color = { 0.0f, 0.0f, 0.0f, 0.0f };
    TVECTOR center = { 0.0f, 0.0f, 0.0f, 1.0f };
    float radius = 0.0f;

    // material
    float specular_exponent = 0.0f;
    float reflective = 0.0f;
};

enum class ELightType : uint8_t {
    Point,
    Spot,
    Directional,
    Ambient
};

struct Light {
    TVECTOR postion = { 0.0f, 0.0f, 0.0f, 1.0f };
    TVECTOR direction = { 1.0f, 1.0f, 1.0f, 0.0f };
    TVECTOR color = { 1.0f, 1.0f, 1.0f, 1.0f };
    float intensity = 1.0f;
    ELightType type = ELightType::Point;
};

struct Scene {
    std::vector<Sphere> spheres;
    std::vector<Light> lights;
    TVECTOR cam_pos{0.f, 0.f, 0.f, 1.f};
} scene;

TVECTOR canvas_to_viewport(int32_t x, int32_t y)
{
    TVECTOR v;
    v.x = x * 1.77777f /*viewport_width*/ / static_cast<float>(canvas_width);
    v.y = y * 0.92f/*viewport_height*/ / static_cast<float>(canvas_height);
    v.z = static_cast<float>(projection_plane_z);
    v.w = 0.0f;

    return v;
}

// Determines if a ray from the origin location intersects the given sphere using the quadratic formula
bool ray_intersect_sphere(const TVECTOR origin, const TVECTOR ray_direction, const Sphere& sphere, float& t1, float& t2)
{
    /*
    *   O == origin, V == viewport, C = center, P = point in space, r = radius
    *   D = V - O // ray/distance vector from origin to viewport
    *   P = O + tD // lerp to find point on ray
    *   dot(P - C,P - C) = r^2 // point on a sphere
    * 
    *   dot(O + tD - C, O + tD - C) = r^2
    *   CO = O - C
    *   dot(CO + tD, CO + tD) = r^2
    *   dot(CO + tD, CO) + dot(CO + tD, tD) = r^2 // dont understand how we got here
    *   // Distributive property of dot product dot(U + V,W) = dot(W,U) + dot(W,V)
    *   dot(CO,CO) + dot(CO,tD) + dot(CO,tD) + dot(tD,tD) = r^2
    *   dot(tD,tD) + 2 * dot(CO,tD) + dot(CO,CO) = r^2
    *   t^2 * dot(D,D) + t * 2 * dot(CO,D) + dot(CO,CO) - r^2 = 0
    *   a = dot(D,D), b = 2 * dot(CO,D), c = dot(CO,CO) - r^2 // redefine some terms
    *   at^2 + bt + c = 0 // Quadratic equation
    *   Because we have a Quadratic equation we can use the quadratic formula to solve for t
    *   [t1, t2] = (-b +/- sqrt(b^2 + 4ac)) / 2a
    *   t1 and t2 being the two points of intersection of a ray through a sphere
    */
    float r = sphere.radius;
    TVECTOR CO_v = Vector_Sub(origin, sphere.center);

    // Quadratic formula
    float a = Vector_Dot(ray_direction, ray_direction);
    float b = 2.0f * Vector_Dot(CO_v, ray_direction);
    float c = Vector_Dot(CO_v, CO_v) - r * r;

    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0.0f) // no intersection/no solution for t
    {
        t1 = t2 = INFINITY;
        return false;
    }

    t1 = (-b + sqrt(discriminant)) / (2.0f * a);
    t2 = (-b - sqrt(discriminant)) / (2.0f * a);
    return true;
}

// returns the closest intersection on the given ray and stores the closest sphere intersected by the ray
float closest_intersection(const TVECTOR ray_origin, const TVECTOR ray_direction, const float near_plane, const float far_plane, Sphere& closest_sphere)
{
    float t1, t2, closest_t = far_plane;

    for (const Sphere& s : scene.spheres)
    {
        if (ray_intersect_sphere(ray_origin, ray_direction, s, t1, t2))
        {
            if (near_plane < t1 && t1 < far_plane && t1 < closest_t)
            {
                closest_t = t1;
                closest_sphere = s;
            }
            if (near_plane < t2 && t2 < far_plane && t2 < closest_t)
            {
                closest_t = t2;
                closest_sphere = s;
            }
        }
    }

    return closest_t;
}

float compute_lighting(TVECTOR point, TVECTOR normal, TVECTOR eye, float specular)
{
    float total_intensity = 0.0f;
    for (const Light& light : scene.lights)
    {
        if (light.type == ELightType::Ambient)
        {
            total_intensity += light.intensity; 
        }
        else
        {
            TVECTOR light_v = { 0.0f, 0.0f, 0.0f, 0.0f };
            float t_max = 0.0f;
            if (light.type == ELightType::Point)
            {
                light_v = Vector_Sub(light.postion, point);
                t_max = 1.0f;
            }
            else // if directional light
            {
                light_v = light.direction;
                t_max = INFINITY;
            }

            //Shadow check
            Sphere shadow_sphere;
            float shadow_t = closest_intersection(point, light_v, 0.001f, t_max, shadow_sphere);
            if (shadow_t != t_max) // if there was an intersection/there is something blocking the light
                continue;

            // Diffuse Reflection
            float n_dot_l = Vector_Dot(normal, light_v);
            if (n_dot_l > 0.0f)
                total_intensity += light.intensity * n_dot_l / (Vector_Length(normal) * Vector_Length(light_v));

            //Specular Reflection
            if (specular != -1.0f)
            {
                TVECTOR Ln = Vector_Scalar_Multiply(normal, n_dot_l); // projection of light vector onto surface normal
                TVECTOR Lp = Vector_Sub(light_v, Ln);
                TVECTOR reflection_v = Vector_Sub(Ln, Lp); //Vector_Sub(Vector_Scalar_Multiply(Ln, 2.0f), light_v);

                // the above is essentially the same as Vector_Reflect()
                reflection_v = Vector_Reflect(light_v, normal);

                float r_dot_v = Vector_Dot(reflection_v, eye);
                if (r_dot_v > 0.0f)
                    total_intensity += light.intensity * powf(r_dot_v / (Vector_Length(reflection_v) * Vector_Length(eye)), specular);
            }

        }
    }

    return total_intensity;
}

// Computes the intersection of the ray with every sphere 
// and returns the color of the sphere at the nearest intersection inside the requested range of t
TVECTOR trace_ray(const TVECTOR ray_origin, const TVECTOR ray_direction, const float near_plane, const float far_plane, const int32_t recursion_depth = 0)
{
    Sphere closest_sphere;
    float closest_t = closest_intersection(ray_origin, ray_direction, near_plane, far_plane, closest_sphere);

    if (closest_t == far_plane) // if there was no sphere intersection
        return BACKGROUND_COLOR;

    TVECTOR point = Vector_Add(ray_origin, Vector_Scalar_Multiply(ray_direction, closest_t));// Compute Intersection (O + t*D)
    TVECTOR normal = Vector_Normalize(Vector_Sub(point, closest_sphere.center)); // Compute Sphere Normal
    TVECTOR eye = Vector_Negate(ray_direction);

    float lighting = compute_lighting(point, normal, eye, closest_sphere.specular_exponent);
    TVECTOR local_color = Vector_Scalar_Multiply(closest_sphere.color, 1.3f*lighting);

    // If we hit the recursion limit or the object is not reflective, we’ re done
    float reflective = closest_sphere.reflective;
    if (recursion_depth <= 0 || reflective <= 0.0f)
        return Vector_Saturate(local_color);

    // Compute the reflected color
    TVECTOR reflection_v = Vector_Reflect(eye, normal);
    TVECTOR reflected_color = trace_ray(point, reflection_v, 0.05f, INFINITY, recursion_depth - 1);

    TVECTOR out_color = Vector_Add(
        Vector_Scalar_Multiply(local_color, (1.0f - reflective)),
        Vector_Scalar_Multiply(reflected_color, reflective));

    return Vector_Saturate(out_color);
}

#ifdef USE_MULTITHREAD
void PaintPixels(ThreadData* threadData, const int32_t width, const int32_t height)
{
    if (threadData == nullptr) {
        return;
    }

    const int32_t x1{ threadData->min[0] };
    const int32_t y1{ threadData->min[1] };
    const int32_t x2{ threadData->max[0] };
    const int32_t y2{ threadData->max[1] };
    
    // for every pixel on the canvas
    for (int32_t y = y1; y < y2; ++y)
    {
        for (int32_t x = x1; x < x2; ++x)
        {
            // convert 2d coordinate to 3d space
            TVECTOR ray = canvas_to_viewport(x, y);
            TVECTOR color = trace_ray(scene.cam_pos, ray, 1.0f, INFINITY, 3);

            threadData->mutex->lock();
            int32_t p_index = (y + height) * canvas_width + (x + width);
            pixels[p_index].red = static_cast<uint8_t>(color.x * 255);
            pixels[p_index].green = static_cast<uint8_t>(color.y * 255);
            pixels[p_index].blue = static_cast<uint8_t>(color.z * 255);
            threadData->mutex->unlock();
        }
    }
}
#endif

int main()
{
    const int32_t canvas_size = canvas_width * canvas_height;
    pixels = std::vector<Pixel>(canvas_size + 1);

    // define the scene
    {
        Sphere s1, s2, s3, s4;
        // red sphere
        s1.color = { 1.0f, 0.0f, 0.0f, 1.0f };
        s1.center = { 0.0f, -1.0f, 3.0f, 1.0f };
        s1.radius = 1.0f;
        s1.specular_exponent = 500.0f; // shiny
        s1.reflective = 0.2f;
        // green sphere
        s2.color = { 0.0f, 1.0f, 0.0f, 1.0f };
        s2.center = { -2.0f, 0.0f, 4.0f, 1.0f};
        s2.radius = 1.0f;
        s2.specular_exponent = 10.0f; // not so shiny
        s2.reflective = 0.4f;
        // blue sphere
        s3.color = { 0.0f, 0.0f, 1.0f, 1.0f };
        s3.center = { 2.0f, 0.0f, 4.0f, 1.0f};
        s3.radius = 1.0f;
        s3.specular_exponent = 500.0f; // shiny
        s3.reflective = 0.3f;
        // yellow sphere
        s4.color = { 1.0f, 1.0f, 0.0f, 1.0f };
        s4.center = { 0.0f, -5001.0f, 0.0f, 1.0f };
        s4.radius = 5000.0f;
        s4.specular_exponent = 1000.0f; // very shiny
        s4.reflective = 0.5f;

        scene.spheres.push_back(s1);
        scene.spheres.push_back(s2);
        scene.spheres.push_back(s3);
        scene.spheres.push_back(s4);

        Light amb_l, point_l, dir_l;
        // Ambient light
        amb_l.type = ELightType::Ambient;
        amb_l.intensity = 0.2f;
        // Point light
        point_l.type = ELightType::Point;
        point_l.intensity = 0.6f;
        point_l.postion = { 1.0f, 0.0f, 1.0f, 1.0f };
        // Directional light
        dir_l.type = ELightType::Directional;
        dir_l.intensity = 0.2f;
        dir_l.direction = { 1.0f, 4.0f, 4.0f, 0.0f };

        scene.lights.push_back(amb_l);
        scene.lights.push_back(point_l);
        scene.lights.push_back(dir_l);
    }
    
    const int32_t width = canvas_width >> 1;
    const int32_t height = canvas_height >> 1;

#ifdef USE_MULTITHREAD
    ThreadData* perThreadData = new ThreadData[4];

    std::vector<std::thread> threads; // 4 threads
    std::mutex mutex;
    std::condition_variable conditional;

    int32_t threadIndex{ 0 };
    for (int32_t i = 0; i < 2; ++i) {
        for (int32_t j = 0; j < 2; ++j)
        {
            int32_t x1 = (i * width) - width;
            int32_t y1 = (j * height) - height;
            int32_t x2 = x1 + width;
            int32_t y2 = y1 + height;

            perThreadData[threadIndex].mutex = &mutex;
            perThreadData[threadIndex].cnd = &conditional;
            perThreadData[threadIndex].min[0] = x1; perThreadData[threadIndex].min[1] = y1;
            perThreadData[threadIndex].max[0] = x2; perThreadData[threadIndex].max[1] = y2;

            threads.push_back(std::thread(PaintPixels, &perThreadData[threadIndex], width, height));
            ++threadIndex;
        }
    }

    for (auto& thread : threads) {
        thread.join();
    }

    delete[] perThreadData;
#else
    // for every pixel on the canvas
    for (int32_t y = -height; y < height; ++y)
    {
        for (int32_t x = -width; x < width; ++x)
        {
            // convert 2d coordinate to 3d space
            TVECTOR ray = canvas_to_viewport(x, y);
            TVECTOR color = trace_ray(scene.cam_pos, ray, 1.0f, INFINITY, 3);
    
            int p_index = (y + height) * canvas_width + (x + width);
            pixels[p_index].red =  static_cast<uint8_t>(color.x * 255);
            pixels[p_index].green =  static_cast<uint8_t>(color.y * 255);
            pixels[p_index].blue =  static_cast<uint8_t>(color.z * 255);
        }
    }
#endif

    // PutPixel
    std::ofstream file("RayTrace.bmp", std::ios::trunc | std::ios::out | std::ios::binary);

    assert(file.is_open());

    if (file.is_open())
    {
        bmpHeader.write_to_file(file);
        bmpInfoHeader.write_to_file(file);

        for (int32_t i = 0; i < canvas_size; ++i) {
            pixels[i].write_to_file(file);
        }
    }
    file.close();

    return 0;
}
