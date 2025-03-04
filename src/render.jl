using LinearAlgebra
using Images
using StaticArrays

# Abstract types for our object hierarchy
abstract type Shape end
abstract type CSGOperation end

# Basic geometric primitives
struct Sphere <: Shape
    center::SVector{3, Float64}
    radius::Float64
end

struct Box <: Shape
    center::SVector{3, Float64}
    dimensions::SVector{3, Float64}
end

# CSG Operations
struct CSGUnion <: CSGOperation
    shapes::Vector{Shape}
end

struct Intersection <: CSGOperation
    shapes::Vector{Shape}
end

struct Difference <: CSGOperation
    a::Shape
    b::Shape
end

# Transform structure
struct Transform
    matrix::SMatrix{4, 4, Float64}
end

# Create basic transformation matrices
function translation(tx::Float64, ty::Float64, tz::Float64)
    Transform(SMatrix{4,4}([
        1.0 0.0 0.0 tx;
        0.0 1.0 0.0 ty;
        0.0 0.0 1.0 tz;
        0.0 0.0 0.0 1.0
    ]))
end

function scaling(sx::Float64, sy::Float64, sz::Float64)
    Transform(SMatrix{4,4}([
        sx  0.0 0.0 0.0;
        0.0 sy  0.0 0.0;
        0.0 0.0 sz  0.0;
        0.0 0.0 0.0 1.0
    ]))
end

function rotation_x(angle::Float64)
    c, s = cos(angle), sin(angle)
    Transform(SMatrix{4,4}([
        1.0 0.0 0.0 0.0;
        0.0 c   -s  0.0;
        0.0 s   c   0.0;
        0.0 0.0 0.0 1.0
    ]))
end

function rotation_y(angle::Float64)
    c, s = cos(angle), sin(angle)
    Transform(SMatrix{4,4}([
        c   0.0 s   0.0;
        0.0 1.0 0.0 0.0;
        -s  0.0 c   0.0;
        0.0 0.0 0.0 1.0
    ]))
end

function rotation_z(angle::Float64)
    c, s = cos(angle), sin(angle)
    Transform(SMatrix{4,4}([
        c   -s  0.0 0.0;
        s   c   0.0 0.0;
        0.0 0.0 1.0 0.0;
        0.0 0.0 0.0 1.0
    ]))
end

# Ray structure for ray casting
struct Ray
    origin::SVector{3, Float64}
    direction::SVector{3, Float64}
end

# Intersection testing functions
function intersect_sphere(ray::Ray, sphere::Sphere)
    oc = ray.origin - sphere.center
    a = dot(ray.direction, ray.direction)
    b = 2.0 * dot(oc, ray.direction)
    c = dot(oc, oc) - sphere.radius^2
    discriminant = b^2 - 4*a*c
    
    if discriminant < 0
        return Inf
    else
        t = (-b - sqrt(discriminant)) / (2.0 * a)
        return t > 0 ? t : Inf
    end
end

function intersect_box(ray::Ray, box::Box)
    # Compute intersection with axis-aligned box
    min_p = box.center - box.dimensions/2
    max_p = box.center + box.dimensions/2
    
    tx1 = (min_p[1] - ray.origin[1]) / ray.direction[1]
    tx2 = (max_p[1] - ray.origin[1]) / ray.direction[1]
    tmin = min(tx1, tx2)
    tmax = max(tx1, tx2)
    
    ty1 = (min_p[2] - ray.origin[2]) / ray.direction[2]
    ty2 = (max_p[2] - ray.origin[2]) / ray.direction[2]
    tmin = max(tmin, min(ty1, ty2))
    tmax = min(tmax, max(ty1, ty2))
    
    tz1 = (min_p[3] - ray.origin[3]) / ray.direction[3]
    tz2 = (max_p[3] - ray.origin[3]) / ray.direction[3]
    tmin = max(tmin, min(tz1, tz2))
    tmax = min(tmax, max(tz1, tz2))
    
    return tmax >= tmin && tmax > 0 ? max(0.0, tmin) : Inf
end

# CSG operation intersection functions
function intersect(ray::Ray, union::CSGUnion)
    t = Inf
    for shape in union.shapes
        t = min(t, intersect(ray, shape))
    end
    return t
end

function intersect(ray::Ray, intersection::Intersection)
    t = 0.0
    for shape in intersection.shapes
        curr_t = intersect(ray, shape)
        if curr_t == Inf
            return Inf
        end
        t = max(t, curr_t)
    end
    return t
end

function intersect(ray::Ray, diff::Difference)
    t_a = intersect(ray, diff.a)
    t_b = intersect(ray, diff.b)
    
    if t_a == Inf
        return Inf
    elseif t_b == Inf || t_a < t_b
        return t_a
    else
        return Inf
    end
end

# Dispatch for different shape types
function intersect(ray::Ray, sphere::Sphere)
    return intersect_sphere(ray, sphere)
end

function intersect(ray::Ray, box::Box)
    return intersect_box(ray, box)
end

# Camera structure
struct Camera
    position::SVector{3, Float64}
    direction::SVector{3, Float64}
    up::SVector{3, Float64}
    fov::Float64
end

# Main rasterizer function
function render(
    shape::Union{Shape, CSGOperation},
    camera::Camera,
    width::Int,
    height::Int,
    color_type::Type{<:Colorant}=RGBA{N0f8}
)
    # Create image buffer
    img = zeros(color_type, height, width)
    
    # Camera parameters
    aspect_ratio = width / height
    tan_half_fov = tan(camera.fov / 2)
    
    # Calculate camera basis vectors
    w = normalize(-camera.direction)
    u = normalize(cross(camera.up, w))
    v = cross(w, u)
    
    # For each pixel
    Threads.@threads for j in 1:height
        for i in 1:width
            # Calculate ray direction
            x = (2 * (i - 0.5) / width - 1) * tan_half_fov * aspect_ratio
            y = (1 - 2 * (j - 0.5) / height) * tan_half_fov
            
            direction = normalize(x*u + y*v - w)
            ray = Ray(camera.position, direction)
            
            # Perform intersection test
            t = intersect(ray, shape)
            
            # Set pixel color based on intersection
            if t != Inf
                # Simple shading based on distance
                intensity = clamp(1.0 - t/10.0, 0.0, 1.0)
                img[j, i] = color_type(intensity, intensity, intensity, 1.0)
            else
                img[j, i] = color_type(0.0, 0.0, 0.0, 1.0)
            end
        end
    end
    
    return img
end

# Example usage
function example()
    # Create some basic shapes
    sphere1 = Sphere(SVector(0.0, 0.0, -5.0), 1.0)
    sphere2 = Sphere(SVector(1.0, 0.0, -5.0), 1.0)
    box = Box(SVector(-1.0, 0.0, -5.0), SVector(1.0, 1.0, 1.0))
    
    # Create a CSG union of shapes
    scene = CSGUnion([sphere1, sphere2, box])
    
    # Set up camera
    camera = Camera(
        SVector(0.0, 0.0, 0.0),  # position
        SVector(0.0, 0.0, -1.0), # direction
        SVector(0.0, 1.0, 0.0),  # up
        π/3                       # fov (60 degrees)
    )
    
    # Render the scene
    img = render(scene, camera, 800, 600, RGBA{N0f8})
    
    # Save the image
    save("output.png", img)
end
