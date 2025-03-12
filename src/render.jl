using LinearAlgebra
using Images
using StaticArrays

# Abstract types for our object hierarchy
abstract type Shape end
abstract type CSGOperation <: Shape end

# Basic geometric primitives
struct Sphere <: Shape
    # Internal
    center::SVector{3, Float64}
    rotation::SVector{3, Float64}
    color::RGBA{N0f8}

    # User-facing
    radius::Float64
end

struct Box <: Shape
    # Internal
    center::SVector{3, Float64}
    rotation::SVector{3, Float64}
    color::RGBA{N0f8}

    # User-facing
    dimensions::SVector{3, Float64}
end

struct Cylinder <: Shape
    # Internal
    center::SVector{3, Float64}
    rotation::SVector{3, Float64}
    color::RGBA{N0f8}

    # User-facing
    radius_top::Float64
    radius_bottom::Float64
    height::Float64
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

# Create basic transformation matrices
function translation(tx::Float64, ty::Float64, tz::Float64)
    SMatrix{4,4}([
        1.0 0.0 0.0 tx;
        0.0 1.0 0.0 ty;
        0.0 0.0 1.0 tz;
        0.0 0.0 0.0 1.0
    ])
end

function scaling(sx::Float64, sy::Float64, sz::Float64)
    SMatrix{4,4}([
        sx  0.0 0.0 0.0;
        0.0 sy  0.0 0.0;
        0.0 0.0 sz  0.0;
        0.0 0.0 0.0 1.0
    ])
end

function rotation_x(angle::Float64)
    c, s = cos(angle), sin(angle)
    SMatrix{4,4}([
        1.0 0.0 0.0 0.0;
        0.0 c   -s  0.0;
        0.0 s   c   0.0;
        0.0 0.0 0.0 1.0
    ])
end

function rotation_y(angle::Float64)
    c, s = cos(angle), sin(angle)
    SMatrix{4,4}([
        c   0.0 s   0.0;
        0.0 1.0 0.0 0.0;
        -s  0.0 c   0.0;
        0.0 0.0 0.0 1.0
    ])
end

function rotation_z(angle::Float64)
    c, s = cos(angle), sin(angle)
    SMatrix{4,4}([
        c   -s  0.0 0.0;
        s   c   0.0 0.0;
        0.0 0.0 1.0 0.0;
        0.0 0.0 0.0 1.0
    ])
end

function rotation_matrix(rotation::SVector{3, Float64})
    return rotation_x(rotation[1]) * rotation_y(rotation[2]) * rotation_z(rotation[3])
end

# Ray structure for ray casting
struct Ray
    origin::SVector{3, Float64}
    direction::SVector{3, Float64}
end

function intersect(ray::Ray, union::CSGUnion)
    t = Inf
    color = RGBA{N0f8}(0,0,0,1)
    # Find closest intersection and get its color
    for shape in union.shapes
        curr_t, curr_color = intersect(ray, shape)
        if curr_t < t
            t = curr_t
            color = curr_color
        end
    end
    return (t, color)
end

function intersect(ray::Ray, intersection::Intersection)
    t = Inf

    color = RGBA{N0f8}(1,1,1,1)
    for shape in intersection.shapes
        curr_t, curr_color = intersect(ray, shape)
        if curr_t != Inf
            t = curr_t
            color = RGBA{N0f8}(
                color.r * curr_color.r,
                color.g * curr_color.g,
                color.b * curr_color.b,
                color.a * curr_color.a
            )
        end
    end
    if t == Inf
        return (Inf, RGBA{N0f8}(0,0,0,1))
    end
    return (t, color)
end

function intersect(ray::Ray, diff::Difference)
    t_a, color_a = intersect(ray, diff.a)
    t_b, _ = intersect(ray, diff.b)

    if t_a == Inf
        # Missed a, so we hit the background
        return (Inf, RGBA{N0f8}(0,0,0,1))
    elseif t_b == Inf || t_a < t_b
        # Hit a
        return (t_a, color_a)
    else
        # "Hit" b, so we hit the background
        return (Inf, RGBA{N0f8}(0,0,0,1))
    end
end

# Intersection testing functions
function intersect(ray::Ray, sphere::Sphere)
    # Transform ray to sphere's local coordinate system
    rot_matrix = rotation_matrix(-sphere.rotation) # Inverse rotation
    # Convert direction to homogeneous coordinates (w=0 for vectors)
    local_dir = SVector{3}((rot_matrix * SVector{4}(ray.direction[1], ray.direction[2], ray.direction[3], 0.0))[1:3])
    # Convert origin to homogeneous coordinates (w=1 for points)
    local_origin = SVector{3}((rot_matrix * SVector{4}((ray.origin - sphere.center)..., 1.0))[1:3])

    # Standard sphere intersection in local coordinates
    a = dot(local_dir, local_dir)
    b = 2.0 * dot(local_origin, local_dir)
    c = dot(local_origin, local_origin) - sphere.radius^2
    discriminant = b^2 - 4*a*c

    if discriminant < 0
        return (Inf, RGBA{N0f8}(0,0,0,1))
    else
        t = (-b - sqrt(discriminant)) / (2.0 * a)
        return (t > 0 ? t : Inf, sphere.color)
    end
end

function intersect(ray::Ray, box::Box)
    # Transform ray to box's local coordinate system
    rot_matrix = rotation_matrix(-box.rotation) # Inverse rotation
    # Convert direction to homogeneous coordinates (w=0 for vectors)
    local_dir = SVector{3}((rot_matrix * SVector{4}(ray.direction[1], ray.direction[2], ray.direction[3], 0.0))[1:3])
    # Convert origin to homogeneous coordinates (w=1 for points)
    local_origin = SVector{3}((rot_matrix * SVector{4}((ray.origin - box.center)..., 1.0))[1:3])

    # Compute intersection with axis-aligned box in local space
    min_p = -box.dimensions/2
    max_p = box.dimensions/2

    tx1 = (min_p[1] - local_origin[1]) / local_dir[1]
    tx2 = (max_p[1] - local_origin[1]) / local_dir[1]
    tmin = min(tx1, tx2)
    tmax = max(tx1, tx2)

    ty1 = (min_p[2] - local_origin[2]) / local_dir[2]
    ty2 = (max_p[2] - local_origin[2]) / local_dir[2]
    tmin = max(tmin, min(ty1, ty2))
    tmax = min(tmax, max(ty1, ty2))

    tz1 = (min_p[3] - local_origin[3]) / local_dir[3]
    tz2 = (max_p[3] - local_origin[3]) / local_dir[3]
    tmin = max(tmin, min(tz1, tz2))
    tmax = min(tmax, max(tz1, tz2))

    if tmax >= tmin && tmax > 0
        return (max(0.0, tmin), box.color)
    else
        return (Inf, RGBA{N0f8}(0,0,0,1))
    end
end

function intersect(ray::Ray, cylinder::Cylinder)
    # Get cylinder/cone parameters
    r1 = cylinder.radius_bottom
    r2 = cylinder.radius_top
    height = cylinder.height

    # Transform ray to cylinder's local coordinate system
    rot_matrix = rotation_matrix(-cylinder.rotation) # Inverse rotation
    # Convert direction to homogeneous coordinates (w=0 for vectors)
    local_dir = SVector{3}((rot_matrix * SVector{4}(ray.direction[1], ray.direction[2], ray.direction[3], 0.0))[1:3])
    # Convert origin to homogeneous coordinates (w=1 for points)
    local_origin = SVector{3}((rot_matrix * SVector{4}((ray.origin - cylinder.center)..., 1.0))[1:3])

    # Calculate cone parameters
    tan_alpha = (r1 - r2) / height  # Slope of cone wall

    # Quadratic equation coefficients for cone
    a = local_dir[1]^2 + local_dir[2]^2 - 
        (tan_alpha * local_dir[3])^2
    b = 2.0 * (local_origin[1]*local_dir[1] + local_origin[2]*local_dir[2] - 
        tan_alpha^2 * local_origin[3]*local_dir[3] - 
        r1 * tan_alpha * local_dir[3])
    c = local_origin[1]^2 + local_origin[2]^2 - 
        (tan_alpha * local_origin[3] + r1)^2

    # Check if ray is parallel to cone surface
    if abs(a) < 1e-6
        return (Inf, RGBA{N0f8}(0,0,0,1))
    end

    # Solve quadratic equation
    discriminant = b^2 - 4*a*c
    if discriminant < 0
        return (Inf, RGBA{N0f8}(0,0,0,1))
    end

    # Find intersection points with infinite cone
    t1 = (-b - sqrt(discriminant)) / (2.0 * a)
    t2 = (-b + sqrt(discriminant)) / (2.0 * a)

    # Sort intersection points
    if t1 > t2
        t1, t2 = t2, t1
    end

    # Find z coordinates of intersection points in local space
    z1 = local_origin[3] + t1*local_dir[3]
    z2 = local_origin[3] + t2*local_dir[3]

    # Check if intersection points are within cylinder height
    half_height = height/2

    # Initialize result as no intersection
    result = Inf

    # Check both intersection points
    for (t, z) in [(t1, z1), (t2, z2)]
        if t > 0 && -half_height <= z <= half_height
            # Calculate radius at intersection point
            r_at_z = r1 + (r2 - r1) * (z + half_height) / height

            # Calculate actual radius at intersection point in local space
            p = local_origin + t * local_dir
            r_intersect = sqrt(p[1]^2 + p[2]^2)

            # Check if intersection point lies on cone surface
            if abs(r_intersect - r_at_z) < 1e-6
                result = t
                break
            end
        end
    end

    if result == Inf
        return (Inf, RGBA{N0f8}(0,0,0,1))
    else
        return (result, cylinder.color)
    end
end

function ast_to_csg(node::NumberNode)
    return node.value
end

function ast_to_csg(node::VectorNode)
    return SVector([ast_to_csg(child) for child in node.values]...)
end

function apply_translation(sphere::Sphere, translation_vector::SVector{3, Float64})
    return Sphere(sphere.center + translation_vector, sphere.rotation, sphere.color, sphere.radius)
end

function apply_translation(box::Box, translation_vector::SVector{3, Float64})
    return Box(box.center + translation_vector, box.rotation, box.color, box.dimensions)
end

function apply_translation(cylinder::Cylinder, translation_vector::SVector{3, Float64})
    return Cylinder(cylinder.center + translation_vector, cylinder.rotation, cylinder.color, cylinder.radius_top, cylinder.radius_bottom, cylinder.height)
end

function apply_rotation(shape::Sphere, rotation_vector::SVector{3, Float64})
    return Sphere(shape.center, rotation_vector, shape.color, shape.radius)
end

function apply_rotation(shape::Box, rotation_vector::SVector{3, Float64})
    return Box(shape.center, rotation_vector, shape.color, shape.dimensions)
end

function apply_rotation(shape::Cylinder, rotation_vector::SVector{3, Float64})
    return Cylinder(shape.center, rotation_vector, shape.color, shape.radius_top, shape.radius_bottom, shape.height)
end

function apply_color(shape::Sphere, color::RGBA{N0f8})
    return Sphere(shape.center, shape.rotation, color, shape.radius)
end

function apply_color(shape::Box, color::RGBA{N0f8})
    return Box(shape.center, shape.rotation, color, shape.dimensions)
end

function apply_color(shape::Cylinder, color::RGBA{N0f8})
    return Cylinder(shape.center, shape.rotation, color, shape.radius_top, shape.radius_bottom, shape.height)
end

function ast_to_csg(node::CallNode)
    if node.name == "union"
        if length(node.children) == 1
            return ast_to_csg(node.children[1])
        else
            return CSGUnion([ast_to_csg(child) for child in node.children])
        end
    elseif node.name == "intersection"
        @assert length(node.children) > 1
        return Intersection([ast_to_csg(child) for child in node.children])
    elseif node.name == "difference"
        @assert length(node.children) == 2
        return Difference(ast_to_csg(node.children[1]), ast_to_csg(node.children[2]))
    elseif node.name == "sphere"
        @assert length(node.children) == 0
        diameter = ast_to_csg(node.arguments[findfirst(arg->arg[1]=="d", node.arguments)][2])
        return Sphere(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), RGBA{N0f8}(1,1,1,1), diameter / 2)
    elseif node.name == "cube"
        @assert length(node.children) == 0
        size = ast_to_csg(node.arguments[1][2]::VectorNode)
        return Box(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), RGBA{N0f8}(1,1,1,1), size)
    elseif node.name == "cylinder"
        @assert length(node.children) == 0
        height = ast_to_csg(node.arguments[1][2])
        radius_top = ast_to_csg(node.arguments[2][2])
        radius_bottom = ast_to_csg(node.arguments[3][2])
        return Cylinder(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), RGBA{N0f8}(1,1,1,1), radius_top, radius_bottom, height)
    elseif node.name == "translate"
        # We immediately apply the translation to the children
        @assert length(node.arguments) == 1
        translation_vector = ast_to_csg(node.arguments[1][2]::VectorNode)
        children_transformed = [apply_translation(ast_to_csg(child), translation_vector) for child in node.children]
        if length(children_transformed) == 1
            return children_transformed[1]
        else
            return CSGUnion(children_transformed)
        end
    elseif node.name == "rotate"
        @assert length(node.arguments) == 1
        rotation_vector = ast_to_csg(node.arguments[1][2]::VectorNode)
        children_transformed = [apply_rotation(ast_to_csg(child), rotation_vector) for child in node.children]
        if length(children_transformed) == 1
            return children_transformed[1]
        else
            return CSGUnion(children_transformed)
        end
    elseif node.name == "color"
        @assert length(node.children) > 0
        if length(node.arguments) == 1
            color = ast_to_csg(node.arguments[1][2])
        elseif length(node.arguments) == 3
            color = RGBA{N0f8}(ast_to_csg(node.arguments[1][2]), ast_to_csg(node.arguments[2][2]), ast_to_csg(node.arguments[3][2]), 1.0)
        else
            error("Unknown color node: $node")
        end
        children_transformed = [apply_color(ast_to_csg(child), color) for child in node.children]
        if length(children_transformed) == 1
            return children_transformed[1]
        else
            return CSGUnion(children_transformed)
        end
    else
        error("Unknown call node name: $(node.name)")
    end
end

function ast_to_csg(node::ModuleNode)
    # FIXME: Register the module name so later CallNodes can reference it
    if length(node.children) == 1
        return ast_to_csg(node.children[1])
    else
        return CSGUnion([ast_to_csg(child) for child in node.children])
    end
end

function ast_to_csg(asts::Vector{ASTNode})
    if length(asts) == 1
        return ast_to_csg(asts[1])
    else
        return CSGUnion([ast_to_csg(ast) for ast in asts])
    end
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
    color_type::Type{<:Colorant}=RGBA{N0f8},
    background_color::Colorant=RGBA{N0f8}(0.0, 0.0, 0.0, 1.0)
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
    Threads.@threads for i in 1:width
        for j in 1:height
            # Calculate ray direction
            x = (2 * (i - 0.5) / width - 1) * tan_half_fov * aspect_ratio
            y = (1 - 2 * (j - 0.5) / height) * tan_half_fov

            direction = normalize(x*u + y*v - w)
            ray = Ray(camera.position, direction)

            # Perform intersection test
            t, color = intersect(ray, shape)

            # Set pixel color based on intersection
            if t != Inf
                # Simple shading based on distance
                intensity = clamp(1.0 - t/10.0, 0.0, 1.0)
                @inbounds img[j, i] = RGBA{N0f8}(color.r * intensity, color.g * intensity, color.b * intensity, color.alpha)
            else
                @inbounds img[j, i] = background_color
            end
        end
    end

    return img
end