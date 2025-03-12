using SCAD3D
using Images
using StaticArrays

# Create some basic shapes
# FIXME: Support "translate() ...;
input = """
union() {
    translate([0, 0, -5]) {
        rotate([0, 0, 0]) {
            color(r = 1, g = 0, b = 1) {
                sphere(d = 2);
            }
        }
    }
    translate([1, 0, -5]) {
        color(r = 0, g = 1, b = 0) {
            sphere(d = 2);
        }
    }
    translate([-1, 0, -5]) {
        rotate([0, 1, 1]) {
            color(r = 1, g = 0, b = 0) {
                cube([1, 1, 1]);
            }
        }
    }
    translate([-2, -1, -6]) {
        rotate([1, 0.5, 0.5]) {
            color(r = 0, g = 0, b = 1) {
                cylinder(h = 5, r1 = 2, r2 = 2);
            }
        }
    }
}
"""

# Parse the input
ast = SCAD3D.parse_scad(input)

# Print the AST
for node in ast
    SCAD3D.print_ast(node)
end

# Show the SCAD code
println(SCAD3D.write_scad(ast))

# Convert the AST to a CSG tree
scene = SCAD3D.ast_to_csg(ast)

# Create a camera
camera = SCAD3D.Camera(
    SVector(0.0, 0.0, 0.0),  # position
    SVector(0.0, 0.0, -1.0), # direction
    SVector(0.0, 1.0, 0.0),  # up
    π/3                       # fov (60 degrees)
)
# Render the scene
img = @time SCAD3D.render(scene, camera, 600, 400, RGBA{N0f8}, RGBA{N0f8}(0.2, 0.2, 0.2, 1.0))

SCAD3D.save("output.png", img)