using SCAD3D

# Example OpenSCAD code
input = """
module hole(diameter = 10, height = 5) {
    cylinder(d = diameter, h = height);
}

difference() {
    cube([20, 20, 10]);
    translate([10, 10, 0]) {
        hole(diameter = 8);
    }
}
"""

# Parse the input
ast = SCAD3D.parse_scad(input)

# Print the AST
for node in ast
    SCAD3D.print_ast(node)
end

println(SCAD3D.write_scad(ast))
