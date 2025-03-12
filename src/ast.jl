# Token types for lexical analysis
@enum TokenType begin
    NUMBER
    IDENTIFIER
    STRING
    OPERATOR
    LPAREN
    RPAREN
    LBRACE
    RBRACE
    LBRACKET
    RBRACKET
    SEMICOLON
    COMMA
    EQUALS
    KEYWORD
    EOF
end

# Represents a token in the OpenSCAD code
struct Token
    type::TokenType
    value::String
    line::Int
    column::Int
end

# Abstract syntax tree node types
abstract type ASTNode end

struct ModuleNode <: ASTNode
    name::String
    parameters::Vector{Pair{String, Union{Nothing, Any}}}
    children::Vector{ASTNode}
end

struct CallNode <: ASTNode
    name::String
    arguments::Vector{Pair{String, Any}}
    children::Vector{ASTNode}
end

struct NumberNode <: ASTNode
    value::Float64
end

struct VectorNode <: ASTNode
    values::Vector{ASTNode}
end

struct VariableNode <: ASTNode
    name::String
end

# Lexer implementation
mutable struct Lexer
    input::String
    position::Int
    line::Int
    column::Int
    current_char::Union{Nothing, Char}
end

function Lexer(input::String)
    lexer = Lexer(input, 1, 1, 1, nothing)
    if length(input) > 0
        lexer.current_char = input[1]
    end
    return lexer
end

function advance!(lexer::Lexer)
    lexer.position += 1
    if lexer.position <= length(lexer.input)
        lexer.current_char = lexer.input[lexer.position]
        if lexer.current_char == '\n'
            lexer.line += 1
            lexer.column = 1
        else
            lexer.column += 1
        end
    else
        lexer.current_char = nothing
    end
end

function skip_whitespace!(lexer::Lexer)
    while lexer.current_char !== nothing && isspace(lexer.current_char)
        advance!(lexer)
    end
end

function skip_comment!(lexer::Lexer)
    if lexer.current_char == '/'
        advance!(lexer)
        if lexer.current_char == '/'
            while lexer.current_char !== nothing && lexer.current_char != '\n'
                advance!(lexer)
            end
        elseif lexer.current_char == '*'
            nested = 1
            advance!(lexer)
            while nested > 0 && lexer.current_char !== nothing
                if lexer.current_char == '/' && peek_next(lexer) == '*'
                    nested += 1
                    advance!(lexer)
                elseif lexer.current_char == '*' && peek_next(lexer) == '/'
                    nested -= 1
                    advance!(lexer)
                end
                advance!(lexer)
            end
        end
    end
end

function peek_next(lexer::Lexer)
    if lexer.position + 1 <= length(lexer.input)
        return lexer.input[lexer.position + 1]
    end
    return nothing
end

function read_number!(lexer::Lexer)
    result = ""
    while lexer.current_char !== nothing && (isdigit(lexer.current_char) || lexer.current_char == '.' || lexer.current_char == '-')
        result *= lexer.current_char
        advance!(lexer)
    end
    return Token(NUMBER, result, lexer.line, lexer.column)
end

function read_identifier!(lexer::Lexer)
    result = ""
    while lexer.current_char !== nothing && (isletter(lexer.current_char) || isdigit(lexer.current_char) || lexer.current_char == '_')
        result *= lexer.current_char
        advance!(lexer)
    end

    # Check if it's a keyword
    keywords = Set(["module", "function"])
    #, "if", "else", "for", "intersection", "union", "difference"])
    if result in keywords
        return Token(KEYWORD, result, lexer.line, lexer.column)
    end

    return Token(IDENTIFIER, result, lexer.line, lexer.column)
end

function read_string!(lexer::Lexer)
    result = ""
    advance!(lexer)  # Skip opening quote
    while lexer.current_char !== nothing && lexer.current_char != '"'
        if lexer.current_char == '\\'
            advance!(lexer)
            if lexer.current_char !== nothing
                result *= lexer.current_char
            end
        else
            result *= lexer.current_char
        end
        advance!(lexer)
    end
    advance!(lexer)  # Skip closing quote
    return Token(STRING, result, lexer.line, lexer.column)
end

function get_next_token(lexer::Lexer)
    while lexer.current_char !== nothing
        if isspace(lexer.current_char)
            skip_whitespace!(lexer)
            continue
        end
        
        if lexer.current_char == '/'
            skip_comment!(lexer)
            continue
        end
        
        if isdigit(lexer.current_char) || (lexer.current_char == '-' && peek_next(lexer) !== nothing && isdigit(peek_next(lexer)))
            return read_number!(lexer)
        end
        
        if isletter(lexer.current_char)
            return read_identifier!(lexer)
        end
        
        if lexer.current_char == '"'
            return read_string!(lexer)
        end
        
        # Single-character tokens
        char_tokens = Dict(
            '(' => LPAREN,
            ')' => RPAREN,
            '{' => LBRACE,
            '}' => RBRACE,
            '[' => LBRACKET,
            ']' => RBRACKET,
            ';' => SEMICOLON,
            ',' => COMMA,
            '=' => EQUALS
        )
        
        if haskey(char_tokens, lexer.current_char)
            token = Token(char_tokens[lexer.current_char], string(lexer.current_char), lexer.line, lexer.column)
            advance!(lexer)
            return token
        end
        
        # Operators
        if lexer.current_char in ['+', '-', '*', '/', '%']
            token = Token(OPERATOR, string(lexer.current_char), lexer.line, lexer.column)
            advance!(lexer)
            return token
        end
        
        error("Invalid character: $(lexer.current_char) at line $(lexer.line), column $(lexer.column)")
    end
    
    return Token(EOF, "", lexer.line, lexer.column)
end

# Parser implementation
mutable struct Parser
    lexer::Lexer
    current_token::Token
end

function Parser(input::String)
    lexer = Lexer(input)
    parser = Parser(lexer, get_next_token(lexer))
    return parser
end

function eat!(parser::Parser, token_type::TokenType)
    if parser.current_token.type == token_type
        parser.current_token = get_next_token(parser.lexer)
    else
        error("Expected token type $token_type, got $(parser.current_token.type)")
    end
end

function parse_number(parser::Parser)
    token = parser.current_token
    eat!(parser, NUMBER)
    return NumberNode(parse(Float64, token.value))
end

function parse_vector(parser::Parser)
    values = ASTNode[]
    eat!(parser, LBRACKET)
    
    if parser.current_token.type != RBRACKET
        push!(values, parse_expression(parser))
        while parser.current_token.type == COMMA
            eat!(parser, COMMA)
            push!(values, parse_expression(parser))
        end
    end
    
    eat!(parser, RBRACKET)
    return VectorNode(values)
end
function parse_variable(parser::Parser)
    name = parser.current_token.value
    eat!(parser, IDENTIFIER)
    return VariableNode(name)
end

function parse_expression(parser::Parser)
    if parser.current_token.type == NUMBER
        return parse_number(parser)
    elseif parser.current_token.type == LBRACKET
        return parse_vector(parser)
    elseif parser.current_token.type == IDENTIFIER
        return parse_variable(parser)
    else
        error("Unexpected token in expression: $(parser.current_token.type)\n$(parser.current_token)")
    end
end

function parse_module(parser::Parser)
    eat!(parser, KEYWORD)  # Consume 'module'
    name = parser.current_token.value
    eat!(parser, IDENTIFIER)

    parameters = Pair{String, Union{Nothing, Any}}[]
    eat!(parser, LPAREN)

    if parser.current_token.type != RPAREN
        param_name = parser.current_token.value
        eat!(parser, IDENTIFIER)

        if parser.current_token.type == EQUALS
            eat!(parser, EQUALS)
            param_value = parse_expression(parser)
            push!(parameters, param_name => param_value)
        else
            push!(parameters, param_name => nothing)
        end

        while parser.current_token.type == COMMA
            eat!(parser, COMMA)
            param_name = parser.current_token.value
            eat!(parser, IDENTIFIER)

            if parser.current_token.type == EQUALS
                eat!(parser, EQUALS)
                param_value = parse_expression(parser)
                push!(parameters, param_name => param_value)
            else
                push!(parameters, param_name => nothing)
            end
        end
    end

    eat!(parser, RPAREN)

    children = ASTNode[]
    eat!(parser, LBRACE)

    while parser.current_token.type != RBRACE
        push!(children, parse_statement(parser))
    end
    
    eat!(parser, RBRACE)
    
    return ModuleNode(name, parameters, children)
end

function parse_call(parser::Parser)
    name = parser.current_token.value
    eat!(parser, IDENTIFIER)

    arguments = Pair{String, Any}[]
    eat!(parser, LPAREN)

    if parser.current_token.type != RPAREN
        if parser.current_token.type == IDENTIFIER && peek_next(parser.lexer) == '='
            arg_name = parser.current_token.value
            eat!(parser, IDENTIFIER)
            eat!(parser, EQUALS)
            arg_value = parse_expression(parser)
            push!(arguments, arg_name => arg_value)
        else
            arg_value = parse_expression(parser)
            push!(arguments, "" => arg_value)
        end
        
        while parser.current_token.type == COMMA
            eat!(parser, COMMA)
            if parser.current_token.type == IDENTIFIER && peek_next(parser.lexer) == '='
                arg_name = parser.current_token.value
                eat!(parser, IDENTIFIER)
                eat!(parser, EQUALS)
                arg_value = parse_expression(parser)
                push!(arguments, arg_name => arg_value)
            else
                arg_value = parse_expression(parser)
                push!(arguments, "" => arg_value)
            end
        end
    end
    
    eat!(parser, RPAREN)
    
    children = ASTNode[]
    if parser.current_token.type == LBRACE
        eat!(parser, LBRACE)
        while parser.current_token.type != RBRACE
            push!(children, parse_statement(parser))
        end
        eat!(parser, RBRACE)
    else
        eat!(parser, SEMICOLON)
    end
    
    return CallNode(name, arguments, children)
end

function parse_statement(parser::Parser)
    if parser.current_token.type == KEYWORD && parser.current_token.value == "module"
        return parse_module(parser)
    elseif parser.current_token.type == IDENTIFIER
        return parse_call(parser)
    else
        error("Unexpected token in statement: $(parser.current_token.type)\n$(parser.current_token)")
    end
end

function parse_scad(input::String)
    parser = Parser(input)
    statements = ASTNode[]

    while parser.current_token.type != EOF
        push!(statements, parse_statement(parser))
    end

    return statements
end

# Helper function to print the AST
function print_ast(node::ASTNode, indent::Int = 0)
    indent_str = " " ^ indent

    if node isa ModuleNode
        println("$(indent_str)Module($(node.name))")
        println("$(indent_str)Parameters:")
        for param in node.parameters
            println("$(indent_str)  $(param.first) = $(param.second === nothing ? "undef" : param.second)")
        end
        println("$(indent_str)Children:")
        for child in node.children
            print_ast(child, indent + 2)
        end
    elseif node isa CallNode
        println("$(indent_str)Call($(node.name))")
        println("$(indent_str)Arguments:")
        for arg in node.arguments
            println("$(indent_str)  $(arg.first == "" ? "" : "$(arg.first) = ")$(arg.second)")
        end
        if !isempty(node.children)
            println("$(indent_str)Children:")
            for child in node.children
                print_ast(child, indent + 2)
            end
        end
    elseif node isa NumberNode
        println("$(indent_str)Number($(node.value))")
    elseif node isa VectorNode
        println("$(indent_str)Vector([")
        for value in node.values
            print_ast(value, indent + 2)
        end
        println("$(indent_str)])")
    elseif node isa VariableNode
        println("$(indent_str)Variable($(node.name))")
    else
        error("Unknown node type: $node")
    end
end

# Add these functions just before the final "end # module":

function write_scad(nodes::Vector{ASTNode}, indent::Int = 0)
    join(map(node -> write_scad(node, indent), nodes), "\n\n")
end

function write_scad(node::ModuleNode, indent::Int = 0)
    indent_str = " " ^ indent
    
    # Format parameters
    params = String[]
    for param in node.parameters
        if param.second === nothing
            push!(params, param.first)
        else
            push!(params, "$(param.first) = $(write_scad(param.second, 0))")
        end
    end
    params_str = join(params, ", ")

    # Format module header
    result = "$(indent_str)module $(node.name)($(params_str)) {\n"

    # Format children with increased indent
    if !isempty(node.children)
        children_str = join(map(child -> write_scad(child, indent + 4), node.children), "\n")
        result *= children_str * "\n"
    end

    result *= "$(indent_str)}"
    return result
end

function write_scad(node::CallNode, indent::Int = 0)
    indent_str = " " ^ indent
    
    # Format arguments
    args = String[]
    for arg in node.arguments
        if arg.first == ""
            push!(args, write_scad(arg.second, 0))
        else
            push!(args, "$(arg.first) = $(write_scad(arg.second, 0))")
        end
    end
    args_str = join(args, ", ")

    result = "$(indent_str)$(node.name)($(args_str))"

    # Handle children if present
    if !isempty(node.children)
        result *= " {\n"
        children_str = join(map(child -> write_scad(child, indent + 4), node.children), "\n")
        result *= children_str * "\n$(indent_str)}"
    else
        result *= ";"
    end

    return result
end

function write_scad(node::NumberNode, indent::Int = 0)
    indent_str = " " ^ indent
    "$(indent_str)$(node.value)"
end

function write_scad(node::VectorNode, indent::Int = 0)
    indent_str = " " ^ indent
    values = map(value -> write_scad(value, 0), node.values)
    "$(indent_str)[$(join(values, ", "))]"
end

function write_scad(node::VariableNode, indent::Int = 0)
    indent_str = " " ^ indent
    "$(indent_str)$(node.name)"
end
