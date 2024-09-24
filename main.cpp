#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <regex>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

using namespace std;

const unsigned int WIDTH = 800;
const unsigned int HEIGHT = 600;

float mouse[2] = {0, 0};
bool mousePressed = false;

// Function to read a shader file
string readShaderFile(const string& filePath) {
    ifstream shaderFile(filePath);
    if (!shaderFile.is_open()) {
		throw runtime_error("Error: Could not open shader file ");
    }
    stringstream buffer;
    buffer << shaderFile.rdbuf();
    return buffer.str();
}

// Function to compile a shader
GLuint compileShader(GLenum type, const string& source) {
    GLuint shader = glCreateShader(type);
    const char* sourceCStr = source.c_str();
    glShaderSource(shader, 1, &sourceCStr, nullptr);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        cerr << "Error: Shader Compilation Failed\n" << infoLog << endl;
        glDeleteShader(shader);
        return 0;
    }
    return shader;
}

// Function to create a shader program
GLuint createProgram(const string& vertexShaderSource, const string& fragmentShaderSource) {
    GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
    GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);

    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return program;
}

// Mouse position callback
void mouseMoveCallback(GLFWwindow* window, double xpos, double ypos) {
    if (mousePressed) {
        mouse[0] = xpos;
        mouse[1] = HEIGHT - ypos; // Invert Y-axis
    }
}

// Mouse button callback
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        mousePressed = (action == GLFW_PRESS);
    }
}


string processIncludes(const string& shaderSource) {
    string processedSource = shaderSource;
    std::regex includeRegex(R"(#include\s*["<]([^">]+)[">])");
    std::smatch match;

    // Loop until no more includes are found
    while (std::regex_search(processedSource, match, includeRegex)) {
        string includePath = "shaders/" + match[1].str();
        string includedSource = readShaderFile(includePath);

        // Replace the include directive with the included source
        processedSource.replace(match.position(0), match.length(0), includedSource);
    }

    return processedSource;
}


// Main function
int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        cerr << "Error: GLFW Initialization Failed" << endl;
        return -1;
    }

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "OpenGL Shader", nullptr, nullptr);
    if (!window) {
        cerr << "Error: Window Creation Failed" << endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        cerr << "Error: GLEW Initialization Failed" << endl;
        return -1;
    }

    // Resize the viewport
    glViewport(0, 0, WIDTH, HEIGHT);

    // Create a fullscreen quad
    float positions[] = {
        -1.0f, -1.0f,
        1.0f, -1.0f,
        -1.0f,  1.0f,
        -1.0f,  1.0f,
        1.0f, -1.0f,
        1.0f,  1.0f
    };

    GLuint positionBuffer;
    glGenBuffers(1, &positionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(positions), positions, GL_STATIC_DRAW);

    // Vertex Shader source
    const string vertexShaderSource = R"(
		attribute vec2 a_position;  // Declare a_position as an attribute
		varying vec2 v_texCoord;     // Declare v_texCoord as a varying

		void main() {
			gl_Position = vec4(a_position, 0.0, 1.0);
			v_texCoord = a_position * 0.5 + 0.5; // Map to 0.0 to 1.0
		}
    )";

    // Fragment Shader source (from file or inline)
    string fragmentShaderSource =  R"(
		#define NONC

		uniform float iTime;
		uniform vec2 iResolution;
		uniform vec4 iMouse;
	 )"
	+ processIncludes(readShaderFile("shaders/main_scene.glsl"))
	+  R"(

		void main() {
			vec4 fragColor;
			mainImage(fragColor, gl_FragCoord.xy); // Call mainImage with output and coordinates
			gl_FragColor = fragColor;
		}
	 )";

    // Create shader program
    GLuint program = createProgram(vertexShaderSource, fragmentShaderSource);
    glUseProgram(program);

    // Set up the attributes
    GLuint positionAttributeLocation = glGetAttribLocation(program, "a_position");
    glEnableVertexAttribArray(positionAttributeLocation);
    glBindBuffer(GL_ARRAY_BUFFER, positionBuffer);
    glVertexAttribPointer(positionAttributeLocation, 2, GL_FLOAT, GL_FALSE, 0, nullptr);

    // Set mouse callbacks
    glfwSetCursorPosCallback(window, mouseMoveCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);

    // Main render loop
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        // Set uniform variables
        glUniform2f(glGetUniformLocation(program, "iResolution"), WIDTH, HEIGHT);
        glUniform1f(glGetUniformLocation(program, "iTime"), glfwGetTime());
        glUniform4f(glGetUniformLocation(program, "iMouse"), mouse[0], mouse[1], 0.0f, 0.0f);

        glDrawArrays(GL_TRIANGLES, 0, 6);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Cleanup
    glDeleteBuffers(1, &positionBuffer);
    glDeleteProgram(program);
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
