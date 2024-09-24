GLEW_DIR	:= /usr/local/Cellar/glew/2.2.0_1

GLFW_DIR	:= /usr/local/Cellar/glfw/3.4

CXXFLAGS	:= 	-I$(GLEW_DIR)/include -L$(GLEW_DIR)/lib \
				-I$(GLFW_DIR)/include -L$(GLFW_DIR)/lib	\
				-lglfw -lglew -framework OpenGL

CXX			:= clang++

all: main

web:
	php -S 127.0.0.1:8080 &\
	open http://127.0.0.1:8080

main: main.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(main)

re:	clean main
