#ifndef NONC
precision mediump float;
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
# define iTime u_time
# define iResolution u_resolution
# define iMouse u_mouse
void mainImage(out vec4 color,in vec2 pxy);
void main() {
	mainImage(gl_FragColor, gl_FragCoord.xy);
}
#endif
