//#version 120
#extension GL_ARB_gpu_shader5 : enable

out INFO {
    
    vec2 texCoords;
    float k1_red;
    float k1_green;
    float k1_blue;
    vec2 fullscr_center;
    vec2 left_center;
    vec2 right_center;

} vs_out;

uniform float k1_red;
uniform float k1_green;
uniform float k1_blue;
uniform vec2 fullscr_center;
uniform vec2 left_center;
uniform vec2 right_center;

void main() {
    
    const vec2 texCoords[4] = vec2[4](  vec2(0.0, 1.0),
                                        vec2(0.0, 0.0),
                                        vec2(1.0, 0.0),
                                        vec2(1.0, 1.0));

    gl_Position = ftransform();
    gl_TexCoord[0] = gl_MultiTexCoord0;
    vs_out.texCoords = texCoords[gl_VertexID];
    
    vs_out.k1_red = k1_red;
    vs_out.k1_green = k1_green;
    vs_out.k1_blue = k1_blue;
    vs_out.fullscr_center = fullscr_center;
    vs_out.left_center = left_center;
    vs_out.right_center = right_center;
}
