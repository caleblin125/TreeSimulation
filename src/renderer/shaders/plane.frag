#version 330 core
out vec4 FragColor;
in vec4 gl_FragCoord;

void main()
{
    float r = sin(0.55f + gl_FragCoord.x*0.1f)*0.5f + 0.5f;
    float g = sin(0.55f + gl_FragCoord.y*0.1f)*0.5f + 0.5f;
    FragColor = vec4(r, g, 1.0f, 1.0f);
}
