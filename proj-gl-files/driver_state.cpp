#include "driver_state.h"
#include <cstring>
#include <algorithm>

using namespace std;

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;
    state.image_color = 0;
    state.image_color = new pixel[width * height];
    state.image_depth = 0;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    for (int i = 0; i < (width * height); i++) {
      state.image_color[i] = make_pixel(0, 0, 0);
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //std::cout<<"TODO: implement rendering."<<std::endl;
    //basically had to change this whole thing -___-
    data_geometry *tri_array = new data_geometry[3];
    data_vertex ver; //{};
    float *p = state.vertex_data;

    switch (type){
        case render_type::triangle:
            for (int i = 0; i < state.num_vertices/3; i++) {
                for (int j = 0; j < 3; j++) {
                    tri_array[j].data = p;
                    p = p + state.floats_per_vertex;
                }

                for (int k = 0; k < 3; k++) {
                    ver.data = tri_array[k].data;
                    state.vertex_shader(ver, tri_array[k], state.uniform_data);
                }

                rasterize_triangle(state, (const data_geometry**) &tri_array);
            }
        break;

        case render_type::indexed:
        break;

        case render_type::fan:
        break;

        case render_type::strip:
        break;

        default:
        break;
    }
    //delete [] tri_array;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    //vec4 left, right, top, bottom, far, near;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    // data_geometry* out = new data_geometry[3];
    // int ax, ay, bx, by, cx, cy;
    // int px, py;
  
    // data_vertex ver;
    //int i, j;
    //int x[3], y[3];
    unsigned image_index;
    int ax, ay, bx, by, cx, cy;
    float area_abc, area_pbc, area_apc, area_abp;
    float alpha, beta, gamma;
    int picCoords[3][2];
    
    for (int k = 0; k < 3; k++) {
        //ver.data = in[k]->data;
        //state.vertext_shader(ver, out[k], state.uniform_data);
        //state.vertex_shader(ver, out[k], state.uniform_data);

        picCoords[k][0] = (state.image_width/2) * (in[k]->gl_Position[0]/in[k]->gl_Position[3]) + (state.image_width/2) - 0.5;
        picCoords[k][1] = (state.image_height/2) * (in[k]->gl_Position[1]/in[k]->gl_Position[3]) + (state.image_height/2) - 0.5;

        // x[k] = i;
        // y[k] = j;
        
        image_index = picCoords[k][0] + picCoords[k][1] * state.image_width;
        //state.image_color[image_index] = make_pixel(255, 255, 255);
        
    }
 
    // area_abc = 0.5 * (((bx * cy) - (cx * by))-((ax * cy) - (cx * ay)) - ((ax * by)-(bx * ay)));
    ax = picCoords[0][0];
    ay = picCoords[0][1];

    bx = picCoords[1][0];
    by = picCoords[1][1];

    cx = picCoords[2][0];
    cy = picCoords[2][1];
    
    area_abc = 0.5 * ((bx * cy - cx * by) - (ax * cy - cx * ay) + (ax * by - bx * ay));
    
    for (int j = 0; j < state.image_height; j++) {
        for (int i = 0; i < state.image_width; i++) {
            area_pbc = 0.5 * ((bx * cy - cx * by) + (by - cy) * i + (cx - bx) * j);
            area_apc = 0.5 * ((cx * ay - ax * cy) + (cy - ay) * i + (ax - cx) * j);
            area_abp = 0.5 * ((ax * by - bx * ay) + (ay - by) * i + (bx - ax) * j);

            alpha = area_pbc/area_abc;
            beta = area_apc/area_abc;
            gamma = area_abp/area_abc;
        
            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                int image_index1 = i + j * state.image_width;
                //state.image_color[image_index] = make_pixel(255, 255, 255);
                auto *data = new float[MAX_FLOATS_PER_VERTEX];
                data_fragment frag{data};
                data_output output;

                for (int k = 0; k < state.floats_per_vertex; k++) {
                    float k_gour;
                    switch(state.interp_rules[k]) {
                        case interp_type::flat:
                            frag.data[k] = in[0]->data[k];
                        break;

                        case interp_type::smooth:
                            k_gour = (alpha/in[0]->gl_Position[3] + beta/in[1]->gl_Position[3] + gamma/in[2]->gl_Position[3]);

                            alpha /= (k_gour * (in[0]->gl_Position[3]));
                            beta /= (k_gour * (in[1]->gl_Position[3]));
                            gamma /= (k_gour * (in[2]->gl_Position[3]));
                        break;

                        case interp_type::noperspective:
                            frag.data[k] = alpha * in[0]->data[i] + beta * in[1]->data[i] + gamma * in[2]->data[i];
                        break;

                        default:
                        break;
                    }
                }
                state.fragment_shader((const data_fragment)frag, output, state.uniform_data);
                output.output_color = output.output_color * 255;
                state.image_color[image_index1] = make_pixel(output.output_color[0], output.output_color[1], output.output_color[2]);
            }
        }
    }
}

