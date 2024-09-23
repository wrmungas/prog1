// Combines all the separate files into one big file 
// I'm trying this to avoid tons of problems with importing since I'm VERY new to JS
// (and the original file did this anyway)
// 
// This way I can look develop individual classes/components separately and then integrate them here
//
// This file will include main() as well as all classes and functions I wrote for the program
//
/* FROM 3D MATH */

class Vector {
    x;
    y;
    z;

    constructor() {
        this.x = 0;
        this.y = 0;
        this.z = 0;
    }

    set(x, y, z) {
        if(typeof(x) != "number" || typeof(y) != "number" || typeof(z) != "number") {
            throw "vector component not a number";
        }
        else {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    }

    copy() {
        var v = new Vector();
        v.set(this.x, this.y, this.z);
        return v;
    }
}

// creates a vector with the given x, y, and z components; shortening of syntax
function create(x, y, z) {
    if(typeof(x) != "number" || typeof(y) != "number" || typeof(z) != "number") {
        throw "vector component not a number";
    }
    var result = new Vector();
    result.set(x, y, z);
    return result;
}

// gives the magnitude of the vector
function mag(v) {
    if(!(v instanceof Vector)) {
        throw "invalid argument: norm() requires a Vector";
    }
    return Math.sqrt(
        v.x * v.x + 
        v.y * v.y + 
        v.z * v.z );
}

// normalizes the vector
function norm(v) {
    if(!(v instanceof Vector)) {
        throw "invalid argument: norm() requires a Vector";
    }

    var result = v.copy();
    var m = mag(result);
    result.x /= m;
    result.y /= m;
    result.z /= m;

    return result;
}

// computes vec1 * vec2 and returns the result
function dot(vec1, vec2) {
    if( !(vec1 instanceof Vector)) {
        throw "invalid argument 1: dot() requires a Vector"
    } 
    if( !(vec2 instanceof Vector)) {
        throw "invalid argument 2: dot() requires a Vector";
    }
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}

// computes vec1 x vec2 and returns the result
function cross(vec1, vec2) {
    if( !(vec1 instanceof Vector && vec2 instanceof Vector)) {
        throw "invalid arguments: cross() is for two Vectors";
    }
    var i = vec1.y * vec2.z - vec1.z * vec2.y;
    var j = vec1.z * vec2.x - vec1.x * vec2.z;
    var k = vec1.x * vec2.y - vec1.y * vec2.x;
    
    var result = new Vector();
    result.set(i, j, k);
    return result;
}

// computes vec1 + vec2 and returns result
function add(vec1, vec2) {
    if( !(vec1 instanceof Vector && vec2 instanceof Vector)) {
        throw "invalid arguments: add() is for two Vectors";
    }
    var result = new Vector();
    result.set(
        vec1.x + vec2.x,
        vec1.y + vec2.y,
        vec1.z + vec2.z );
    
    return result;
}


// computes vec1 - vec2 and returns result: (vec1 - vec2), or a vector from vec2 to vec1
function sub(vec1, vec2) {
    if( !(vec1 instanceof Vector && vec2 instanceof Vector)) {
        throw "invalid arguments: sub() is for two Vectors";
    }
    var result = new Vector();
    result.set(
        vec1.x - vec2.x, 
        vec1.y - vec2.y, 
        vec1.z - vec2.z);
    
    return result;
}

// multiplies the components of vec1 by constant a and stores in vecresult
function mul(vec1, a) {
    if( !(vec1 instanceof Vector) ) {
        throw "invalid arg1: mul is for a Vector";
    }
    if(typeof(a) !== "number") {
        throw "invalid arg2: mul requires a number";
    }
    var result = new Vector();
    result.set(
        vec1.x * a, 
        vec1.y * a, 
        vec1.z * a);

    return result;
}


// END 3d Math
// BEGIN Color


// Color constructor
class Color {
    constructor(r,g,b,a) {
        try {
            if ((typeof(r) !== "number") || (typeof(g) !== "number") || (typeof(b) !== "number") || (typeof(a) !== "number"))
                throw "color component not a number";
            else if ((r<0) || (g<0) || (b<0) || (a<0)) 
                throw "color component less than 0";
            else if ((r>255) || (g>255) || (b>255) || (a>255)) 
                throw "color component bigger than 255";
            else {
                this.r = r; this.g = g; this.b = b; this.a = a; 
            }
        } // end try
        
        catch (e) {
            console.log(e);
        }
    } // end Color constructor

        // Color change method
    change(r,g,b,a) {
        try {
            if ((typeof(r) !== "number") || (typeof(g) !== "number") || (typeof(b) !== "number") || (typeof(a) !== "number"))
                throw "color component not a number";
            else if ((r<0) || (g<0) || (b<0) || (a<0)) 
                throw "color component less than 0";
            else if ((r>255) || (g>255) || (b>255) || (a>255)) 
                throw "color component bigger than 255";
            else {
                this.r = r; this.g = g; this.b = b; this.a = a; 
            }
        } // end throw
        
        catch (e) {
            console.log(e);
        }
    } // end Color change method
} // end color class


// END Color
// BEGIN Material


// defines a class for material properties of objects in the scene
class Material {
    ambient;
    diffuse;
    specular;
    n;
    constructor(ambient, diffuse, specular, n) {
        if(!(ambient instanceof Color) || !(diffuse instanceof Color) || !(specular instanceof Color)) {
            throw "material is defined by ambient, diffuse, and specular colors";
        }
        if(typeof(n) !== "number") {
            throw "material needs numerical n";
        }

        this.ambient = ambient;
        this.diffuse = diffuse;
        this.specular = specular;
        this.n = n;
    }
}


// END Material
// BEGIN Light


// define a class for lights

class Light {
    constructor(position, ambient, diffuse, specular) {
        if(!position instanceof Vector) {
            throw "Light must have a 3d vector position";
        }
        if(!ambient instanceof Color || !diffuse instanceof Color || !specular instanceof Color) {
            throw "Light must be composed of ambient, diffuse, and specular color components";
        }
        this.position = position;
        this.ambient = ambient;
        this.diffuse = diffuse;
        this.specular = specular;
    }
}


// END Light
// BEGIN Triangle


// Define a triangle class and functions to read a list of triangles

class Triangle {
    constructor(vertex1, vertex2, vertex3, material) {
        this.vertex1 = vertex1;
        this.vertex2 = vertex2;
        this.vertex3 = vertex3;
        this.material = material;
        this.normal = cross( sub(vertex1, vertex2), sub(vertex3, vertex2) );
    }
}

function getInputTriangles(url) {
    // load the triangles file
    var httpReq = new XMLHttpRequest(); // a new http request
    httpReq.open("GET",url,false); // init the request
    httpReq.send(null); // send the request
    var startTime = Date.now();
    while ((httpReq.status !== 200) && (httpReq.readyState !== XMLHttpRequest.DONE)) {
        if ((Date.now()-startTime) > 3000)
            break;
    } // until its loaded or we time out after three seconds
    if ((httpReq.status !== 200) || (httpReq.readyState !== XMLHttpRequest.DONE)) {
        console.log*("Unable to open input triangles file!");
        return String.null;
    } else
        return JSON.parse(httpReq.response); 
} // end get input triangles

function read_tris(url) {
    var data = getInputTriangles(url); 

    if(data == String.null) {
        return [];
    }
    
    var triangles_list = [];

    // loop over all data sections
    for(var i = 0; i < data.length; i++) {

        //get the colors and n for the current object's material
        var diff = new Color(data[i].material.diffuse[0], data[i].material.diffuse[1], data[i].material.diffuse[2], 255);
        var spec = new Color(data[i].material.specular[0], data[i].material.specular[1], data[i].material.specular[2], 255);
        var amb = new Color(data[i].material.ambient[0], data[i].material.ambient[1], data[i].material.ambient[2], 255);
        var n = data[i].material.n;

        // combine into a single material
        var current_material = new Material(amb, diff, spec, n);

        // loop over the all triangles
        var num_tris = data[i].triangles.length;
        for(var j = 0; j < num_tris; j++) {

            // prep the vectors for each vertex

            var vertices = [];

            var indices = [data[i].triangles[j][0], data[i].triangles[j][1], data[i].triangles[j][2]];
           
            for(var k = 0; k < 3; k++) {
                vertices.push(create(data[i].vertices[indices[k]][0], data[i].vertices[indices[k]][1], data[i].vertices[indices[k]][2]));
            }

            // create the triangle from the three vertices
            var triangle = new Triangle(vertices[0], vertices[1], vertices[2], current_material);

            // add it to the list of triangles
            triangles_list.push(triangle);
        }
    }

    // return the list of triangles
    return triangles_list;
}

/* END TRIANGLE */
//
//
//
/* FROM RAY */


// class defining a ray
class Ray {
    dir;
    origin;

    constructor(dir, origin) {
        if( !( dir instanceof Vector && origin instanceof Vector)) {
            throw "error: ray must be composed of vectors";
        }
        this.dir = norm(dir); // useful for intersection distances
        this.origin = origin;
    }
}


/* END RAY */
//
//
//
/* HELPER FUNCTIONS */ 

// I found pairing rays with corresponding 2d pixel data useful
class Pixel_Ray {
    ray;
    x;
    y;

    constructor(ray, x, y) {
        if( !(ray instanceof Ray) || typeof(x) !== "number" || typeof(y) !== "number") {
            throw "pixel needs a ray and corresponding x/y screen coordinates";
        }
        this.ray = ray;
        this.x = x;
        this.y = y;
    }
}

function side(n, i, v1, v2) {
    var result = dot(n, cross(sub(i, v1), sub(v2, v1)));
    if(result < 0) {
        return '-';
    }
    else {
        return '+';
    }
}

// checks whether a ray and a triangle intersect (returns -1 if not)
function check_intersection(ray, tri) {
    //validation
    if( !(ray instanceof Ray) ) {
        throw "intersection test requires a Ray";
    }
    if( !(tri instanceof Triangle)) {
        throw "intersection test requires a Triangle"
    }

    //variables setup
    var epsilon = Number.epsilon; // for float errors
    
    var n = tri.normal;
    var point = tri.vertex1;

    var a = tri.vertex1;
    var b = tri.vertex2;
    var c = tri.vertex3;

    var origin = ray.origin;
    var dir = ray.dir; 

    // first check that the line intersects the plane of the triangle
    if(Math.abs( dot(dir, n) ) < epsilon ) {
        return -1;
    }

    // then find the point of ray-plane intersection
    var d = dot(n, point);
    var t = (d - dot(n, origin)) / (dot(n, dir));
    
    // intersection point
    var i = add(origin, mul(dir, t));

    // check that the intersection is within 'world' bounds
    // if( i.x > 1 || i.x < 0 || i.y > 1 || i.y < 0 || i.z > 1 || i.z < 0 ) {
    //     return -1;
    // }

    // check that this point is on the same side of all of the triangle's edges;
    // if it is, it is on the interior of the triangle
    var side_ab = side(n, i, a, b);
    var side_bc = side(n, i, b, c);
    var side_ca = side(n, i, c, a);
    if(side_ab === side_bc && side_bc === side_ca && side_ca === side_ab) {
        return t;
    }
    // otherwise it is outside the triangle
    return -1;
}

// draw a pixel at x,y using color
function drawPixel(imagedata,x,y,color) {
    try {
        if ((typeof(x) !== "number") || (typeof(y) !== "number"))
            throw "drawpixel location not a number";
        else if ((x<0) || (y<0) || (x>=imagedata.width) || (y>=imagedata.height))
            throw "drawpixel location outside of image";
        else if (color instanceof Color) {
            var pixelindex = (y*imagedata.width + x) * 4;
            imagedata.data[pixelindex] = color.r;
            imagedata.data[pixelindex+1] = color.g;
            imagedata.data[pixelindex+2] = color.b;
            imagedata.data[pixelindex+3] = color.a;
        } else 
            throw "drawpixel color is not a Color";
    } // end try
    
    catch(e) {
        console.log(e);
    }
} // end drawPixel


// function to interpret a view and create the pixels/rays for it
function create_pixel_ray_pairs(eye, facing, up, canvas) {
    var view_center = add(eye, mul(norm(facing), 0.5));
    var view_i = norm(cross(facing, up));
    var view_j = norm(cross(view_i, facing));

    // create a ray for every pixel in the display
    var pairs= [];
    for(var i = 0; i < canvas.width; i++) {
        for(var j = 0; j < canvas.height; j++) {
            var origin = eye.copy();

            // point on the plane that this pixel's ray passes through
            var point = add(
                view_center, 
                add(
                    mul(view_i, 0.5 - ((2 * i + 1) / (canvas.width * 2)) ), // 'x' coordinate on the view plane
                    mul(view_j, 0.5 - ((2 * j + 1) / (canvas.height * 2)) ) // 'y' coordinate on the view plane
                ) 
            );
            
            // get the direction that this pixel's ray points (normalize for ease)
            var dir = norm(sub(point, eye));

            var ray = new Ray(dir, origin);

            // get pixels on the screen
            var pair = new Pixel_Ray(ray, i, j);
            
            pairs.push(pair);
        }       
    }

    return pairs;
}

function blinn_phong(ray, tri, dist, light) {
    // ray.dir should already be normalized
    // setup
    var point = add(ray.origin, mul(ray.dir, dist)); // point of intersection

    var v = norm( sub(point, ray.origin) ); // view vector, V
    var l = norm( sub(point, light.position) ); // vector toward light, L
    var n = norm( tri.normal ); // normal of the triangle surface, N
    var h = norm( add(v, l) ); // 'half vector' of V and L (for specular reflections)

    // adjust normal to ensure we have the right one
    if(dot(n, v) < 0) {
        n = mul(n, -1);
    }

    // pull out repeated terms and simplify names
    var N_L = dot(n, l);
    var N_H_n = Math.pow(dot(n, h), tri.material.n);
    var K_a = tri.material.ambient;
    var K_d = tri.material.diffuse;
    var K_s = tri.material.specular;
    var L_a = light.ambient;
    var L_d = light.diffuse;
    var L_s = light.specular;

    // components:[   ambient   ]   [      diffuse      ]   [      specular      ]
    var color_r = (K_a.r * L_a.r) + (K_d.r * L_d.r * N_L) + (K_s.r * L_s.r * N_H_n);
    var color_g = (K_a.g * L_a.g) + (K_d.g * L_d.g * N_L) + (K_s.g * L_s.g * N_H_n);
    var color_b = (K_a.b * L_a.b) + (K_d.b * L_d.b * N_L) + (K_s.b * L_s.b * N_H_n);

    // clamp color multipliers to between 0 and 1
    color_r = Math.min(Math.max(color_r, 0), 1);
    color_g = Math.min(Math.max(color_g, 0), 1);
    color_b = Math.min(Math.max(color_b, 0), 1);

    var result = new Color(
        255 * color_r,
        255 * color_g,
        255 * color_b,
        255 // alpha (unaffected by blinn-phong lighting)
    );

    return result;
}

// function to cast the created rays
function cast_rays_and_draw_pixels(imagedata, pairs, tris, light) {
    // for test
    console.log("Pixel-Ray pairs created: " + pairs.length);
    // bookkeeping setup 
    // record ray-triangle intersection data
    var total_intersections = 0;
    var rays_that_hit = 0;
    var intersections_per_tri = [];
    var hits_per_tri = [];
    for(var i = 0; i < tris.length; i++) {
        intersections_per_tri.push(0);
        hits_per_tri.push(0);
    }

    // loop over all pixel/ray pairs
    for (var i = 0; i < pairs.length; i++) {
        var intersected = false; 
        var closest_idx = -100; // index of the closest hit triangle in the list of triangles
        var closest_dist = -100; // distance along the ray to the closest intersection

        // loop over all triangles to see if the ray intersects any
        for(var j = 0; j < tris.length; j++) {
            var result = check_intersection(pairs[i].ray, tris[j]);
           
            if(result > 0) {
                total_intersections++;
                intersections_per_tri[j]++;

                if(!intersected) {
                    rays_that_hit++;
                    intersected = true;
                    closest_dist = result;
                    closest_idx = j;
                }
                else if(result < closest_dist) {
                    closest_dist = result;
                    closest_idx = j;
                }
            }
        }

        // find the color to draw at this pixel (default black)
        var c = new Color(0, 0, 0, 255);
        // if there was an intersection calculate Blinn-Phong lighting for the hit point
        if(intersected) {
           c = blinn_phong(pairs[i].ray, tris[closest_idx], closest_dist, light);
           hits_per_tri[closest_idx]++;
        }    

        // draw the color
        var pixel_x = pairs[i].x;
        var pixel_y = pairs[i].y;
        drawPixel(imagedata, pixel_x, pixel_y, c);
    }

    // end function with some stats
    console.log("Total line intersections " + total_intersections);
    console.log("Total rays that hit: " + rays_that_hit);
    console.log("Total times each triangle was intersected by a line: \n" + intersections_per_tri);
    console.log("Total times each triangle was hit by a ray: \n" + hits_per_tri);
}

/* END HELPERS */
//
//
//
/* FINALLY: MAIN FUNCTION*/

function main() {
    // Get the canvas and context
    var canvas = document.getElementById("viewport"); 
    var context = canvas.getContext("2d");
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);

    // read the triangles
    const filename = "https://ncsucgclass.github.io/prog1/triangles2.json";
    var tris = read_tris(filename);

    var light = new Light(
        create(-3.0, 1.0, -0.5),
        new Color(1.0, 1.0, 1.0, 255),
        new Color(1.0, 1.0, 1.0, 255),
        new Color(1.0, 1.0, 1.0, 255)
    );

    // for test
    console.log("Triangles read: " + tris.length);
    console.log("Triangle data: \n");
    for(var i = 0; i < tris.length; i++) {
        console.log("Triangles[" + i + "]: " + "{\n" +
            "   vertex1: < " + tris[i].vertex1.x + ", " + tris[i].vertex1.y + ", " + tris[i].vertex1.z + " >\n" +
            "   vertex2: < " + tris[i].vertex2.x + ", " + tris[i].vertex2.y + ", " + tris[i].vertex2.z + " >\n" +
            "   vertex3: < " + tris[i].vertex3.x + ", " + tris[i].vertex3.y + ", " + tris[i].vertex3.z + " >\n" +
            "   material: {\n" +
            "       ambient: [ " + tris[i].material.ambient.r + ", " + tris[i].material.ambient.g + ", " + tris[i].material.ambient.b + " ]\n" +
            "       diffuse: [ " + tris[i].material.diffuse.r + ", " + tris[i].material.diffuse.g + ", " + tris[i].material.diffuse.b + " ]\n" +
            "       specular: [ " + tris[i].material.specular.r + ", " + tris[i].material.specular.g + ", " + tris[i].material.specular.b + " ]\n" +
            "       exponent: " + tris[i].material.n + "\n" +
            "   }\n" +
            "}"
        )
    }

    // setup view
    var eye = create(0.5, 0.5, -0.5); // world coordinates
    var facing = create(0, 0, 1); // direction vector
    var up = create(0, 1, 0); // direction vector

    var pairs = create_pixel_ray_pairs(eye, facing, up, canvas);
    cast_rays_and_draw_pixels(imagedata, pairs, tris, light);

    context.putImageData(imagedata, 0, 0);
}