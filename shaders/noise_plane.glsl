// Modeling - 2024.09.15
// Eric Galin
#include "_common.glsl"
#line 5

struct Ray{
    vec3 o;// Origin
    vec3 d;// Direction
};

struct Val {
  float v; // Signed distance
  int c; // Cost
};

// Compute point on ray
// ray : The ray
//   t : Distance
vec3 Point(Ray ray,float t)
{
    return ray.o+t*ray.d;
}

// Random direction in a hemisphere
// seed : Integer seed, from 0 to N
//    n : Direction of the hemisphere
vec3 Hemisphere(int seed,vec3 n)
{
    float a=fract(sin(176.19*float(seed)));// Uniform randoms
    float b=fract(sin(164.19*float(seed)));

    float u=2.*3.1415*a;// Random angle
    float v=acos(2.*b-1.);// Arccosine distribution to compensate at poles

    vec3 d=vec3(cos(u)*cos(v),sin(u)*cos(v),sin(v));// Direction
    if(dot(d,n)<0.){d=-d;}// Hemisphere

    return d;
}

// Camera -------------------------------------------------------------------------------

// Rotation matrix around z axis
// a : Angle
mat3 Rz(float a)
{
  float sa=sin(a);float ca=cos(a);
  return mat3(ca,sa,0.,-sa,ca,0.,0.,0.,1.);
}

// Compute the ray
//      m : Mouse position
//      p : Pixel
Ray CreateRay(vec2 m,vec2 p)
{
  float a=3.*3.14*m.x;
  float le=3.5;

  // Origin
  vec3 ro=vec3(37.,0.,15.);
  ro*=Rz(a);

  // Target point
  vec3 ta=vec3(0.,0.,1.);

  // Orthonormal frame
  vec3 w=normalize(ta-ro);
  vec3 u=normalize(cross(w,vec3(0.,0.,1.)));
  vec3 v=normalize(cross(u,w));
  vec3 rd=normalize(p.x*u+p.y*v+le*w);
  return Ray(ro,rd);
}

// Texture --------------------------------------------------------------------------------------

// Hashing function
// Returns a random number in [-1,1]
// p : Vector in space
float Hash(in vec3 p)
{
    p  = fract( p*0.3199+0.152 );
	p *= 17.0;
    return fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
}

// Procedural value noise with cubic interpolation
// x : Point
float Noise(in vec3 p)
{
    vec3 i = floor(p);
    vec3 f = fract(p);

    f = f*f*(3.0-2.0*f);
    // Could use quintic interpolation instead of cubic
    // f = f*f*f*(f*(f*6.0-15.0)+10.0);

    return mix(mix(mix( Hash(i+vec3(0,0,0)),
                        Hash(i+vec3(1,0,0)),f.x),
                   mix( Hash(i+vec3(0,1,0)),
                        Hash(i+vec3(1,1,0)),f.x),f.y),
               mix(mix( Hash(i+vec3(0,0,1)),
                        Hash(i+vec3(1,0,1)),f.x),
                   mix( Hash(i+vec3(0,1,1)),
                        Hash(i+vec3(1,1,1)),f.x),f.y),f.z);
}

// turbulence (somme de fonctions de bruit) pour perturber la surface
float Turbulence(vec3 p)
{
    float deformation = 0.0;

    deformation += 1.0 * Noise(p / 2.0);
    deformation += 0.5 * Noise(p / 1.0);
    deformation += 0.25 * Noise(p / 0.5);

    deformation /= 1.6;  // pour lisser un peu

    return deformation;
}


// Primitives -------------------------------------------------------------------------------

// Sphere
// p : point
// c : center of skeleton
// r : radius
Val Sphere(vec3 p,vec3 c,float r)
{
  return Val(length(p-c)-r,1);
}

// Plane
// p : point
// n : Normal of plane
// o : Point on plane
Val Plane(vec3 p, vec3 n, vec3 o)
{
    return Val(dot((p-o),n),1);
}

// Operators

// Union
// a,b : field function of left and right sub-trees
Val Union(Val a,Val b)
{
  return Val(min(a.v,b.v),a.c+b.c+1);
}

// Potential field of the object
// p : point
Val object(vec3 p)
{
    Val montagne = Plane(p,vec3(0.,0.,1.),vec3(0.0,0.0,4.0));

  montagne.v += 5.0 * Turbulence(p);

  return montagne;
}


// Analysis of the scalar field -----------------------------------------------------------------

const int Steps=200;// Number of steps
const float Epsilon=.01;// Marching epsilon

// Object normal
// p : point
vec3 ObjectNormal(vec3 p)
{
  const float eps=.001;
  vec3 n;
  Val val=object(p);
  float v=val.v;
  n.x=object(vec3(p.x+eps,p.y,p.z)).v-v;
  n.y=object(vec3(p.x,p.y+eps,p.z)).v-v;
  n.z=object(vec3(p.x,p.y,p.z+eps)).v-v;
  return normalize(n);
}

// Trace ray using ray marching
// ray : The ray
//   e : Maximum distance
//   h : hit
//   s : Number of steps
//   c : cost
bool SphereTrace(Ray ray,float e,out float t,out int s,out int c)
{
  bool h=false;

  // Start at the origin
  t=0.0;
  c=0;

  for(int i=0;i<Steps;i++)
  {
    s=i;
    vec3 p=Point(ray,t);
    Val val=object(p);
    float v=val.v;
    c+=val.c;
    // Hit object
    if(v<0.)
    {
      h=true;
      break;
    }
    // Move along ray
    t+=max(Epsilon,v);
    // Escape marched too far away
    if(t>e)
    {
      break;
    }
  }
  return h;
}

// Lighting -------------------------------------------------------------------------------

// Background color
// ray : Ray
vec3 background(Ray ray)
{
  return mix(vec3(0.4, 0.6, 1.0), vec3(0.6, 0.8, 1.0), ray.d.z * 0.5 + 0.5);
}

// Shadowing
// p : Point
// n : Normal
// l : Light direction
float Shadow(vec3 p,vec3 n,vec3 l)
{
  float t;
  int s;
  int c;
  bool hit=SphereTrace(Ray(p+Epsilon*n,l),100.,t,s,c);
  if(!hit)
  {
    return 1.;
  }
  return 0.;
}

// Shading and lighting
//   p : Point
//   n : Normal at point
// eye : Eye direction
vec3 Shade(vec3 p, vec3 n, Ray eye)
{
  // Point light
  const vec3 lp = vec3(5., 10., 25.);

  // Light direction to point light
  vec3 l = normalize(lp - p);

  // Ambient color
  vec3 ambient = .25 + .25 * background(Ray(p, n));

  // Shadow computation
  float shadow = Shadow(p, n, l);

  // Phong diffuse =>
  vec3 mountainColor = vec3(0.6, 2., 0.1);
  vec3 diffuse = .35 * clamp(dot(n, l), 0., 1.) * mountainColor;


  // Specular
  vec3 r = reflect(eye.d, n);
  vec3 specular = .15 * pow(clamp(dot(r, l), 0., 1.), 35.) * vec3(1., 1., 1.);

  vec3 c = ambient + shadow * (diffuse + specular);

  return c;
}

// Shading according to the number of steps in sphere tracing
// n : Number of steps
vec3 ShadeSteps(int n,int m)
{
  float t=float(n)/(float(m));
  return.5+mix(vec3(.05,.05,.5),vec3(.65,.39,.65),t);
}

// Image
void mainImage(out vec4 color,in vec2 pxy)
{
  // Convert pixel coordinates
    vec2 pixel=(-iResolution.xy+2.*pxy)/iResolution.y;

  // Mouse
  vec2 m=iMouse.xy/iResolution.xy;

  // Camera
  Ray ray=CreateRay(m,pixel);

  // Trace ray

  // Hit and number of steps
  float t=0.0;
  int s=0;
  int c;
  bool hit=SphereTrace(ray,100.,t,s,c);

  // Shade background
  vec3 rgb=background(ray);

  if(hit)
  {
    // Position
    vec3 p=Point(ray,t);

    // Compute normal
    vec3 n=ObjectNormal(p);

    // Shade object with light
    rgb=Shade(p,n,ray);
  }

  color=vec4(rgb,1.);
}

