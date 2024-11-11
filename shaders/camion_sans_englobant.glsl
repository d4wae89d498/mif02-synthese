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

// Object transformations ------------------------------------------------------------------

// Translate a point
// p : the point to translate
// t : translation vector
vec3 Translate(vec3 p, vec3 t) {
    return p - t;
}

// Scale (Homothety)
// p : the point to scale
// s : scaling factor (can be a vector for non-uniform scaling)
vec3 Scale(vec3 p, vec3 s) {
    return p / s;
}

// Rotate around X axis
// p : the point to rotate
// a : angle (in radians) to rotate around the X-axis
vec3 RotateX(vec3 p, float a) {
    float sa = sin(a);
    float ca = cos(a);
    return vec3(p.x, ca * p.y - sa * p.z, sa * p.y + ca * p.z);
}

// Rotate around Y axis
// p : the point to rotate
// a : angle (in radians) to rotate around the Y-axis
vec3 RotateY(vec3 p, float a) {
    float sa = sin(a);
    float ca = cos(a);
    return vec3(ca * p.x + sa * p.z, p.y, -sa * p.x + ca * p.z);
}

// Rotate around Z axis
// p : the point to rotate
// a : angle (in radians) to rotate around the Z-axis
vec3 RotateZ(vec3 p, float a) {
    float sa = sin(a);
    float ca = cos(a);
    return vec3(ca * p.x - sa * p.y, sa * p.x + ca * p.y, p.z);
}

// Transform a point with translation, scaling (homothety), and rotation
// p : the point to transform
// t : translation vector
// s : scaling vector (non-uniform scaling possible)
// r : rotation angles for X, Y, and Z axes (in radians)
vec3 Transform(vec3 p, vec3 t, vec3 s, vec3 r) {
    p = Translate(p, t);
    p = Scale(p, s);
    p = RotateX(p, r.x);
    p = RotateY(p, r.y);
    p = RotateZ(p, r.z);
    return p;
}

// Operators -------------------------------------------------------------------------------

// Union
// a,b : field function of left and right sub-trees
Val Union(Val a,Val b) // cout de la somme des deux objets + 1
{
  return Val(min(a.v,b.v),a.c+b.c+1);
}

// Difference of two signed distance functions
Val Difference(Val a, Val b) { //cout de la somme des deux objets + 1
    return Val(max(a.v, -b.v), a.c + b.c + 1); // Subtract b from a
}

// Intersection of two signed distance functions
Val Intersection(Val a, Val b) { // cout de la somme des deux objets + 1
    return Val(max(a.v, b.v), a.c + b.c + 1); // Keep the larger distance
}


// Primitives -------------------------------------------------------------------------------

// Sphere => renvoie la distance signée d'un point à une sphère
// p : point
// c : center of skeleton
// r : radius
Val Sphere(vec3 p,vec3 c,float r)
{
  return Val(length(p-c)-r,1); // => cout à 1
}

// Plane
// p : point
// n : Normal of plane
// o : Point on plane
Val Plane(vec3 p, vec3 n, vec3 o)
{
    return Val(dot((p-o),n),1); // => cout à 1
}

// Ellipsoid
// p : point to test
// c : center of the ellipsoid
// r : radii along the ellipsoid axes (as a vec3)
Val Ellipsoid(vec3 p, vec3 c, vec3 r) {
    return Val(length((p - c) / r) - 1.0, 5); // => cout à 5
}

// Box
// p : point to test
// c : center of the box
// r : half-extents (size in each direction from the center) of the box
Val Box(vec3 p, vec3 c, vec3 r) {
    vec3 d = abs(p - c) - r;
    return Val(length(max(d, 0.0)) + min(max(d.x, max(d.y, d.z)), 0.0), 3); // => cout à 3
}

// Cylinder
// p : point to test
// c : center of the base of the cylinder
// r : radius of the cylinder
// h : height of the cylinder
Val Cylinder(vec3 p, vec3 c, float r, float h) {
    vec2 d = abs(vec2(length(p.xz - c.xz), p.y - c.y)) - vec2(r, h);
    return Val(min(max(d.x, d.y), 0.0) + length(max(d, 0.0)), 7); // => cout à 7
}

// Capsule
// p : point to test
// a : starting point of the capsule's line segment
// b : ending point of the capsule's line segment
// r : radius of the capsule
Val Capsule(vec3 p, vec3 a, vec3 b, float r) {
    vec3 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return Val(length(pa - ba * h) - r, 9); // => cout à 9
}

// Torus
// p : point to test
// c : center of the torus
// t : 2D vector representing the major (t.x) and minor (t.y) radii of the torus
Val Torus(vec3 p, vec3 c, vec2 t) {
    vec3 q = p - c;
    return Val(length(vec2(length(q.xz) - t.x, q.y)) - t.y, 11); // => cout à 11
}

// BoundingBox
// p : point à tester
// minCorner : coin minimum de la boîte englobante
// maxCorner : coin maximum de la boîte englobante
Val BoundingBox(vec3 p, vec3 minCorner, vec3 maxCorner) {
    vec3 d = max(minCorner - p, p - maxCorner);
    float dist = length(max(d, 0.0));
    return Val(dist, 1); // => cout à 1
}


// Transformed primitives ------------------------------------------------------------

Val object(vec3 p) {
    // la boite englobante pour le camion
    Val boundingBoxTruck = BoundingBox(p, vec3(-12.0, -6.5, -7.0), vec3(7.0, 6.5, 12.0));


    // si le rayon ne touche pas la boite, on passe
    if (boundingBoxTruck.v > 0.0) {
        return boundingBoxTruck;
    }

    // si le rayon touche la boite on teste les objets à l interieur
    vec3 translationDown = vec3(0.0, 0.0, -5.0);
    vec3 translatedP = Translate(p, translationDown);

    // coprs du camion
    Val remorque = Box(translatedP, vec3(-1.5, 0.0, 0.75) * 2.0, vec3(3.0, 1.2, 0.75) * 2.0);
    Val cabine = Box(translatedP, vec3(2.0, 0.0, 1.0) * 2.0, vec3(1.0, 1.0, 1.0) * 2.0);

    // roues
    Val roueAvantGauche = Cylinder(translatedP, vec3(1.5, -1.0, 0.3) * 2.0, 0.5 * 2.0, 0.25 * 2.0);
    Val roueAvantDroite = Cylinder(translatedP, vec3(1.5, 1.0, 0.3) * 2.0, 0.5 * 2.0, 0.25 * 2.0);
    Val roueArriereGauche = Cylinder(translatedP, vec3(-2.0, -1.0, 0.3) * 2.0, 0.5 * 2.0, 0.25 * 2.0);
    Val roueArriereDroite = Cylinder(translatedP, vec3(-2.0, 1.0, 0.3) * 2.0, 0.5 * 2.0, 0.25 * 2.0);

    // ce qui a en dessus du camion
    Val spokeRotated90 = Cylinder(RotateZ(translatedP, radians(90.0)), vec3(0.0, -3.0, 10.0), 6.0, 0.5);
    Val ring = Torus(translatedP, vec3(-3.5, 0.0, 10.0), vec2(5.5, 0.8));

    // construction du camion complet
    Val camion = Union(remorque,
                       Union(cabine,
                             Union(roueAvantGauche,
                                   Union(roueAvantDroite,
                                         Union(roueArriereGauche,
                                               Union(roueArriereDroite,
                                                   Union(spokeRotated90, ring)))))));

    return camion;
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
  return mix(vec3(.45,.55,.99),vec3(.65,.69,.99),ray.d.z*.5+.5);
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
vec3 Shade(vec3 p,vec3 n,Ray eye)
{
  // Point light
  const vec3 lp=vec3(5.,10.,25.);

  // Light direction to point light
  vec3 l=normalize(lp-p);

  // Ambient color
  vec3 ambient=.25+.25*background(Ray(p,n));

  // Shadow computation
  float shadow=Shadow(p,n,l);

  // Phong diffuse
  vec3 diffuse=.35*clamp(dot(n,l),0.,1.)*vec3(1.,1.,1.);

  // Specular
  vec3 r=reflect(eye.d,n);
  vec3 specular=.15*pow(clamp(dot(r,l),0.,1.),35.)*vec3(1.,1.,1.);
  vec3 c=ambient+shadow*(diffuse+specular);
  return c;
}

// Shading according to the number of steps in sphere tracing
// n : Number of steps
vec3 ShadeSteps(int n,int m)
{
  float t=float(n)/(float(m));
  return.5+mix(vec3(.05,.05,.5),vec3(.65,.39,.65),t);
}


// Main image function
void mainImage(out vec4 color, in vec2 pxy) {
    // Convert pixel coordinates
    vec2 pixel = (-iResolution.xy + 2.0 * pxy) / iResolution.y;

    // Mouse
    vec2 m = iMouse.xy / iResolution.xy;

    // Camera
    Ray ray = CreateRay(m, pixel);

    // Trace ray
    float t = 0.0;
    int s = 0;
    int c = 0;  // Cost
    bool hit = SphereTrace(ray, 100.0, t, s, c); // cout calculé ici

    // Background color
    vec3 rgb = background(ray);

    if (hit) {
        vec3 p = Point(ray, t);
        vec3 n = ObjectNormal(p);
        rgb = Shade(p, n, ray);
    }

    // Uncomment this line to shade cost
    rgb=ShadeSteps(c,500);


    color = vec4(rgb, 1.0);
}
